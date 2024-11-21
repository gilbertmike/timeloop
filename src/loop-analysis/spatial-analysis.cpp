#include "loop-analysis/spatial-analysis.hpp"

#include <boost/range/adaptor/reversed.hpp>
#include <isl/cpp.h>
#include <isl/map.h>
#include <isl/polynomial_type.h>
#include <barvinok/isl.h>

#include "isl-wrapper/ctx-manager.hpp"
#include "isl-wrapper/isl-functions.hpp"

namespace analysis
{

/******************************************************************************
 * Local declarations
 *****************************************************************************/

isl::map MakeMeshConnectivity(size_t num_spatial_dims);

std::vector<size_t> GetSpatialTagsIdxs(const std::vector<SpaceTime>& tags,
                                       BufferId buf_id);

std::optional<int> GetLastTemporalTagIdx(const std::vector<SpaceTime>& tags);

std::vector<size_t> MakeConnectivityDimPermutation(
  const std::vector<size_t>& spatial_idxs,
  size_t n_dims
);

std::vector<bool> MakeMulticastDimRemoveMask(
  const std::vector<SpaceTime>& tags,
  BufferId buf_id
);

/******************************************************************************
 * Global function implementations
 *****************************************************************************/

FillProvider::FillProvider(const LogicalBuffer& buf, const Occupancy& occ) : 
  buf(buf), occupancy(occ)
{
}


SpatialReuseAnalysisInput::SpatialReuseAnalysisInput(
  LogicalBuffer buf,
  const Fill& children_fill,
  bool count_hops
) : buf(buf), children_fill(children_fill), count_hops(count_hops)
{
}

SpatialReuseInfo SpatialReuseAnalysis(const SpatialReuseAnalysisInput& input)
{
  auto transfer_infos = std::vector<TransferInfo>();

  auto fill_to_be_provided = input.children_fill;
  for (const auto& fill_provider : input.fill_providers)
  {
    transfer_infos.emplace_back(
      fill_provider.spatial_reuse_model->Apply(input.buf.buffer_id,
                                               fill_to_be_provided,
                                               fill_provider.occupancy)
    );
    fill_to_be_provided = transfer_infos.back().unfulfilled_fill;
  }

  return SpatialReuseInfo{.transfer_infos = std::move(transfer_infos)};
}


SimpleLinkTransferModel::SimpleLinkTransferModel()
{
}

TransferInfo SimpleLinkTransferModel::Apply(
  BufferId buf_id,
  const Fill& fill,
  const Occupancy& occupancy
) const
{
  if (fill.dim_in_tags.size() != occupancy.dim_in_tags.size())
  {
    throw std::logic_error("fill and occupancy have different sizes");
  }

  auto n = fill.dim_in_tags.size();

  auto spatial_dim_idxs = GetSpatialTagsIdxs(fill.dim_in_tags, buf_id);
  auto n_spatial_dims = spatial_dim_idxs.size();

  auto last_temporal_opt = GetLastTemporalTagIdx(fill.dim_in_tags);

  auto transfer_info = TransferInfo();
  transfer_info.is_link_transfer = true;

  if (!last_temporal_opt || n_spatial_dims == 0) // Cannot fulfill via link transfers
  {
    transfer_info.fulfilled_fill =
      Transfers(fill.dim_in_tags, fill.map.subtract(fill.map)); // empty map
    transfer_info.parent_reads = 
      Reads(occupancy.dim_in_tags,
            occupancy.map.subtract(occupancy.map)); // empty map
    transfer_info.unfulfilled_fill = fill;
    auto domain = fill.map.domain();
    auto p_domain_space = domain.space().release();
    auto p_hops = isl_pw_qpolynomial_from_qpolynomial(
      isl_qpolynomial_zero_on_domain(p_domain_space)
    );
    transfer_info.p_hops = p_hops;

    return transfer_info;
  }

  auto connectivity = MakeMeshConnectivity(n_spatial_dims);
  auto padded_connectivity =
    isl::insert_equal_dims(connectivity, 0, 0, n - n_spatial_dims - 1);
  auto permutation = MakeConnectivityDimPermutation(spatial_dim_idxs, n);
  auto p_reorder_map = isl::reorder_projector(GetIslCtx().get(), permutation);
  auto complete_connectivity = isl::manage(isl_map_apply_range(
    isl_map_apply_range(
      isl_map_copy(p_reorder_map),
      padded_connectivity.release()
    ),
    isl_map_reverse(p_reorder_map)
  ));

  auto available_from_neighbors =
    complete_connectivity.apply_range(occupancy.map);
  auto fill_set = fill.map.intersect(available_from_neighbors);
  auto remaining_fill = fill.map.subtract(fill_set);

  auto p_fill_set = isl_map_wrap(fill_set.copy());
  auto p_hops = isl_pw_qpolynomial_from_qpolynomial(
    isl_qpolynomial_one_on_domain(isl_set_get_space(p_fill_set))
  );
  p_hops = isl_pw_qpolynomial_intersect_domain(p_hops, p_fill_set);

  transfer_info.fulfilled_fill = Transfers(fill.dim_in_tags, fill_set);
  transfer_info.parent_reads= Reads(
    fill.dim_in_tags,
    fill_set.subtract(fill_set) // Empty set since there are no parent reads
  );
  transfer_info.unfulfilled_fill = Fill(fill.dim_in_tags, remaining_fill);
  transfer_info.p_hops = p_hops;

  return transfer_info;
}


struct HopsAccesses
{
  double hops;
  double accesses;

  /**
   * @brief Simple weighted average for hops and accumulation for accesses.
   */
  void InsertHopsAccesses(double extra_hops, double extra_accesses)
  {
    accesses += extra_accesses;
    hops += extra_hops;
  }
};

struct Accumulator
{
  std::map<uint64_t, HopsAccesses> multicast_to_hops_accesses;
  isl_pw_qpolynomial* p_time_data_to_hops;
};

/**
 * @brief Accumulates scatter, hops, and accesses for many multicast factors.
 *
 * @param p_domain A set with signature $\{ [st_{n-1},t_n] -> data \}$ where
 *                 data is some set of data.
 * @param p_multicast_factor A qpolynomial assumed to be constant that equals
 *                           the multicast factor.
 * @param p_voided_accumulator Voided pointer that is cast into Accumulator.
 */
isl_stat ComputeMulticastScatterHops(isl_set* p_domain,
                                     isl_qpolynomial* p_multicast_factor,
                                     void* p_voided_accumulator)
{
  auto& accumulator = *static_cast<Accumulator*>(p_voided_accumulator);
  // WARNING: assumes constant multicast factor over piecewise domain.
  // It is unclear what conditions may cause this to break.
  auto multicast_factor = isl::val_to_double(
    isl_qpolynomial_eval(p_multicast_factor,
                         isl_set_sample_point(isl_set_copy(p_domain)))
  );
  auto& hops_accesses = accumulator.multicast_to_hops_accesses[multicast_factor];

  auto p_hops_pw_qp = isl_set_apply_pw_qpolynomial(
    isl_set_copy(p_domain),
    isl_pw_qpolynomial_copy(accumulator.p_time_data_to_hops)
  );
  if (isl_pw_qpolynomial_isa_qpolynomial(p_hops_pw_qp) == isl_bool_false)
  {
    throw std::runtime_error("accesses is not a single qpolynomial");
  }
  auto hops = isl::val_to_double(isl_qpolynomial_get_constant_val(
    isl_pw_qpolynomial_as_qpolynomial(p_hops_pw_qp)
  ));

  auto p_time_to_data = isl_set_unwrap(p_domain);
  auto p_accesses_pw_qp = isl_pw_qpolynomial_sum(isl_map_card(p_time_to_data));
  if (isl_pw_qpolynomial_isa_qpolynomial(p_accesses_pw_qp) == isl_bool_false)
  {
    throw std::runtime_error("accesses is not a single qpolynomial");
  }
  auto accesses = isl::val_to_double(isl_qpolynomial_get_constant_val(
    isl_pw_qpolynomial_as_qpolynomial(p_accesses_pw_qp)
  ));

  hops_accesses.InsertHopsAccesses(hops, accesses);

  return isl_stat_ok;
}

SimpleMulticastModel::SimpleMulticastModel(bool count_hops) :
  count_hops_(count_hops)
{
}

TransferInfo
SimpleMulticastModel::Apply(
  BufferId buf_id,
  const Fill& fill,
  const Occupancy& occ
) const
{
  (void) occ;
  auto transfer_info = TransferInfo();
  transfer_info.is_multicast = true;

  auto n = isl::dim(fill.map, isl_dim_in);

  auto spatial_dim_idxs = GetSpatialTagsIdxs(fill.dim_in_tags, buf_id);
  auto n_spatial_dims = spatial_dim_idxs.size();
  auto permutation = MakeConnectivityDimPermutation(spatial_dim_idxs, n);
  auto p_reorder_map = isl::reorder_projector(GetIslCtx().get(), permutation);

  auto p_fill = fill.map.copy();
  p_fill = isl_map_apply_range(isl_map_reverse(p_reorder_map), p_fill);

  // Creates [[t] -> [data]] -> [t, x, y] when n_spatial_dims == 2
  // Creates [[t] -> [data]] -> [t, x] when n_spatial_dims == 1
  auto p_wrapped_fill = isl_map_uncurry(isl_map_project_out(
    isl_map_reverse(isl_map_range_map(isl_map_reverse(p_fill))),
    isl_dim_in,
    n-n_spatial_dims,
    n_spatial_dims
  ));
  auto wrapped_fill = isl::manage(p_wrapped_fill);

  auto p_multicast_factor = isl_map_card(wrapped_fill.copy());

  isl_pw_qpolynomial* p_hops = nullptr;
  if (count_hops_)
  {
    isl_map* p_data_to_x = nullptr;
    if (n_spatial_dims == 2)
    {
      auto p_y_hops_cost = isl_pw_qpolynomial_from_qpolynomial(
        isl_qpolynomial_add(
          isl_qpolynomial_var_on_domain(
            isl_space_range(isl_map_get_space(wrapped_fill.get())),
            isl_dim_set,
            n-1
          ),
          isl_qpolynomial_one_on_domain(
            isl_space_range(isl_map_get_space(wrapped_fill.get()))
          )
        )
      );
      auto p_y_hops = isl_map_apply_pw_qpolynomial(wrapped_fill.copy(),
                                                    p_y_hops_cost);
      p_hops = p_y_hops;
      p_data_to_x = isl_map_lexmax(
        isl_map_project_out(wrapped_fill.copy(), isl_dim_out, n-1, 1)
      );
    }
    else
    {
      p_data_to_x = wrapped_fill.copy();
    }

    // Remove y, leaving only x
    auto p_x_hops_cost = isl_pw_qpolynomial_from_qpolynomial(
      isl_qpolynomial_add(
        isl_qpolynomial_var_on_domain(
          isl_space_range(isl_map_get_space(p_data_to_x)),
          isl_dim_set,
          n-n_spatial_dims
        ),
        isl_qpolynomial_one_on_domain(
          isl_space_range(isl_map_get_space(p_data_to_x))
        )
      )
    );
    auto p_x_hops = isl_map_apply_pw_qpolynomial(p_data_to_x,
                                                 p_x_hops_cost);

    if (p_hops == nullptr)
    {
      p_hops = p_x_hops;
    }
    else
    {
      p_hops = isl_pw_qpolynomial_add(p_hops, p_x_hops);
    }
  }
  else
  {
    p_hops = isl_pw_qpolynomial_from_qpolynomial(
      isl_qpolynomial_zero_on_domain(isl_pw_qpolynomial_get_domain_space(
        p_multicast_factor
      ))
    );
  }

  // auto accumulator = Accumulator();
  // accumulator.p_time_data_to_hops = p_hops;
  // isl_pw_qpolynomial_foreach_piece(p_multicast_factor,
  //                                  &ComputeMulticastScatterHops,
  //                                  static_cast<void*>(&accumulator));

  // isl_pw_qpolynomial_free(p_multicast_factor);

  // for (const auto& [multicast, hops_accesses] :
  //       accumulator.multicast_to_hops_accesses)
  // {
  //   auto& stats =
  //     transfer_info.compat_access_stats[std::make_pair(multicast, 1)];
  //   stats.accesses = hops_accesses.accesses;
  //   stats.hops = hops_accesses.hops / hops_accesses.accesses;
  // }
  // Remove y, leaving only x

  auto total_accesses = isl::val_to_double(isl_qpolynomial_get_constant_val(
    isl_pw_qpolynomial_as_qpolynomial(
      isl_set_card(isl_pw_qpolynomial_domain(
        isl_pw_qpolynomial_copy(p_multicast_factor)
      ))
    )
  ));

  auto p_total_hops = isl_pw_qpolynomial_sum(isl_pw_qpolynomial_sum(
    isl_pw_qpolynomial_copy(p_hops)
  ));
  auto total_hops = isl::val_to_double(isl_qpolynomial_get_constant_val(
    isl_pw_qpolynomial_as_qpolynomial(p_total_hops)
  ));

  auto p_total_multicast = isl_pw_qpolynomial_sum(isl_pw_qpolynomial_sum(
    p_multicast_factor
  ));
  auto total_multicast = isl::val_to_double(isl_qpolynomial_get_constant_val(
    isl_pw_qpolynomial_as_qpolynomial(p_total_multicast)
  ));

  auto avg_multicast = total_multicast / total_accesses;
  auto avg_hops = total_hops / total_accesses;

  auto& stats = transfer_info.compat_access_stats[std::make_pair(
    avg_multicast,
    1
  )];
  stats.accesses = total_accesses;
  stats.hops = avg_hops;

  transfer_info.fulfilled_fill = Transfers(fill.dim_in_tags, fill.map);

  auto project_out_mask = MakeMulticastDimRemoveMask(fill.dim_in_tags, buf_id);
  auto domain = fill.map.domain();
  auto p_space = isl_set_get_space(domain.release());
  auto p_projector = isl::dim_projector(p_space, project_out_mask);
  auto parent_reads_tags = std::vector<SpaceTime>();
  size_t i = 0;
  for (const auto should_remove : project_out_mask)
  {
    if (!should_remove)
    {
      parent_reads_tags.emplace_back(fill.dim_in_tags.at(i));
    }
    ++i;
  }
  transfer_info.parent_reads = Reads(
    parent_reads_tags,
    isl::manage(isl_map_apply_range(p_projector, fill.map.copy()))
  );
  transfer_info.unfulfilled_fill = Fill(
    fill.dim_in_tags,
    fill.map.subtract(fill.map)  // empty map
  );
  transfer_info.p_hops = p_hops;

  return transfer_info;
}


DistributedMulticastHypercubeModel::DistributedMulticastHypercubeModel(
  bool count_hops, isl::map dist_func
): count_hops_(count_hops), dist_func_(dist_func)
{
}


DistributedMulticastOrderedExtentsDORModel::DistributedMulticastOrderedExtentsDORModel(
  bool count_hops, isl::map dist_func
): count_hops_(count_hops), dist_func_(dist_func)
{
}


/******************************************************************************
 * Local function implementations
 *****************************************************************************/


isl::map MakeMeshConnectivity(size_t num_spatial_dims)
{
  if (num_spatial_dims == 2)
  {
    return isl::map(GetIslCtx(),
                    "{ [t, x, y] -> [t-1, x', y'] : "
                    " (y'=y and x'=x-1) or (y'=y and x'=x+1) "
                    " or (x'=x and y'=y-1) or (x'=x and y'=y+1) }");
  }
  else if (num_spatial_dims == 1)
  {
    return isl::map(GetIslCtx(),
                    "{ [t, x] -> [t-1, x'] : "
                    " (x'=x-1) or (x'=x+1) }");
  }

  auto err_msg = std::stringstream();
  err_msg << "Cannot make mesh with " << num_spatial_dims << "dims";
  throw std::logic_error(err_msg.str());
}

std::vector<size_t> GetSpatialTagsIdxs(const std::vector<SpaceTime>& tags,
                                    BufferId buf_id)
{
  std::vector<size_t> spatial_dim_idxs;
  int i = 0;
  for (const auto& tag : tags)
  {
    if (std::holds_alternative<Spatial>(tag))
    {
      const auto& spatial_tag = std::get<Spatial>(tag);
      if (spatial_tag.target == buf_id)
      {
        spatial_dim_idxs.emplace_back(i);
      }
    }
    ++i;
  }
  return spatial_dim_idxs;
}

std::optional<int> GetLastTemporalTagIdx(const std::vector<SpaceTime>& tags)
{
  if (tags.size() == 0)
  {
    return std::nullopt;
  }

  size_t idx = tags.size()-1;
  for (const auto& tag : tags | boost::adaptors::reversed)
  {
    if (analysis::IsTemporal(tag))
    {
      return idx;
    }
    --idx;
  }

  return std::nullopt;
}

std::vector<size_t> MakeConnectivityDimPermutation(
  const std::vector<size_t>& spatial_idxs,
  size_t n_dims
)
{
  auto permutation = std::vector<size_t>();

  size_t cur_spatial_i = 0;
  for (size_t i = 0; i < n_dims; ++i)
  {
    if (cur_spatial_i < spatial_idxs.size()
        && i == spatial_idxs.at(cur_spatial_i))
    {
      ++cur_spatial_i;
    }
    else
    {
      permutation.emplace_back(i);
    }
  }

  for (const auto spatial_idx : spatial_idxs)
  {
    permutation.emplace_back(spatial_idx);
  }

  return permutation;
}

std::vector<bool> MakeMulticastDimRemoveMask(
  const std::vector<SpaceTime>& tags,
  BufferId buf_id
)
{
  std::vector<bool> mask;
  for (const auto& tag : tags)
  {
    if (std::holds_alternative<Spatial>(tag))
    {
      const auto& spatial_tag = std::get<Spatial>(tag);
      mask.emplace_back(spatial_tag.target == buf_id);
    }
    else
    {
      mask.emplace_back(false);
    }
  }
  return mask;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Multicast Model Pit of Temporary Spaghetti
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// NOTES FOR NON-TREE MULTICAST SCENARIO
// - Load balancing issues for multiple minimally distant sources.
// - Compose minimal distances with the other set to remove non-minimal pairs then
// move on with the rest of the algorithm.
__isl_give const isl::map identify_mesh_casts( 
    __isl_take const isl::map src_occupancy, 
    __isl_take const isl::map dst_fill, 
    __isl_take const isl::map dist_func
) {
  /* Makes [[dst -> data] -> dst] -> [data] */
  isl::set wrapped_dst_fill = dst_fill.wrap();
  isl::map wrapped_fill_identity = isl::manage(isl_map_identity(
    isl_space_copy(wrapped_dst_fill.get_space().map_from_set().get())
  ));
  wrapped_fill_identity = wrapped_fill_identity.intersect_domain(wrapped_dst_fill);

  /* Makes [dst -> data] -> [dst -> data] */
  isl::map uncurried_fill_identity = wrapped_fill_identity.uncurry();

  /* Inverts src_occupancy such that data implies source.
  * i.e. {src -> data} becomes {data -> src} */
  isl::map data_presence = src_occupancy.reverse();

  isl::map fills_to_dst_TO_src = uncurried_fill_identity.apply_range(
    data_presence
  );
  isl::map fills_to_matches = fills_to_dst_TO_src.curry();

  // Calculates the distance of all the dst-src pairs with matching data.
  isl::map distances_map = fills_to_matches.apply_range(dist_func);
  isl::map fills_to_matches_TO_matches = fills_to_matches.range_map().as_map();
  isl::map fills_to_matches_TO_dist = fills_to_matches_TO_matches.apply_range(dist_func);

  // Gets the minimal distance pairs.
  isl::map lexmin_distances = distances_map.lexmin();
  isl::map assoc_dist_with_src = lexmin_distances.apply_range(
    fills_to_matches_TO_dist.reverse()
  );
  // Isolates the relevant minimal pairs.
  isl::map minimal_pairs = assoc_dist_with_src.range().unwrap();
  // Isolates the multicast networks.
  isl::map multicast_networks = minimal_pairs.curry().range().unwrap().uncurry().lexmin().curry();

  return multicast_networks;
}

////////////////////////////////////////////////////////////////////////////////
// Extent Based Models
////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Gets the extent of each dimension per [data -> src].
 * @param mesh_cast_networks The map of [data -> src] -> dst.
 * @param dist_func The map of [data -> src] -> distance.
 * 
 * @return A vector of piecewise affines mapping [data -> src] -> extent.
 */
std::vector<__isl_give isl::pw_aff> calculate_extents(
  __isl_take isl::map mesh_cast_networks,
  __isl_take isl::map dist_func
) {
  /// @brief Makes mesh_cash_networks from data -> [dst -> src] to [data -> src] -> dst
  isl::map potential_sources = mesh_cast_networks.range_reverse().uncurry();
  /// @note Sources are part of the extents, so we union it with the destinations.
  isl::map sources = potential_sources.domain().unwrap().range_map().as_map();
  isl::map casting_extents = isl::manage(
    isl_map_union(sources.copy(), potential_sources.copy())
  );

  /// @brief Projects away all dimensions but one to find their extent for hypercube.
  unsigned dimensions = potential_sources.range_tuple_dim();
  long min_cost = LONG_MAX;
  // Creates a mask of what to project out.
  std::vector<bool> project_out = std::vector<bool>(dimensions, true);
  std::vector<isl::pw_aff> dim_extents = std::vector<isl::pw_aff>(dimensions);

  /// @brief Gets the extents of all the dimensions.
  for (unsigned noc_dim = 0; noc_dim < dimensions; noc_dim++) {
    /// @brief Projects out all the dimensions of the output besides noc_dim.
    project_out[noc_dim] = false;
    isl::map extent_mapper = isl::manage(isl::dim_projector<std::vector<bool>>(
      isl_space_copy(casting_extents.range().get_space().get()), project_out
    )).reverse();
    isl::map dim_extent_space = casting_extents.apply_range(extent_mapper);
    project_out[noc_dim] = true;

    // Finds max(noc_dim) - min(noc_dim) for each [dara -> src]
    isl::pw_aff max_extent = isl::manage(isl_map_dim_max(dim_extent_space.copy(), 0));
    isl::pw_aff min_extent = isl::manage(isl_map_dim_min(dim_extent_space.copy(), 0));

    // Subtracts the max from the min to get the extent per [data -> src]
    dim_extents[noc_dim] = max_extent.sub(min_extent).coalesce();
  }

  return dim_extents;
}

/**
 * @brief Calculates the cost of casting from a hypercube.
 * @param mesh_cast_networks The map of [data -> src] -> dst.
 * @param dist_func The map of [data -> src] -> distance.
 * 
 * @return A piecewise qpolynomial representing the cost of hypercube casting.
 */
isl_pw_qpolynomial *cost_mesh_cast_hypercube(
    __isl_take const isl::map mesh_cast_networks,
    __isl_take const isl::map dist_func
) { 
  auto dim_extents = calculate_extents(mesh_cast_networks, dist_func);

  // Tracks the total cost of the hypercube cast per src -> data.
  isl::pw_aff one = isl::manage(isl_pw_aff_val_on_domain(
    isl_pw_aff_domain(dim_extents[0].copy()),
    isl_val_int_from_si(GetIslCtx().get(), 1)
  ));
  isl_pw_qpolynomial *hypercube_costs = isl_pw_qpolynomial_from_pw_aff(one.copy());
  // 
  /**
   * @brief Calculates the cost of the hypercube.
   * @note Hypercube cost = 
   *       = \sum_{i=0}^{D} ((extent_i - 1) * \prod_{j=0}^{i-1} extent_j)
   *       = \prod_{i=0}^{D} extent_i
   */
  for (auto& dim_extent : dim_extents) {
    // Adds the dim_extent times the casting volume to the hypercube cost.
    auto dim_plus = isl_pw_qpolynomial_from_pw_aff(
      dim_extent.add(one).coalesce().release()
    );
    hypercube_costs = isl_pw_qpolynomial_coalesce(
      isl_pw_qpolynomial_mul(hypercube_costs, dim_plus)
    );
  }
  hypercube_costs = isl_pw_qpolynomial_coalesce(isl_pw_qpolynomial_sub(
    hypercube_costs, isl_pw_qpolynomial_from_pw_aff(one.release())
  ));

  // returns the hypercube cost as a piecewise polynomial.
  return isl_pw_qpolynomial_sum(isl_pw_qpolynomial_sum(hypercube_costs));
}


TransferInfo DistributedMulticastHypercubeModel::Apply(
  BufferId buf_id,
  const Fill& fills,
  const Occupancy& occupancy
) const
{
  (void) buf_id;

  isl::map mcs = identify_mesh_casts(
    occupancy.map, fills.map, this->dist_func_
  );
  isl_pw_qpolynomial *res = cost_mesh_cast_hypercube(
    mcs, this->dist_func_
  );

  // TODO:: Read once from all buffers, assert that card(mcs) == tensor_size * D
  return TransferInfo{
    .fulfilled_fill=Transfers(fills.dim_in_tags, fills.map),
    .parent_reads=Reads(occupancy.dim_in_tags, mcs),
    .unfulfilled_fill=Fill(fills.dim_in_tags, fills.map.subtract(fills.map)),
    .p_hops=res,
  };
}



isl_pw_qpolynomial *cost_mesh_cast_extent_first(
    __isl_take const isl::map mesh_cast_networks,
    __isl_take const isl::map dist_func
) {
  // Gets the [data -> src] -> extent of each dimension.
  const std::vector<isl::pw_aff> dim_extents = calculate_extents(
    mesh_cast_networks, dist_func
  );
  /// @brief Argsorts the dimensions by their extents.
  std::vector<size_t> sorted_dims = std::vector<size_t>(dim_extents.size());
  std::iota(sorted_dims.begin(), sorted_dims.end(), 0);
  std::sort(sorted_dims.begin(), sorted_dims.end(),
    ///@brief Lambda function that reverse argsorts by total extent.
    [&dim_extents](size_t a, size_t b) -> bool
    {
      // Calculates the sum of the extents for each [data -> src].
      isl_pw_qpolynomial *a_total = isl_pw_qpolynomial_sum(
        isl_pw_qpolynomial_from_pw_aff(dim_extents[a].copy())
      );
      isl_pw_qpolynomial *b_total = isl_pw_qpolynomial_sum(
        isl_pw_qpolynomial_from_pw_aff(dim_extents[b].copy())
      );
      // Gets the constant value of the sum.
      long a_val = isl::manage(
        isl_pw_qpolynomial_eval(a_total, 
          isl_point_zero(
            isl_pw_qpolynomial_get_domain_space(a_total)
          )
        )
      ).get_num_si();
      long b_val = isl::manage(
        isl_pw_qpolynomial_eval(b_total, 
          isl_point_zero(
            isl_pw_qpolynomial_get_domain_space(b_total)
          )
        )
      ).get_num_si();
      return a_val > b_val;
    }
  );

  // Casts in descending order of extents.
  isl::val one = isl::manage(isl_val_int_from_si(GetIslCtx().get(), 1));
  
  return nullptr;
}


TransferInfo DistributedMulticastOrderedExtentsDORModel::Apply(
  BufferId buf_id,
  const Fill& fills,
  const Occupancy& occupancy
) const
{
  (void) buf_id;

  // Defines the distance function string.
  isl::map mcs = identify_mesh_casts(
    occupancy.map, 
    fills.map, 
    this->dist_func_
  );

  // Calculates the cost of the extent-first model.
  isl_pw_qpolynomial *res = cost_mesh_cast_extent_first(mcs, this->dist_func_);

  // TODO:: Read once from all buffers, assert that card(mcs) == tensor_size * D
  return TransferInfo{
    .fulfilled_fill=Transfers(fills.dim_in_tags, fills.map),
    .parent_reads=Reads(occupancy.dim_in_tags, mcs),
    .unfulfilled_fill=Fill(fills.dim_in_tags, fills.map.subtract(fills.map)),
    .p_hops=res,
  };
}
} // namespace analysis