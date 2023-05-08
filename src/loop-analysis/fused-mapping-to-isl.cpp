#include "loop-analysis/isl-ir.hpp"
#include "loop-analysis/temporal-analysis.hpp"
#include "isl-wrapper/ctx-manager.hpp"

namespace analysis
{

/******************************************************************************
 * Local declarations
 *****************************************************************************/

LogicalBufTiling LogicalBufTilingFromMapping(mapping::FusedMapping& mapping,
                                             BranchTilings branch_tiling);

std::map<DataSpaceID, size_t>
DspaceTopIdxFromMapping(mapping::FusedMapping& mapping);

std::optional<size_t> BranchIdxFromMapping(mapping::FusedMapping& mapping);

BranchTilings
LoopBoundsInference(BranchTilings tilings,
                    mapping::FusedMapping& mapping,
                    const problem::FusedWorkload& workload,
                    size_t pipeline_tiling_idx,
                    const std::map<DataSpaceID, size_t>& dspace_top_idx);

std::vector<std::pair<LogicalBuffer, size_t>>
BufferIterLevelsFromMapping(mapping::FusedMapping& mapping);

BranchTilings TilingFromMapping(mapping::FusedMapping& mapping,
                                const problem::FusedWorkload& workload);

std::map<problem::DimensionId, std::map<problem::DimensionId, int>>
EinsumDimensionStridesFromWorkload(const problem::FusedWorkload& workload);

void GatherSubsequentEinsumStrides(
  std::map<problem::DimensionId, std::map<problem::DimensionId, int>>& strides,
  problem::EinsumId prod_einsum,
  problem::EinsumId cons_einsum,
  problem::DataSpaceId dspace,
  const problem::FusedWorkload& workload
);

/******************************************************************************
 * Global function implementations
 *****************************************************************************/
LogicalBufOccupancies
OccupanciesFromMapping(mapping::FusedMapping& mapping,
                       const problem::FusedWorkload& workload)
{
  auto branch_tiling = analysis::TilingFromMapping(mapping, workload);
  auto branch_idx = analysis::BranchIdxFromMapping(mapping);
  auto dspace_indices = analysis::DspaceTopIdxFromMapping(mapping);
  branch_tiling = analysis::LoopBoundsInference(std::move(branch_tiling),
                                                mapping,
                                                workload,
                                                *branch_idx,
                                                dspace_indices);

  LogicalBufOccupancies occupancies;
  auto tilings = LogicalBufTilingFromMapping(mapping, branch_tiling);
  for (const auto& [buf, tiling] : tilings)
  {
    std::cout << buf << ": " << tiling << std::endl;
  }
  // auto buf_skew = LogicalBufSkewsFromMapping(maping);

  return occupancies;
}

/******************************************************************************
 * Local function implementations
 *****************************************************************************/

LogicalBufTiling LogicalBufTilingFromMapping(mapping::FusedMapping& mapping,
                                             BranchTilings branch_tiling)
{
  auto buf_to_iter_level = BufferIterLevelsFromMapping(mapping);

  LogicalBufTiling result;
  for (auto& [buf, level] : buf_to_iter_level)
  {
    result.emplace(std::make_pair(
      buf,
      project_dim_in_after(isl::map(branch_tiling.at(buf.branch_leaf_id)),
                           level)
    ));
  }

  return result;
}

BranchTilings LoopBoundsInference(BranchTilings tilings,
                                  mapping::FusedMapping& mapping,
                                  const problem::FusedWorkload& workload,
                                  size_t pipeline_tiling_idx,
                                  const std::map<DataSpaceID, size_t>& dspace_top_idx)
{
  BranchTilings inferred_tilings(tilings.begin(), tilings.end());

  for (size_t i = 0; i < inferred_tilings.size(); ++i)
  {
    for (auto& [leaf_id, tiling] : inferred_tilings)
    {
      bool complete = true;
      auto domain = tiling.domain();
      for (auto i = 0; i < isl_set_dim(domain.get(), isl_dim_set); ++i)
      {
        complete = complete &&
                   isl_set_dim_has_lower_bound(domain.get(), isl_dim_set, i);
      }
      if (!complete)
      {
        continue;
      }

      auto einsum_id =
        std::get<mapping::Compute>(mapping.NodeAt(leaf_id)).kernel;
      for (const auto& read_tensor : workload.TensorsReadByEinsum(einsum_id))
      {
        auto producer_einsum_opt = workload.WriterEinsum(read_tensor);
        if (!producer_einsum_opt)
        {
          // Not an intermediate tensor
          continue;
        }

        // Decide how much the consumer (einsum_id) needs
        auto pruned_tiling = project_dim_in_after(tiling, pipeline_tiling_idx);
        auto read_accesses = workload.ReadAccesses(einsum_id, read_tensor);
        auto required_data = pruned_tiling.apply_range(read_accesses);

        auto top_idx = dspace_top_idx.at(read_tensor);
        auto shifter = 
          isl::MapToPriorData(pipeline_tiling_idx, top_idx);
        auto buffered_data = shifter.apply_range(required_data);

        auto computed_data = required_data.subtract(buffered_data).coalesce();
        computed_data =
          computed_data.intersect_range(workload.DataSpaceBound(read_tensor));

        auto producer_write_dep =
          workload.WriteAccesses(*producer_einsum_opt, read_tensor);
        auto required_ops =
          computed_data.apply_range(producer_write_dep.reverse());

        for (const auto& [leaf_id, producer_tiling] : inferred_tilings)
        {
          if (std::get<mapping::Compute>(mapping.NodeAt(leaf_id)).kernel !=
              *producer_einsum_opt)
          {
            continue;
          }
          auto required_iters = ConstraintDimEquals(
            required_ops.apply_range(producer_tiling.reverse()),
            pipeline_tiling_idx
          );

          auto inferred_prod_tiling =
            producer_tiling
            .intersect_domain(required_iters.range())
            .coalesce();
          std::cout << inferred_prod_tiling << std::endl;
          inferred_tilings.at(leaf_id) = inferred_prod_tiling;
        }
      }
    }
  }
  return inferred_tilings;
}

std::optional<size_t> BranchIdxFromMapping(mapping::FusedMapping& mapping)
{
  auto idx = std::optional<size_t>();
  for (auto path : GetPaths(mapping))
  {
    size_t potential_idx = 0;
    for (const auto& node : path)
    {
      std::visit(
        [&idx, &potential_idx] (auto&& node) {
          using NodeT = std::decay_t<decltype(node)>;
          if constexpr (mapping::IsBranchV<NodeT>)
          {
            idx = potential_idx;
          }
          else if constexpr (mapping::IsLoopV<NodeT>)
          {
            ++potential_idx;
          }
        },
        node
      );
      if (idx)
      {
        break;
      }
    }
    if (idx)
    {
      break;
    }
  }

  return idx;
}

std::map<DataSpaceID, size_t>
DspaceTopIdxFromMapping(mapping::FusedMapping& mapping)
{
  std::map<DataSpaceID, size_t> dspace_to_idx;

  for (auto path : GetPaths(mapping))
  {
    size_t loop_idx = 0;
    for (const auto& node : path)
    {
      std::visit(
        [&dspace_to_idx, &loop_idx] (auto&& node) {
          using NodeT = std::decay_t<decltype(node)>;
          if constexpr (std::is_same_v<NodeT, mapping::Storage>)
          {
            auto dspace = node.dspace;
            if (dspace_to_idx.find(dspace) == dspace_to_idx.end())
            {
              dspace_to_idx[dspace] = loop_idx;
            }
          }
          else if constexpr (mapping::IsLoopV<NodeT>)
          {
            ++loop_idx;
          }
        },
        node
      );
    }
  }

  return dspace_to_idx;
}

std::vector<std::pair<LogicalBuffer, size_t>>
BufferIterLevelsFromMapping(mapping::FusedMapping& mapping)
{
  std::vector<std::pair<LogicalBuffer, size_t>> result;
  for (auto path : GetPaths(mapping))
  {
    size_t iter_idx = 0;
    std::vector<std::pair<LogicalBuffer, size_t>> new_results;
    for (const auto& node : path)
    {
      std::visit(
        [&new_results, &iter_idx] (auto&& node) {
          using NodeT = std::decay_t<decltype(node)>;

          if constexpr (std::is_same_v<NodeT, mapping::Storage>)
          {
            auto buffer = LogicalBuffer(node.buffer, node.dspace, 0);
            new_results.emplace_back(
              std::make_pair(std::move(buffer), iter_idx)
            );
          } else if constexpr (mapping::IsLoopV<NodeT>)
          {
            ++iter_idx;
          } else if constexpr (std::is_same_v<NodeT, mapping::Compute>)
          {
            for (auto& [buf, _] : new_results)
            {
              buf.branch_leaf_id = node.id;
            }
          }
        },
        node
      );
    }
    result.insert(result.end(), new_results.begin(), new_results.end());
  }

  return result;
}

BranchTilings TilingFromMapping(mapping::FusedMapping& mapping,
                                const problem::FusedWorkload& workload)
{
  BranchTilings result;
  for (auto path : GetPaths(mapping))
  {
    auto strides = EinsumDimensionStridesFromWorkload(workload);

    std::map<problem::DimensionId, std::vector<std::pair<size_t, int>>>
    prob_id_to_expr;

    size_t cur_dim_idx = 0;
    problem::EinsumId einsum_id;
    mapping::NodeID leaf_id;
    for (const auto& node : path)
    {
      std::visit(
        [&prob_id_to_expr, &cur_dim_idx, &einsum_id, &leaf_id] (auto&& node) {
          using NodeT = std::decay_t<decltype(node)>;
          if constexpr (std::is_same_v<NodeT, mapping::For>
                        || std::is_same_v<NodeT, mapping::ParFor>)
          {
            if (node.tile_size)
            {
              prob_id_to_expr[node.op_dim].emplace_back(std::make_pair(
                cur_dim_idx,
                *node.tile_size
              ));
            }
            ++cur_dim_idx;
          } else if constexpr (std::is_same_v<NodeT, mapping::Compute>)
          {
            leaf_id = node.id;
            einsum_id = node.kernel;
          }
        },
        node
      );
    }

    auto iter_space = isl::space_set_alloc(GetIslCtx(), 0, cur_dim_idx);
    auto eq_maff = isl::multi_aff::zero(isl::space_from_domain_and_range(
      isl::space_set_alloc(GetIslCtx(), 0, cur_dim_idx),
      workload.EinsumOspaceBound(einsum_id).space()
    ));
    auto iter_set = isl::set::universe(eq_maff.domain().space());
    auto zero_aff = isl::aff::zero_on_domain(iter_set.space());
    for (const auto& [prob_idx, expr] : prob_id_to_expr)
    {
      const auto& einsum_dim_to_idx = workload.EinsumDimToIdx(einsum_id);
      if (einsum_dim_to_idx.find(prob_idx) != einsum_dim_to_idx.end())
      {
        auto einsum_dim_idx = einsum_dim_to_idx.at(prob_idx);
        auto eq_aff = eq_maff.get_at(einsum_dim_idx);

        std::optional<decltype(expr)::value_type::second_type> last_coef;
        for (const auto& [iter_id, coef] : expr)
        {
          eq_aff = isl::set_coefficient_si(eq_aff, isl_dim_in, iter_id, coef);

          auto var_aff =
            isl::set_coefficient_si(zero_aff, isl_dim_in, iter_id, 1);
          iter_set = iter_set.intersect(var_aff.ge_set(zero_aff));

          if (last_coef)
          {
            auto upper_aff = isl::set_constant_si(zero_aff, *last_coef);
            iter_set = iter_set.intersect(var_aff.lt_set(upper_aff));
          }

          last_coef = coef;
        }

        eq_maff = eq_maff.set_at(einsum_dim_idx, eq_aff);

      }
      else
      {
        for (auto [einsum_dim, einsum_dim_idx] : einsum_dim_to_idx)
        {
          auto stride = strides[prob_idx][einsum_dim];
          auto eq_aff = eq_maff.get_at(einsum_dim_idx);

          std::optional<decltype(expr)::value_type::second_type> last_coef;
          for (const auto& [iter_id, coef] : expr)
          {
            eq_aff = isl::set_coefficient_si(eq_aff,
                                             isl_dim_in,
                                             iter_id,
                                             stride*coef);

            auto var_aff =
              isl::set_coefficient_si(zero_aff, isl_dim_in, iter_id, 1);
            iter_set = iter_set.intersect(var_aff.ge_set(zero_aff));

            if (last_coef)
            {
              auto upper_aff = isl::set_constant_si(zero_aff, *last_coef);
              iter_set = iter_set.intersect(var_aff.lt_set(upper_aff));
            }

            last_coef = coef;
          }

          eq_maff = eq_maff.set_at(einsum_dim_idx, eq_aff);
        }
      }
    }

    auto map = isl::map_from_multi_aff(eq_maff);
    map = map.intersect_domain(iter_set);

    result.emplace(std::make_pair(leaf_id, map));
  }

  for (const auto& [leaf_id, map] : result)
  {
    const auto& node_var = mapping.NodeAt(leaf_id);
    const auto& compute_node = std::get<mapping::Compute>(node_var);
    const auto& ospace_bound = workload.EinsumOspaceBound(compute_node.kernel);
    result[leaf_id] = map.intersect_range(ospace_bound);
  }

  return result;
}

std::map<problem::DimensionId, std::map<problem::DimensionId, int>>
EinsumDimensionStridesFromWorkload(const problem::FusedWorkload& workload)
{
  std::map<problem::DimensionId, std::map<problem::DimensionId, int>> strides;

  for (size_t i = 0; i < workload.EinsumNameToId().size(); ++i)
  {
    for (const auto& [_, cons_einsum] : workload.EinsumNameToId())
    {
      for (const auto& dspace : workload.TensorsReadByEinsum(cons_einsum))
      {
        auto prod_einsum_opt = workload.WriterEinsum(dspace);
        if (!prod_einsum_opt)
        {
          continue;
        }

        const auto& prod_einsum = *prod_einsum_opt;
        GatherSubsequentEinsumStrides(strides,
                                      prod_einsum,
                                      cons_einsum,
                                      dspace,
                                      workload);
      }
    }
  }

  return strides;
}

void GatherSubsequentEinsumStrides(
  std::map<problem::DimensionId, std::map<problem::DimensionId, int>>& strides,
  problem::EinsumId prod_einsum,
  problem::EinsumId cons_einsum,
  problem::DataSpaceId dspace,
  const problem::FusedWorkload& workload
)
{
  for (auto [prod_dim, prod_dim_i] : workload.EinsumDimToIdx(prod_einsum))
  {
    for (auto [dspace_dim, dspace_dim_i] : workload.DspaceDimToIdx(dspace))
    {
      for (auto [cons_dim, cons_dim_i] : workload.EinsumDimToIdx(cons_einsum))
      {
        strides[prod_dim][cons_dim] =
          isl_val_get_num_si(isl_aff_get_coefficient_val(
            workload.ReadAccessesAff(cons_einsum, dspace)
              .get_at(dspace_dim_i)
              .get(),
            isl_dim_in,
            cons_dim_i
          ));
      }
    }
  }
}

} // namespace analysis