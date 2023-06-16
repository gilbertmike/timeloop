#include "loop-analysis/isl-to-legacy-adaptor.hpp"

#include "isl-wrapper/isl-functions.hpp"

namespace analysis
{

std::pair<CompoundComputeNest, CompoundDataMovementNest>
GenerateLegacyNestAnalysisOutput(IslAnalysisOutput)
{
  CompoundDataMovementNest working_sets;

  size_t num_compute_units = 1;
  for (const auto& state : nest_state_)
  {
    if (loop::IsSpatial(state.descriptor.spacetime_dimension))
    {
      num_compute_units *= state.descriptor.end;
    }
  }
  // Insert innermost level with number of iterations divided by spatial elements
  BufferID innermost_buf_id = storage_tiling_boundaries_.size();

  uint64_t max_temporal_iterations = 1;
  for (auto& state : nest_state_)
  {
    if (!loop::IsSpatial(state.descriptor.spacetime_dimension))
      max_temporal_iterations *= state.descriptor.end;
  }

  for (auto& [buf, occupancy] : occupancies)
  {
    if (buf.buffer_id == innermost_buf_id)
    {
      auto compute_info = ComputeInfo();
      compute_info.replication_factor = num_compute_units;
      std::cout << occupancy.map.domain() << std::endl;
      compute_info.accesses = isl::val_to_double(
        isl::get_val_from_singular_qpolynomial(
          isl::set_card(occupancy.map.domain())
        )
      ) / num_compute_units;
      compute_info.max_temporal_iterations = max_temporal_iterations;
      compute_info_sets_.push_back(compute_info);
      break;
    }
  }
  for (decltype(nest_state_)::size_type i = 0; i < nest_state_.size() - 1; ++i)
  {
    auto compute_info = ComputeInfo();
    compute_info_sets_.push_back(compute_info);
  }

  BufferID cur_buffer_id = storage_tiling_boundaries_.size();
  bool first_loop = true;
  bool last_boundary_found = false;
  bool should_dump = false;
  for (const auto& cur : nest_state_)
  {
    bool valid_level = !loop::IsSpatial(cur.descriptor.spacetime_dimension) || master_spatial_level_[cur.level];
    if (!valid_level)
    {
      continue;
    }

    auto is_master_spatial = master_spatial_level_[cur.level];
    auto is_boundary = storage_boundary_level_[cur.level];

    if (first_loop)
    {
      first_loop = false;
      last_boundary_found = false;
      should_dump = true;
    }
    else if (is_boundary && !last_boundary_found)
    {
      last_boundary_found = true;
      should_dump = false;
    }
    else if (is_boundary && last_boundary_found)
    {
      last_boundary_found = true;
      should_dump = true;
    }
    else if (is_master_spatial && last_boundary_found)
    {
      last_boundary_found = false;
      should_dump = true;
    }

    for (unsigned dspace_id = 0;
        dspace_id < workload_->GetShape()->NumDataSpaces;
        ++dspace_id)
    {
      DataMovementInfo tile;
      tile.link_transfers = 0;
      tile.replication_factor = num_spatial_elems_[cur.level];
      tile.fanout = logical_fanouts_[cur.level];
      tile.is_on_storage_boundary = storage_boundary_level_[cur.level];
      tile.is_master_spatial = master_spatial_level_[cur.level];

      if (should_dump)
      {
        auto& occ = eff_occupancies.at(LogicalBuffer(cur_buffer_id,
                                                      dspace_id,
                                                      0));
        const auto& key_to_access_stats =
          result.multicast_info.compat_access_stats.at(
            LogicalBuffer(cur_buffer_id, dspace_id, 0)
          );

        for (const auto& [key, access_stats] : key_to_access_stats)
        {
          tile.access_stats.stats[key] = AccessStats{
            .accesses = access_stats.accesses / tile.replication_factor,
            .hops = access_stats.hops
          };
        }

        for (const auto& [buf_ab, transfers] :
            result.link_transfer_info.link_transfers)
        {
          const auto& buf = buf_ab.first;
          if (buf.buffer_id == cur_buffer_id && buf.dspace_id == dspace_id) 
          {
            auto p_val = isl::get_val_from_singular_qpolynomial(
              isl::sum_map_range_card(transfers.map)
            );
            p_val = isl_val_div(
              p_val,
              isl_val_int_from_si(GetIslCtx().get(),
                                  num_spatial_elems_[cur.level])
            );
            tile.link_transfers = isl::val_to_double(p_val);
          }
        }

        auto p_occ_map = occ.map.copy();
        std::cout << "occ: " << isl_map_to_str(p_occ_map) << std::endl;
        auto p_occ_count = isl::get_val_from_singular_qpolynomial_fold(
          isl_pw_qpolynomial_bound(
            isl_map_card(
              isl_map_project_out(
                p_occ_map,
                isl_dim_in,
                isl_map_dim(p_occ_map, isl_dim_in)-3,
                2
              )
            ),
            isl_fold_max,
            nullptr
          )
        );
        tile.size = isl::val_to_double(p_occ_count);
      }

      working_sets[dspace_id].push_back(tile);
    }
    
    if (should_dump)
    {
      should_dump = false;
      cur_buffer_id--;
    }
  }
}

}; // namespace analysis