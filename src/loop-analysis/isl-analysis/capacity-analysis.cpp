#include "loop-analysis/isl-analysis/capacity-analysis.hpp"

#include <barvinok/isl.h>
#include <boost/range/adaptor/reversed.hpp>

#include "isl-wrapper/isl-functions.hpp"

namespace analysis
{

struct BufferInfo
{
  BufferId buffer_id;
  DataSpaceID dataspace_id;
  bool exploits_reuse;
  bool right_above_branch;
};

struct IntermediateResult
{
  std::set<EinsumID> einsums;
};

IntermediateResult ComputeCapacityFromMapping(
  mapping::NodeID cur_node_id,
  mapping::FusedMapping& mapping,
  const std::map<LogicalBuffer, Occupancy>& occupancies
);

std::map<mapping::BufferId, isl_pw_qpolynomial*>
ComputeCapacityFromMapping(
  mapping::FusedMapping& mapping,
  const std::map<LogicalBuffer, Occupancy>& occupancies
)
{
  ComputeCapacityFromMapping(mapping.GetRoot().id, mapping, occupancies);
  return std::map<mapping::BufferId, isl_pw_qpolynomial*>();
}

IntermediateResult ComputeCapacityFromMapping(
  mapping::NodeID cur_node_id,
  mapping::FusedMapping& mapping,
  const std::map<LogicalBuffer, Occupancy>& occupancies
)
{
  (void) occupancies;
  IntermediateResult result;
  std::vector<BufferInfo> buffers;
  auto keep_going = true;
  while (keep_going)
  {
    const auto& node = mapping.NodeAt(cur_node_id);
    std::visit(
      [&cur_node_id, &buffers, &keep_going, &result, &mapping, &occupancies]
      (auto&& node)
      {
        using T = std::decay_t<decltype(node)>;

        if constexpr (std::is_same_v<T, mapping::Storage>)
        {
          buffers.emplace_back(
            BufferInfo{node.buffer, node.dspace, node.exploits_reuse, true}
          );
        }
        else if constexpr (std::is_same_v<T, mapping::Compute>)
        {
          result.einsums.insert(node.kernel);

          for (const auto& b : buffers)
          {
            // Do something here
          }
        }
        else if constexpr (mapping::IsBranchV<T>)
        {
          for (auto child_id : node.children)
          {
            auto child_res =
              ComputeCapacityFromMapping(child_id, mapping, occupancies);
            result.einsums.merge(std::move(child_res.einsums));
          }

          if constexpr (std::is_same_v<T, mapping::Sequential>)
          {
            for (const auto& b : buffers)
            {
              if (!b.exploits_reuse && b.right_above_branch)
              {
                // Do something here
              }
              else
              {
                // Do something here
              }
            }
          }
          else if constexpr (std::is_same_v<T, mapping::Pipeline>)
          {
            for (const auto& b : buffers)
            {
              // Do something here
            }
          }
          else
          {
            throw std::logic_error("unknown mapping branch node");
          }
        }
        else if constexpr (mapping::IsLoopV<T>)
        {
          for (auto& b : buffers)
          {
            b.right_above_branch = false;
          }
        }

        if constexpr(HasOneChildV<T>)
        {
          cur_node_id = *node.child;
        }
        else
        {
          keep_going = false;
        }
      },
      node
    );
  }

  return result;
}
}; // namespace analysis