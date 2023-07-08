#include "loop-analysis/isl-analysis/capacity-analysis.hpp"

#include <barvinok/isl.h>
#include <boost/range/adaptor/reversed.hpp>

#include "isl-wrapper/isl-functions.hpp"

namespace analysis
{

std::map<mapping::BufferId, isl_pw_qpolynomial*>
ComputeCapacityFromMapping(
  mapping::NodeID cur_node_id,
  mapping::FusedMapping& mapping,
  const std::map<LogicalBuffer, Occupancy>& occupancies
);

std::map<mapping::BufferId, isl_pw_qpolynomial*>
ComputeCapacityFromMapping(
  mapping::FusedMapping& mapping,
  const std::map<LogicalBuffer, Occupancy>& occupancies,
  const std::map<mapping::NodeID, std::vector<LogicalBuffer>>& node_to_lbuf
)
{
  using namespace boost::adaptors;
  using namespace mapping;

  (void) occupancies;

  auto result = std::map<mapping::BufferId, isl_pw_qpolynomial*>();

  std::map<mapping::NodeID, std::vector<isl_pw_qpolynomial*>> branch_to_caps;

  auto dfs_order = IterateInDfsOrder(mapping);
  auto dfs_node_id_list = std::vector<std::tuple<NodeID, NodeID, size_t>>();
  for (auto node_parent : IterateInDfsOrder(mapping))
  {
    dfs_node_id_list.emplace_back(node_parent);
  }

  for (auto [node_id, parent_id, n_loops] : dfs_node_id_list | reversed )
  {
    auto it = node_to_lbuf.find(node_id);
    if (it == node_to_lbuf.end())
    {
      continue;
    }
    const auto& lbufs = it->second;
    for (const auto& lbuf : lbufs)
    {
      (void) lbuf;
      // const auto& occ = occupancies.at(lbuf).map;
    }

    // Aggregate children capacities

    // Aggregate with current capacities
  }

  return result;
}

std::map<mapping::BufferId, isl_pw_qpolynomial*>
ComputeCapacityFromMapping(
  mapping::FusedMapping& mapping,
  const std::map<LogicalBuffer, Occupancy>& occupancies
)
{
  return
    ComputeCapacityFromMapping(mapping.GetRoot().id, mapping, occupancies);
}

struct BufferInfo
{
  BufferId buffer_id;
  DataSpaceID dataspace_id;
  bool exploits_reuse;
};

std::map<mapping::BufferId, isl_pw_qpolynomial*>
ComputeCapacityFromMapping(
  mapping::NodeID cur_node_id,
  mapping::FusedMapping& mapping,
  const std::map<LogicalBuffer, Occupancy>& occupancies
)
{
  (void) occupancies;
  std::vector<BufferInfo> buffers;
  auto keep_going = true;
  while (keep_going)
  {
    const auto& node = mapping.NodeAt(cur_node_id);
    std::visit(
      [&cur_node_id, &buffers, &keep_going, &mapping, &occupancies](auto&& node)
      {
        using T = std::decay_t<decltype(node)>;

        if constexpr (std::is_same_v<T, mapping::Storage>)
        {
          buffers.emplace_back(
            BufferInfo{node.buffer, node.dspace, node.exploits_reuse}
          );
        }
        else if constexpr (std::is_same_v<T, mapping::Compute>)
        {
          // Do something here
          keep_going = false;
        }
        else if constexpr (mapping::IsBranchV<T>)
        {
          for (auto child_id : node.children)
          {
            ComputeCapacityFromMapping(child_id, mapping, occupancies);
          }

          if constexpr (std::is_same_v<T, mapping::Sequential>)
          {
            // Do something here
          }
          else if constexpr (std::is_same_v<T, mapping::Pipeline>)
          {
            // Do something here
          }
          else
          {
            throw std::logic_error("unknown mapping branch node");
          }
          keep_going = false;
        }
        else if constexpr (mapping::HasOneChildV<T>)
        {
          cur_node_id = *node.child;
        }
        else 
        {
          throw std::logic_error("unknown mapping node");
        }
      },
      node
    );
  }
}
}; // namespace analysis