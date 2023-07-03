#include "loop-analysis/isl-analysis/capacity-analysis.hpp"

#include <barvinok/isl.h>
#include <boost/range/adaptor/reversed.hpp>

#include "isl-wrapper/isl-functions.hpp"

namespace analysis
{

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

}; // namespace analysis