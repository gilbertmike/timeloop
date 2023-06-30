#include "loop-analysis/mapping-to-isl/fused-mapping-to-isl.hpp"

namespace analysis
{

void ComputeLatency::CalculateLatency(LatencyAggregator&)
{
}

void SequentialLatency::CalculateLatency(LatencyAggregator& agg)
{
  auto n_dims = std::optional<size_t>{};
  isl_pw_qpolynomial* p_latency = nullptr;
  for (auto child_id : children)
  {
    auto& child = agg.AggregatorAt(child_id);
    auto p_child_latency = std::visit(
      [&agg](auto&& child_agg) {
        child_agg.CalculateLatency(agg);
        return isl_pw_qpolynomial_copy(child_agg.latency);
      },
      child
    );

    auto child_n_dims = isl_pw_qpolynomial_dim(p_child_latency, isl_dim_in);
    if (n_dims && *n_dims != child_n_dims)
    {
      throw std::logic_error("mismatched latency isl_dim_in");
    }
    else
    {
      n_dims = child_n_dims;
    }
    if (p_latency)
    {
      p_latency = isl_pw_qpolynomial_
    }
    else
    {
      p_latency = p_child_latency;
    }
  }
}

AggregatorTypes& LatencyAggregator::AggregatorAt(LatencyId id)
{
  return aggregators.at(id);
}

void LatencyAggregator::CalculateLatency()
{
  std::get<RootLatency>(aggregators.at(root)).CalculateLatency(*this);
}

LatencyId LatencyAggregator::GetRootId() const
{
  return root;
}

void LatencyAggregator::SetComputeLatency(mapping::NodeID compute,
                                          LatencyId latency)
{
  compute_to_aggregator.emplace(compute, latency);
}

void LatencyAggregator::SetLatency(mapping::NodeID compute,
                                   __isl_take isl_pw_qpolynomial* latency)
{
  auto compute_latency_id = compute_to_aggregator.at(compute);
  auto& compute_latency =
    std::get<ComputeLatency>(AggregatorAt(compute_latency_id));
  compute_latency.latency = latency;
}

LatencyAggregator
CreateLatencyAggregatorFromMapping(mapping::FusedMapping& mapping)
{
  using namespace problem;

  auto aggregator = LatencyAggregator();
  auto agg_root = aggregator.GetRootId();

  auto root = mapping.GetRoot().id;
  auto dfs_stack = std::vector<std::pair<mapping::NodeID, LatencyId>>();
  dfs_stack.emplace_back(root, agg_root);

  // auto start_idx = 0;

  while (dfs_stack.size() > 0)
  {
    auto node_id = dfs_stack.back().first;
    auto cur_agg_id = dfs_stack.back().second;
    dfs_stack.pop_back();
    const auto& node = mapping.NodeAt(node_id);

    std::visit(
      [&dfs_stack, &cur_agg_id, &aggregator] (auto&& node)
      {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, mapping::Pipeline>)
        {
          auto new_aggregator =
            aggregator.AddChild<PipelineLatency>(cur_agg_id);
          for (const auto& child : node.children)
          {
            dfs_stack.emplace_back(child, new_aggregator.id);
          }
        }
        else if constexpr (std::is_same_v<T, Sequential>)
        {
          auto new_aggregator =
            aggregator.AddChild<SequentialLatency>(cur_agg_id);
          for (const auto& child : node.children)
          {
            dfs_stack.emplace_back(child, new_aggregator.id);
          }
        }
        else if constexpr (std::is_same_v<T, mapping::Compute>)
        {
          auto new_agg = aggregator.AddChild<ComputeLatency>(cur_agg_id);
          aggregator.SetComputeLatency(node.id, new_agg.id);
        }
        else if constexpr (mapping::HasOneChildV<T>)
        {
          dfs_stack.emplace_back(*node.child, cur_agg_id);
        }
        else
        {
          throw std::logic_error("unknown mapping node type");
        }
      },
      node
    );
  }

  return aggregator;
}

}; // namespace analysis