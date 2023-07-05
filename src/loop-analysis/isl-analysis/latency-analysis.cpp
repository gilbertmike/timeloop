#include "loop-analysis/isl-analysis/latency-analysis.hpp"

#include "isl-wrapper/isl-functions.hpp"

namespace analysis
{

void ComputeLatency::CalculateLatency(LatencyAggregator&)
{
}

void RootLatency::CalculateLatency(LatencyAggregator& agg)
{
  latency = nullptr;
  auto& child = agg.AggregatorAt(child_id);
  auto p_child_latency = std::visit(
    [&agg](auto&& child_agg) {
      child_agg.CalculateLatency(agg);
      return isl_pw_qpolynomial_copy(child_agg.latency);
    },
    child
  );

  auto n_dims = isl_pw_qpolynomial_dim(p_child_latency, isl_dim_in);
  latency = p_child_latency;

  auto p_projector = isl::dim_projector(
    isl_pw_qpolynomial_get_domain_space(latency),
    0,
    n_dims
  );

  latency = isl_map_apply_pw_qpolynomial(p_projector, latency);
  std::cout << "latency: " << isl_pw_qpolynomial_to_str(latency) << std::endl;
}

void PipelineLatency::CalculateLatency(LatencyAggregator& agg)
{
  auto n_dims = std::optional<int>{};
  latency = nullptr;
  for (auto child_id : children_id)
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
    if (latency)
    {
      latency = isl_pw_qpolynomial_add(latency, p_child_latency);
    }
    else
    {
      latency = p_child_latency;
    }
  }

  auto p_domain = isl_pw_qpolynomial_domain(isl_pw_qpolynomial_copy(latency));

  auto p_identity = isl_map_identity(
    isl_space_map_from_set(isl_set_get_space(p_domain))
  );
  p_identity = isl_map_intersect_domain(p_identity, isl_set_copy(p_domain));

  // Latency should be in the form [..., (PipelineSequential)*, PipelineSpatial]
  auto p_next_sequential = isl::map_to_next(
    isl_set_copy(p_domain),
    start_idx,
    *n_dims-start_idx
  );

  std::cout << isl_map_to_str(p_next_sequential) << std::endl;

  auto p_projector = isl::dim_projector(
    isl_pw_qpolynomial_get_domain_space(latency),
    start_idx,
    *n_dims - start_idx
  );

  latency = isl_map_apply_pw_qpolynomial(p_projector, latency);
  std::cout << "latency: " << isl_pw_qpolynomial_to_str(latency) << std::endl;
}

void SequentialLatency::CalculateLatency(LatencyAggregator& agg)
{
  auto n_dims = std::optional<int>{};
  latency = nullptr;
  for (auto child_id : children_id)
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
    if (latency)
    {
      latency = isl_pw_qpolynomial_add(latency, p_child_latency);
    }
    else
    {
      latency = p_child_latency;
    }
  }

  auto p_projector = isl::dim_projector(
    isl_pw_qpolynomial_get_domain_space(latency),
    start_idx,
    *n_dims - start_idx
  );

  latency = isl_map_apply_pw_qpolynomial(p_projector, latency);
  std::cout << "latency: " << isl_pw_qpolynomial_to_str(latency) << std::endl;
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

LatencyAggregator::LatencyAggregator() :
  root(0)
{
  aggregators.emplace_back(RootLatency());
}

LatencyAggregator
CreateLatencyAggregatorFromMapping(mapping::FusedMapping& mapping)
{
  using namespace problem;

  auto aggregator = LatencyAggregator();
  auto agg_root = aggregator.GetRootId();

  auto root = mapping.GetRoot().id;
  auto dfs_stack =
    std::vector<std::tuple<mapping::NodeID, LatencyId, size_t>>();
  dfs_stack.emplace_back(root, agg_root, 0);

  // auto start_idx = 0;

  while (dfs_stack.size() > 0)
  {
    auto node_id = std::get<0>(dfs_stack.back());
    auto cur_agg_id = std::get<1>(dfs_stack.back());
    auto cur_start_idx = std::get<2>(dfs_stack.back());
    dfs_stack.pop_back();
    const auto& node = mapping.NodeAt(node_id);

    std::visit(
      [&dfs_stack, &cur_agg_id, &aggregator, &cur_start_idx] (auto&& node)
      {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, mapping::Pipeline>)
        {
          auto new_aggregator =
            aggregator.AddChild<PipelineLatency>(cur_agg_id);
          new_aggregator.start_idx = cur_start_idx;
          for (const auto& child : node.children)
          {
            dfs_stack.emplace_back(child, new_aggregator.id, cur_start_idx+1);
          }
        }
        else if constexpr (std::is_same_v<T, mapping::Sequential>)
        {
          auto new_aggregator =
            aggregator.AddChild<SequentialLatency>(cur_agg_id);
          new_aggregator.start_idx = cur_start_idx;
          for (const auto& child : node.children)
          {
            dfs_stack.emplace_back(child, new_aggregator.id, cur_start_idx+1);
          }
        }
        else if constexpr (std::is_same_v<T, mapping::Compute>)
        {
          auto new_agg = aggregator.AddChild<ComputeLatency>(cur_agg_id);
          aggregator.SetComputeLatency(node.id, new_agg.id);
        }
        else if constexpr (mapping::IsLoopV<T>)
        {
          dfs_stack.emplace_back(*node.child, cur_agg_id, cur_start_idx+1);
        }
        else if constexpr (mapping::HasOneChildV<T>)
        {
          dfs_stack.emplace_back(*node.child, cur_agg_id, cur_start_idx);
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