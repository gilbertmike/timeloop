#pragma once

#include <isl/polynomial.h>

#include "loop-analysis/isl-ir.hpp"

namespace analysis
{

using LatencyId = size_t;
struct LatencyAggregator;

struct LatencyBase
{
  virtual void CalculateLatency(LatencyAggregator& agg) = 0;

  LatencyId id;
  isl_pw_qpolynomial* latency;
};

struct ComputeLatency : public LatencyBase
{
  void CalculateLatency(LatencyAggregator&) override;
};

struct RootLatency : public LatencyBase
{
  LatencyId child;
};

struct BranchLatencyBase : public LatencyBase
{
  size_t start_idx;
  std::vector<LatencyId> children;
};

struct PipelineLatency : public BranchLatencyBase
{
  void CalculateLatency(LatencyAggregator&) override;
};

struct SequentialLatency : public BranchLatencyBase
{
  void CalculateLatency(LatencyAggregator&) override;
};

template<typename T>
inline constexpr bool HasOneChildV = std::is_same_v<T, RootLatency>;

template<typename T>
inline constexpr bool HasManyChildrenV = std::is_same_v<T, PipelineLatency>
                                       | std::is_same_v<T, SequentialLatency>;

using AggregatorTypes = std::variant<ComputeLatency,
                                     RootLatency,
                                     PipelineLatency,
                                     SequentialLatency>;

struct LatencyAggregator
{
  AggregatorTypes& AggregatorAt(LatencyId id);

  void CalculateLatency();

  LatencyId GetRootId() const;

  template<typename T>
  T& AddChild(LatencyId parent)
  {
    LatencyId child_id = aggregators.size();
    auto& child = std::get<T>(aggregators.emplace_back(T()));
    child.id = child_id;

    std::visit(
      [&child_id](auto&& agg)
      {
        using ParentT = std::decay_t<decltype(agg)>;
        if constexpr (HasOneChildV<ParentT>)
        {
          agg.child = child_id;
        }
        else if constexpr (HasManyChildrenV<ParentT>)
        {
          agg.children.push_back(child_id);
        }
        else
        {
          throw std::logic_error("cannot add child to this node");
        }
      },
      AggregatorAt(parent)
    );
    return child;
  }

  void SetComputeLatency(mapping::NodeID compute, LatencyId latency);

  void SetLatency(mapping::NodeID compute,
                  __isl_take isl_pw_qpolynomial* latency);

 private:
  LatencyId root;
  std::vector<AggregatorTypes> aggregators;
  std::map<mapping::NodeID, LatencyId> compute_to_aggregator;
};

LatencyAggregator
CreateLatencyAggregatorFromMapping(mapping::FusedMapping& mapping);

/**
 * @brief Results of mapping analysis that will become input into reuse
 *   analysis.
 */
struct MappingAnalysisResult
{
  /**
   * @brief The occupancy of every logical buffer as defined in the mapping.
   */
  std::map<LogicalBuffer, Occupancy> lbuf_to_occupancy;
  /**
   * @brief Tiling of each branch. The tiling is a relation between tiling
   *   variables and operations.
   * 
   * An uncompletely tiled branch will have multiple-valued isl::map.
   */
  std::map<mapping::NodeID, isl::map> branch_tiling;
  /**
   * @brief We can assume an amount of parallelism to quickly calculate approx.
   *   compute latency by simply dividing number of operations with assumed
   *   parallelism.
   */
  std::map<mapping::NodeID, double> compute_to_assumed_parallelism;
  /**
   * @brief An aggregator to calculate compute latency given branch latencies.
   */
  LatencyAggregator compute_latency_aggregator;
};

/**
 * @brief Compute occupancy, tiling, and other intermediate representations
 *   from the mapping.
 * 
 * @see analysis::MappingAnalysisResult
 */
MappingAnalysisResult
OccupanciesFromMapping(mapping::FusedMapping& mapping,
                       const problem::FusedWorkload& workload);

}; // namespace analysis
