#include <vector>

#include <barvinok/isl.h>

#include "loop-analysis/isl-ir.hpp"
#include "mapping/fused-mapping.hpp"
#include "workload/fused-workload.hpp"

namespace analysis 
{

using CapacityId = size_t;
struct CapacityAggregator;

struct CapacityBase
{
  virtual void CalculateCapacity(CapacityAggregator& agg) = 0;

  CapacityId id;
  std::map<std::pair<mapping::BufferId, problem::DataSpaceId>,
           isl_map*> occupancy;
};

struct ComputeCapacity : public CapacityBase
{
  size_t start_idx;

  void CalculateCapacity(CapacityAggregator&) override;
};

struct RootCapacity : public CapacityBase
{
  virtual void CalculateCapacity(CapacityAggregator& agg) override;
  CapacityId child_id;
  std::map<mapping::BufferId, isl_pw_qpolynomial*> final_capacity;
};

struct BranchCapacityBase : public CapacityBase
{
  size_t start_idx;
  std::vector<CapacityId> children_id;
};

struct PipelineCapacity : public BranchCapacityBase
{
  void CalculateCapacity(CapacityAggregator&) override;
};

struct SequentialCapacity : public BranchCapacityBase
{
  void CalculateCapacity(CapacityAggregator&) override;
};

template<typename T>
inline constexpr bool HasOneChildV = std::is_same_v<T, RootCapacity>;

template<typename T>
inline constexpr bool HasManyChildrenV = std::is_same_v<T, PipelineCapacity>
                                       | std::is_same_v<T, SequentialCapacity>;

using AggregatorTypes = std::variant<ComputeCapacity,
                                     RootCapacity,
                                     PipelineCapacity,
                                     SequentialCapacity>;

struct CapacityAggregator
{
  CapacityAggregator();

  AggregatorTypes& AggregatorAt(CapacityId id);

  void CalculateCapacity();

  CapacityId GetRootId() const;

  template<typename T>
  T& AddChild(CapacityId parent)
  {
    CapacityId child_id = aggregators.size();
    auto& child = std::get<T>(aggregators.emplace_back(T()));
    child.id = child_id;

    std::visit(
      [&child_id](auto&& agg)
      {
        using ParentT = std::decay_t<decltype(agg)>;
        if constexpr (HasOneChildV<ParentT>)
        {
          agg.child_id = child_id;
        }
        else if constexpr (HasManyChildrenV<ParentT>)
        {
          agg.children_id.push_back(child_id);
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

 private:
  CapacityId root;
  std::vector<AggregatorTypes> aggregators;
};

CapacityAggregator
CreateCapacityAggregatorFromMapping(mapping::FusedMapping& mapping);

std::map<mapping::BufferId, isl_pw_qpolynomial*>
ComputeCapacityFromMapping(
  mapping::FusedMapping& mapping,
  const std::map<LogicalBuffer, Occupancy>& occupancies
);

};