#include "loop-analysis/isl-nest-analysis.hpp"

#include "loop-analysis/mapping-to-isl.hpp"
#include "loop-analysis/spatial-analysis.hpp"
#include "loop-analysis/temporal-analysis.hpp"

namespace analysis
{

ReuseAnalysisOutput ReuseAnalysis(
  const loop::Nest& nest,
  const problem::Workload& workload
)
{
  auto occupancies = OccupanciesFromMapping(nest, workload);

  auto spatial_reuse_models = SpatialReuseModels(SimpleLinkTransferModel(),
                                                  SimpleMulticastModel());

  for (const auto& [buf, occ] : occupancies)
  {
    auto temp_reuse_out = TemporalReuseAnalysis(
      TemporalReuseAnalysisInput(
        buf,
        occ,
        BufTemporalReuseOpts{.exploit_temporal_reuse=true}
      )
    );

    auto spatial_reuse_out = SpatialReuseAnalysis(
      SpatialReuseAnalysisInput(buf,
                                temp_reuse_out.fill,
                                static_cast<const Occupancy&>(nullptr)),
      spatial_reuse_models
    );
  }

}

}; // namespace analysis