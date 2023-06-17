#include "loop-analysis/isl-nest-analysis.hpp"

#include "loop-analysis/mapping-to-isl.hpp"
#include "loop-analysis/spatial-analysis.hpp"
#include "loop-analysis/temporal-analysis.hpp"

namespace analysis
{

LogicalBufferStats::LogicalBufferStats(const LogicalBuffer& buf) : buf(buf)
{
}

ReuseAnalysisOutput ReuseAnalysis(
  const loop::Nest& nest,
  const problem::Workload& workload
)
{
  auto occupancies = OccupanciesFromMapping(nest, workload);

  auto spatial_reuse_models = SpatialReuseModels();
  spatial_reuse_models.EmplaceLinkTransferModel<SimpleLinkTransferModel>();
  spatial_reuse_models.EmplaceMulticastModel<SimpleMulticastModel>();

  auto output = ReuseAnalysisOutput();

  for (const auto& [buf, occ] : occupancies)
  {
    auto stats = LogicalBufferStats(buf);

    auto temp_reuse_out = TemporalReuseAnalysis(
      TemporalReuseAnalysisInput(
        buf,
        occ,
        BufTemporalReuseOpts{.exploit_temporal_reuse=true}
      )
    );

    std::cout << buf << std::endl;
    std::cout << temp_reuse_out.fill.map << std::endl;
    std::cout << temp_reuse_out.effective_occupancy.map << std::endl;

    auto spatial_reuse_out = SpatialReuseAnalysis(
      SpatialReuseAnalysisInput(buf,
                                temp_reuse_out.fill,
                                temp_reuse_out.effective_occupancy),
      spatial_reuse_models
    );

    stats.occupancy = occ;
    stats.effective_occupancy = temp_reuse_out.effective_occupancy;
    stats.fill = temp_reuse_out.fill;
    stats.link_transfer = spatial_reuse_out.link_transfer_info.link_transfer;
    stats.parent_reads = spatial_reuse_out.multicast_info.reads;
    stats.compat_access_stats =
      spatial_reuse_out.multicast_info.compat_access_stats;

    output.buf_to_stats.emplace(std::make_pair(buf, std::move(stats)));
  }

  return output;
}

}; // namespace analysis