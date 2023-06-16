#pragma once

#include <loop-analysis/isl-ir.hpp>

namespace analysis
{

struct LogicalBufferStats
{
  const LogicalBuffer& buf;
  const LogicalBuffer& parent_buf;
  Occupancy occupancy;
  Occupancy effective_occupancy;
  Fill fill;
  Transfers link_transfer;
  Reads parent_reads;

  /***************** Compatibility with Timeloop v2.0 ************************/
  struct AccessStats
  {
    double accesses;
    double hops;
  };
  std::map<std::pair<uint64_t, uint64_t>, AccessStats> compat_access_stats;
  /***************************************************************************/
};

struct ReuseAnalysisOutput
{
  std::map<LogicalBuffer, LogicalBufferStats> buf_to_stats;
};

ReuseAnalysisOutput ReuseAnalysis(const loop::Nest&, const problem::Workload&);

}; // namespace analysis