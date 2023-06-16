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
};

struct IslAnalysisOutput
{
  std::map<LogicalBuffer, LogicalBufferStats> buf_to_stats;
};

}; // namespace analysis