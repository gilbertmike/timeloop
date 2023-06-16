#pragma once

#include <utility>

#include "loop-analysis/nest-analysis-tile-info.hpp"
#include "loop-analysis/isl-nest-analysis.hpp"

namespace analysis
{

std::pair<CompoundComputeNest, CompoundDataMovementNest>
GenerateLegacyNestAnalysisOutput(IslAnalysisOutput);

}; // namespace analysis