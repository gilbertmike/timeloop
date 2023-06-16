#include "loop-analysis/isl-ir.hpp"

namespace analysis {

Skew::Skew(const std::vector<spacetime::Dimension>& dim_in_tags,
                     isl::map map) :
  dim_in_tags(dim_in_tags), map(std::move(map))
{
}

Occupancy::Occupancy(const std::vector<spacetime::Dimension>& dim_in_tags,
                     isl::map map) :
  dim_in_tags(dim_in_tags), map(std::move(map))
{
}

Fill::Fill(const std::vector<spacetime::Dimension>& dim_in_tags,
           isl::map map) :
  dim_in_tags(dim_in_tags), map(std::move(map))
{
}

Transfers::Transfers(const std::vector<spacetime::Dimension>& dim_in_tags,
                     isl::map map) :
  dim_in_tags(dim_in_tags), map(std::move(map))
{
}

}; // namespace analysis