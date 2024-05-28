#include "binding/binding.hpp"

IslBinding BindingToIsl(const BindingSpec& spec)
{
  IslBinding binding = IslBinding{.virtual_buffer=spec.virtual_buffer,
                                  .physical_buffer=spec.physical_buffer};
  
  for (const auto& desc : spec.descriptors)
  {
    // add physical spatial dimensions to mapping
    if (std::holds_alternative<SpatialPartition>(desc))
    {
      // add spatial dimensions to mapping
    }
  }

  return binding;
}

std::ostream& operator<<(std::ostream& os, const SpatialPartition& sp)
{
  os << "Partition(" << ".rank=" << sp.rank << ", ";
  os << ".stride=" << sp.stride << ", ";
  os << ".tile_size=" << sp.tile_size << ", ";
  os << ".num_tiles=" << sp.num_tiles << ", ";
  os << ".spatial_dimension=" << sp.spatial_dimension << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const SpatialDuplicate& sd)
{
  os << "Duplicate(.num_duplication=" << sd.num_duplication << ", ";
  os << ".spatial_dimension=" << sd.spatial_dimension << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const BindingDescriptor& desc)
{
  std::visit([&os](auto&& arg) { os << arg; }, desc);
  return os;
}

std::ostream& operator<<(std::ostream& os, const BindingSpec& spec)
{
  os << "Binding(.virtual_buffer=" << spec.virtual_buffer << ", ";
  os << ".physical_buffer=" << spec.physical_buffer << ", ";
  os << ".descriptors={";

  bool first = true;
  for (const auto& desc : spec.descriptors)
  {
    if (!first)
    {
      os << ", ";
    }
    os << desc;
  }
  os << "})";

  return os;
}