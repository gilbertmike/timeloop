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
}