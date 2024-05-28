#pragma once

#include <variant>
#include <isl/cpp.h>
#include <iostream>
#include "loop-analysis/isl-ir.hpp"

struct IslBinding
{
  analysis::BufferId virtual_buffer;
  analysis::BufferId physical_buffer;
  isl::map map;
};

struct SpatialPartition
{
  int rank;

  int stride;
  int tile_size;
  int num_tiles;

  int spatial_dimension;
};

struct SpatialDuplicate
{
  int num_duplication;
  int spatial_dimension;
};

using BindingDescriptor = std::variant<SpatialPartition, SpatialDuplicate>;

struct BindingSpec
{
  analysis::BufferId virtual_buffer;
  analysis::BufferId physical_buffer;

  std::vector<BindingDescriptor> descriptors;
};

IslBinding BindingToIsl(const BindingSpec& spec);

std::ostream& operator<<(std::ostream&, const SpatialPartition&);
std::ostream& operator<<(std::ostream&, const SpatialDuplicate&);
std::ostream& operator<<(std::ostream&, const BindingDescriptor&);
std::ostream& operator<<(std::ostream&, const BindingSpec&);