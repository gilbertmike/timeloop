#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <isl/cpp.h>

#include "compound-config/compound-config.hpp"
#include "yaml.h"

namespace distributed {
isl::pw_aff noc_from_yaml(std::string name, const YAML::Node& root);
std::vector<std::string> get_dims(const YAML::Node& topology);
isl::set get_constraints(
    const std::string& name, const std::vector<std::string>& dims, 
    const std::vector<std::string>& constraints
);
isl::pw_aff mesh_noc(const std::vector<std::string> dims, const isl::set& constraints);
isl::pw_aff torus_noc(const std::vector<std::string> dims, const isl::set& constraints);

struct Noc;

struct TopologySpec {
    const std::vector<std::string> dims;
    isl::pw_aff noc_cost;
    isl::set domain;
};

struct PlacementSpec {
    const Noc &parent;
    const isl::map placement;   // ArrayID -> NoC
};

struct Noc {
    const std::string name;
    TopologySpec topology;
    PlacementSpec placement;
};

struct Buffer {
    PlacementSpec placement;
    isl::map capacity;
};
};  // namespace distributed