#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <isl/cpp.h>

#include "compound-config/compound-config.hpp"
#include "yaml.h"

namespace distributed {
struct Noc;

struct NocSpec {
    const std::string name;
    const std::vector<std::string> dims;
    isl::pw_multi_aff noc_cost;
    isl::set domain;
};

struct PlacementSpec {
    const NocSpec &parent;
    const isl::map placement;   // ArrayID -> NoC
};

NocSpec noc_from_yaml(std::string name, const YAML::Node& root);
std::vector<std::string> get_dims(const YAML::Node& topology);
isl::set get_constraints(
    const std::string& name, const std::vector<std::string>& dims, 
    const std::vector<std::string>& constraints
);
isl::pw_multi_aff mesh_noc(const std::vector<std::string> dims, const isl::set& constraints);
isl::pw_multi_aff torus_noc(const std::vector<std::string> dims, const isl::set& constraints);
};  // namespace distributed