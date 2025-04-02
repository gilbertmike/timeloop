#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <isl/cpp.h>

#include "compound-config/compound-config.hpp"
#include "yaml.h"

namespace distributed {

struct NocSpec { // NoC specification
    const std::string name;
    const std::vector<std::string> dims;
    isl::pw_multi_aff noc_cost;
    isl::set domain;
};

// Placement specification of a component on a NoC.
typedef std::pair<const std::shared_ptr<NocSpec>, isl::map> PlacementSpec;
struct PhysicalComponent { // Component of physical spec attached to NOC
    const std::string name;
    const std::string type;
    PlacementSpec placement;
};

struct PhysicalSpec { // Physical specification
    const std::string name;
    const std::vector<PhysicalComponent> components;
};

struct LogicalComponent { // Logical component of the accelerator.
    const std::string name;
    const std::string type;
    isl::set domain;
};

struct LogicalSpec {
    const std::string name;
    const std::vector<LogicalComponent> components;
};

struct BindingSpec {
    const std::string name;
    const PhysicalSpec physical_spec;
    const LogicalSpec logical_spec;
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