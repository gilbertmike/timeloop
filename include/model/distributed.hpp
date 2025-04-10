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
    const isl::pw_multi_aff noc_cost;
    const isl::set domain;
};

struct PlacementSpec { // Placement specification of a component on a NoC
    const NocSpec &noc_spec;
    const isl::map placement_map;

    /// @brief Validates the placement map upon construction.
    PlacementSpec(
        const NocSpec &noc_spec, const isl::map &placement_map
    ) : noc_spec(noc_spec), placement_map(placement_map) {
        // Extracts the domain from the placement map.
        isl::set domain = noc_spec.domain;
        // Checks if the range of the placement map is a subset of the domain.
        if (!placement_map.range().is_subset(domain)) {
            throw std::runtime_error("Placement map range is not a subset of the domain.");
        }
    }
};

struct PhysicalComponent { // Component of physical spec attached to NOC
    const std::string name;
    const PlacementSpec placement;
};
// Physical component aliases for static analysis.
typedef PhysicalComponent PhysicalStorage;
typedef PhysicalComponent PhysicalProcessingElement;
// Ports connecting different NoC components.
struct Port { // Port of the physical component
    const std::shared_ptr<NocSpec> parent;
    const std::shared_ptr<NocSpec> child;
    isl::map port_map;
};

struct PhysicalSpec { // Physical specification
    const std::string name;
    const std::vector<PhysicalComponent> components;
    const std::vector<Port> ports;
};

struct LogicalComponent { // Logical component of the accelerator.
    const std::string name;
    isl::set domain;
};
// Logical component aliases for static analysis.
typedef LogicalComponent LogicalStorage;
typedef LogicalComponent LogicalProcessingElement;

struct LogicalSpec {
    const std::string name;
    const std::vector<LogicalComponent> components;
};

struct BindingSpec {
    const std::string name;
    const std::shared_ptr<PhysicalSpec> physical_spec;
    const std::shared_ptr<LogicalSpec> logical_spec;
    const std::vector<isl::map> bindings;
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