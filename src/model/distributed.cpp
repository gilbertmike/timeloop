#include "isl-wrapper/ctx-manager.hpp"
#include "model/distributed.hpp"

#include <map>

namespace distributed {
NocSpec noc_from_yaml(std::string name, const YAML::Node& root) {
    // Gets the topology type.
    const std::string topo = root["topo"].as<std::string>();
    // Gets the dimensions of the topology.
    std::vector<std::string> dims = get_dims(root);
    // Gets the constraints of the topology.
    isl::set constraints = get_constraints(
        name, dims, root["constraints"].as<std::vector<std::string>>()
    );
    // Formulates the noc cost function based on the topology.
    isl::pw_multi_aff noc_cost;
    if (topo == "Mesh") {
        noc_cost = mesh_noc(dims, constraints);
    } else if (topo == "Torus") {
        noc_cost = torus_noc(dims, constraints);
    } else {
        throw std::runtime_error("Unknown topology type: " + topo);
    }

    return NocSpec{name, dims, noc_cost, constraints};
}

std::vector<std::string> get_dims(const YAML::Node& topology) {
    // Parses the topology data to get the dimensions.
    std::vector<std::string> dims;
    for (const auto& dim : topology["dims"]) {
        dims.push_back(dim.as<std::string>());
    }
    return dims;
}

// Constrains a domain based on a list of dimensions and constraints.
isl::set get_constraints(
    const std::string& domain,
    const std::vector<std::string>& dims, 
    const std::vector<std::string>& constraints
) {
    // Creates a string that models the constraints.
    std::string domain_str = "{" + domain + "[";
    for (size_t i = 0; i < dims.size(); ++i) {
        domain_str += dims[i];
        if (i < dims.size() - 1) {
            domain_str += ", ";
        } else {
            domain_str += "] : ";
        }
    }

    // Adds the constraints to the domain string.
    for (size_t i = 0; i < constraints.size(); ++i) {
        domain_str += constraints[i];
        if (i < constraints.size() - 1) {
            domain_str += " and ";
        } else {
            domain_str += "}";
        }
    }

    return isl::set(GetIslCtx(), domain_str);
}

/**
 * Defines the n-dimensional Manhattan distance function. This is done programatically
 * as ISL does not have an absolute value function.
 * 
 * @param dims      A vector of strings representing the dimensions.
 * 
 * @return          A piecewise affine function string representing the Manhattan distance.
 */
isl::pw_multi_aff mesh_noc(const std::vector<std::string> dims, const isl::set& constraints)
{
    // Allocates computer memory for the isl space where dist calculations are done.
    isl::space noc_space = constraints.get_space().map_from_set();
    isl::pw_aff noc_dist = noc_space.wrap().zero_aff_on_domain();
    // Converts domain into a local space.
    isl_local_space *p_dist_local = isl_local_space_from_space(isl_space_wrap(noc_space.copy()));
    // Constructs all the absolute value affines per dimension and adds to metric.
    for (int i = 0; i < dims.size(); ++i)
    {
        // Constructs the affine for the first point.
        isl::pw_aff p1_aff = isl::manage(isl_pw_aff_var_on_domain(
            isl_local_space_copy(p_dist_local), isl_dim_set, i
        ));
        // Constructs the affine for the second point.
        isl::pw_aff p2_aff = isl::manage(isl_pw_aff_var_on_domain(
            isl_local_space_copy(p_dist_local), isl_dim_set, dims.size() + i
        ));

        // Subtracts for the dst. 
        isl::pw_aff diff = p1_aff.sub(p2_aff);
        // Grabs the negation.
        isl::pw_aff neg_diff = diff.neg();
        // Constructs the affine for the absolute value.
        isl::pw_aff abs_diff = diff.max(neg_diff);

        // Adds the absolute value affine to the vector.
        noc_dist = noc_dist.add(abs_diff);
    }
    
    // Labels the scalar output as hops.
    isl::id hops_id = isl::id(GetIslCtx(), "hops");

    return noc_dist.set_range_tuple(hops_id);
}

isl::pw_multi_aff torus_noc (const std::vector<std::string> dims, const isl::set& constraints)
{
    throw std::runtime_error("Not implemented yet.");
}

PhysicalSpec physical_spec_from_yaml(
    std::string name, const YAML::Node& root
) {
    // Gets the physical spec name.
    std::string spec_name = root["name"].as<std::string>();
    
    // Reads in all the network specifications first.
    std::unordered_map<std::string, std::shared_ptr<const NocSpec>> noc_specs;
    for (const auto& node: root)
    {
        // Checks that it has the !Network tag.
        if (node.Tag() == "!Network")
        {
            // Gets the noc spec.
            const NocSpec noc_spec = noc_from_yaml(
                node.first.as<std::string>(), node.second
            );
            // Ensures each NocSpec has a unique name.
            if (noc_specs.find(noc_spec.name) != noc_specs.end())
            {
                throw std::runtime_error("Noc spec already exists: " + noc_spec.name);
            }
            // Adds the noc spec to the map.
            noc_specs[noc_spec.name] = std::make_shared<const NocSpec>(noc_spec);
        }
    }
    
    // Reads in the rest of the components.
    std::vector<PhysicalComponent> components;
    for (const auto& node: root)
    {
        // Checks that it has the !Component tag.
        if (node.Tag() == "!Component")
        {
            // Gets the component name.
            std::string comp_name = node["name"].as<std::string>();
            // Fetches the network information.
            auto network_node = node["network"];
            std::string noc_name = network_node["name"].as<std::string>();
            isl::map placement_map = isl::map(
                GetIslCtx(), network_node["placement"].as<std::string>()
            );

            // Gets the noc spec.
            auto noc_it = noc_specs.find(noc_name);
            if (noc_it == noc_specs.end()) {
                throw std::runtime_error("Noc spec not found: " + noc_name);
            }
            const NocSpec noc_spec = *noc_it->second;

            // Creates the placement spec.
            PlacementSpec placement_spec(noc_spec, placement_map);

            // Creates the component.
            PhysicalComponent component{
                comp_name,
                placement_spec
            };
            // Adds the component to the vector.
            components.push_back(component);
        }
    }
}
}  // namespace distributed