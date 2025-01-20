#include <yaml-cpp/yaml.h>
#include <iostream>
#include <format>
#include <string>

#include <boost/test/unit_test.hpp>

#include "isl-wrapper/ctx-manager.hpp"
#include "model/distributed.hpp"

BOOST_AUTO_TEST_CASE(TestNocFromYaml)
{
  using namespace distributed;
  YAML::Node file = YAML::LoadFile("./src/unit-test/distributed-spec/distributed.yaml");
  for (const auto& root: file)
  {
    // Gets this topology spec.
    NocSpec spec = topology_from_yaml(root.first.as<std::string>(), root.second);
    BOOST_ASSERT(
      spec.name == root.first.as<std::string>()
    );

    // Verifies the dimensions of the topology.
    std::vector<std::string> dims = root.second["dims"].as<std::vector<std::string>>();
    BOOST_ASSERT(
      dims == spec.dims
    );

    // Computes the pw_multi_aff for the NoC cost.
    isl::pw_multi_aff aff = isl::pw_multi_aff(GetIslCtx(), root.second["affine"].as<std::string>());
    BOOST_ASSERT(
      isl_pw_multi_aff_is_equal(
        spec.noc_cost.copy(), aff.copy()
      )
    );

    // Verifies the domain restrictions.
    isl::set input_domain  = spec.domain.product(spec.domain);
    isl::pw_multi_aff restricted = spec.noc_cost.intersect_domain(input_domain);
    isl::pw_multi_aff map = isl::pw_multi_aff(GetIslCtx(), root.second["domain_restriction"].as<std::string>());
    BOOST_ASSERT(
      isl_pw_multi_aff_is_equal(
        restricted.copy(), map.copy()
      )
    );
  }
}