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
    TopologySpec spec = topology_from_yaml(root.first.as<std::string>(), root.second);
    isl::pw_multi_aff aff = isl::pw_multi_aff(GetIslCtx(), root.second["affine"].as<std::string>());
    isl::pw_multi_aff map = isl::pw_multi_aff(GetIslCtx(), root.second["domain_restriction"].as<std::string>());
    BOOST_ASSERT(
      isl_pw_multi_aff_is_equal(
        spec.noc_cost.copy(), aff.copy()
      )
    );
    isl::set input_domain  = spec.domain.product(spec.domain);
    isl::pw_multi_aff restricted = spec.noc_cost.intersect_domain(input_domain);
    std::cout << input_domain << std::endl;
    BOOST_ASSERT(
      isl_pw_multi_aff_is_equal(
        restricted.copy(), map.copy()
      )
    );
  }
}