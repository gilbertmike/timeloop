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
    isl::pw_multi_aff exp = isl::pw_multi_aff(GetIslCtx(), root.second["equiv"].as<std::string>());
    BOOST_ASSERT(
        spec.noc_cost.plain_is_equal(exp)
    );
  }
}