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
    isl::pw_aff noc = noc_from_yaml(root.first.as<std::string>(), root.second);
    std::cout << "Computed: " << noc << std::endl;
    isl::pw_aff exp = isl::pw_aff(GetIslCtx(), root.second["equiv"].as<std::string>());
    std::cout << "Expected: " << exp << std::endl;
    BOOST_ASSERT(
        isl_pw_aff_is_equal(noc.get(), exp.get())
    );
  }
}