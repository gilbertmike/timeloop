#include <ctime>
#include <iostream>
#include <string>

#include <boost/test/unit_test.hpp>

#include "isl-wrapper/ctx-manager.hpp"
#include "loop-analysis/spatial-analysis.hpp"

#define DPRINT(x) if (std::string(std::getenv("DEBUG")) == "1") std::cout << #x << ": " << x << std::endl;

BOOST_AUTO_TEST_CASE(TestSimpleMulticastModel_0)
{
  using namespace analysis;

  auto fill = Fill(
    {Temporal(), Spatial(0, 0), Spatial(1, 0)},
    isl::map(
      GetIslCtx(),
      "{ [t, x, y] -> [t + x + y] : 0 <= x < 2 and 0 <= y < 2 and 0 <= t < 4 }"
    )
  );

  auto occ = Occupancy(fill.dim_in_tags, fill.map);

  auto multicast_model = SimpleMulticastModel(false);

  auto info = multicast_model.Apply(0, fill, occ);

  BOOST_CHECK(info.fulfilled_fill.map.is_equal(
    isl::map(
      GetIslCtx(),
      "{ [t, x, y] -> [t + x + y] : 0 <= t < 4 and 0 <= x < 2 and 0 <= y < 2 }"
    )
  ));

  BOOST_CHECK(info.parent_reads.map.is_equal(
    isl::map(
      GetIslCtx(),
      "{ [t] -> [d] : 0 <= t < 4 and t <= d < t + 3 }"
    )
  ));

  BOOST_CHECK(info.compat_access_stats.size() == 1);
  for (const auto& [multicast_scatter, stats] : info.compat_access_stats)
  {
    auto [multicast, scatter] = multicast_scatter;

    BOOST_CHECK(multicast == 1);
    BOOST_CHECK(scatter == 1);
    BOOST_CHECK(stats.accesses == 12);
    BOOST_CHECK(stats.hops == 0);
  }

  multicast_model = SimpleMulticastModel(true);

  info = multicast_model.Apply(0, fill, occ);

  BOOST_CHECK(info.fulfilled_fill.map.is_equal(
    isl::map(
      GetIslCtx(),
      "{ [t, x, y] -> [t + x + y] : 0 <= t < 4 and 0 <= x < 2 and 0 <= y < 2 }"
    )
  ));

  BOOST_CHECK(info.parent_reads.map.is_equal(
    isl::map(
      GetIslCtx(),
      "{ [t] -> [d] : 0 <= t < 4 and t <= d < t + 3 }"
    )
  ));

  BOOST_CHECK(info.compat_access_stats.size() == 1);
  for (const auto& [multicast_scatter, stats] : info.compat_access_stats)
  {
    auto [multicast, scatter] = multicast_scatter;

    BOOST_CHECK(multicast == 1);
    BOOST_CHECK(scatter == 1);
    BOOST_CHECK(stats.accesses == 12);
    BOOST_TEST(stats.hops == 3.667, boost::test_tools::tolerance(0.001));
  }
}

BOOST_AUTO_TEST_CASE(TestSimpleMulticastModel_SpatialPC)
{
  using namespace analysis;

  auto fill = Fill(
    {Temporal(), Spatial(0, 0), Spatial(1, 0)},
    isl::map(
      GetIslCtx(),
      "{ [t, x, y] -> [d, y] : 0 <= x < 4 and 0 <= y < 2 and 0 <= t < 4 and x <= d < x+2 }"
    )
  );

  auto occ = Occupancy(fill.dim_in_tags, fill.map);

  auto multicast_model = SimpleMulticastModel(true);

  auto info = multicast_model.Apply(0, fill, occ);

  BOOST_CHECK(info.fulfilled_fill.map.is_equal(
    isl::map(
      GetIslCtx(),
      "{ [t, x, y] -> [d, y] : 0 <= x < 4 and 0 <= y < 2 and 0 <= t < 4 and x <= d < x+2 }"
    )
  ));

  BOOST_CHECK(info.parent_reads.map.is_equal(
    isl::map(
      GetIslCtx(),
      "{ [t] -> [d, y] : 0 <= y < 2 and 0 <= t < 4 and 0 <= d < 5 }"
    )
  ));

  BOOST_CHECK(info.compat_access_stats.size() == 1);
  for (const auto& [multicast_scatter, stats] : info.compat_access_stats)
  {
    auto [multicast, scatter] = multicast_scatter;

    BOOST_CHECK(multicast == 1);
    BOOST_CHECK(scatter == 1);
    BOOST_CHECK(stats.accesses == 40);
    BOOST_TEST(stats.hops == 5.2, boost::test_tools::tolerance(0.001));
  }
}


std::vector<analysis::SpaceTime> construct_space_time(const YAML::Node &dims)
{
  std::vector<analysis::SpaceTime> space_time;
  for (YAML::Node dim : dims)
  {
    if (dim["type"].as<std::string>() == "Temporal")
    {
      space_time.push_back(analysis::Temporal());
    }
    else if (dim["type"].as<std::string>() == "Spatial")
    {
      space_time.push_back(
        analysis::Spatial(dim["spatial_dim"].as<int>(), dim["target"].as<int>())
      );
    }
  }
  return space_time;
}


BOOST_AUTO_TEST_CASE(TestDistributedMulticastHyperCubeModel)
{
  using namespace analysis;

  std::string TEST_CASES_FILE = "./src/unit-test/multicast/test_cases.yaml";
  YAML::Node test_cases = YAML::LoadFile(TEST_CASES_FILE);
  
  std::cout << "Running DistributedMulticastHyperCubeModel Test" << std::endl;
  for (auto test : test_cases.as<std::vector<YAML::Node>>()) {
    ///@brief Reads test case parameters.
    int buf_id = 0;
    std::string fill_str = test["fill"].as<std::string>();
    std::string occ_str = test["occ"].as<std::string>();
    std::string dist_func_str = test["dist_func"].as<std::string>();
    ///@brief Construct the necessary Timeloop objects.
    auto dims = construct_space_time(test["dims"]);
    DPRINT(fill_str);
    Fill fill = Fill(
      dims,
      isl::map(GetIslCtx(), fill_str)
    );
    DPRINT(occ_str);
    Occupancy occ = Occupancy(
      dims,
      isl::map(GetIslCtx(), occ_str)
    );
    DPRINT(dist_func_str);
    auto multicast_model = DistributedMulticastHypercubeModel(
      true, isl::map(GetIslCtx(), dist_func_str)
    );
    ///@brief Apply the model
    TransferInfo info = multicast_model.Apply(buf_id, fill, occ);
    ///@brief Check the results
    isl::val sum_extract = isl::manage(isl_pw_qpolynomial_eval(
      info.p_hops, isl_point_zero(isl_pw_qpolynomial_get_domain_space(info.p_hops))
    ));
    long ret = sum_extract.get_num_si();
    std::cout << "Returned Value: " << ret << std::endl;
    
    ///@note The if block is used for debugging test cases not yet implemented.
    if (test["expected"]["hypercube_hops"].IsNull()) {
      std::cout << "Test case in progress" << std::endl;
      std::cout << "Fill: " << fill_str << std::endl;
      std::cout << "Occupancy: " << occ_str << std::endl;
      std::cout << "Dist Func: " << dist_func_str << std::endl;
    } else {
      BOOST_CHECK(ret == test["expected"]["hypercube_hops"].as<long>());
    }
  }
  std::cout << "DistributedMulticastHyperCubeModel Test Passed" << std::endl;
}
/*
BOOST_AUTO_TEST_CASE(TestDistributedMulticastOrderedExtentsDORModel)
{
  using namespace analysis;

  std::string TEST_CASES_FILE = "./src/unit-test/multicast/test_cases.yaml";
  YAML::Node test_cases = YAML::LoadFile(TEST_CASES_FILE);
  
  std::cout << "Running DistributedMulticastOrderedExtentsDORModel Test" << std::endl;
  for (auto test : test_cases.as<std::vector<YAML::Node>>()) {
    // Read test case parameters
    int buf_id = 0;
    auto dims = construct_space_time(test["dims"]);
    Fill fill = Fill(
      dims,
      isl::map(GetIslCtx(), test["fill"].as<std::string>())
    );
    Occupancy occ = Occupancy(
      dims,
      isl::map(GetIslCtx(), test["occ"].as<std::string>())
    );
    // Apply the model
    auto multicast_model = DistributedMulticastOrderedExtentsDORModel(
      true, isl::map(GetIslCtx(), test["dist_func"].as<std::string>())
    );
    TransferInfo info = multicast_model.Apply(buf_id, fill, occ);
    // Check the results
    isl::val sum_extract = isl::manage(isl_pw_qpolynomial_eval(
      info.p_hops, isl_point_zero(isl_pw_qpolynomial_get_domain_space(info.p_hops))
    ));
    long ret = sum_extract.get_num_si();
    
    std::cout << "Returned Value: " << ret << std::endl;
    BOOST_CHECK(ret == test["expected"]["extent_DOR_hops"].as<long>());
  }
  std::cout << "DistributedMulticastOrderedExtentsDORModel Test Passed" << std::endl;
}
*/