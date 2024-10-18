#include <ctime>
#include <iostream>
#include <string>

#include <boost/test/unit_test.hpp>

#include "isl-wrapper/ctx-manager.hpp"
#include "loop-analysis/spatial-analysis.hpp"

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
      throw std::runtime_error("Temporal tag not supported");
    }
    else if (dim["type"].as<std::string>() == "Spatial")
    {
      space_time.push_back(analysis::Spatial(dim["spatial_dim"].as<int>(), dim["target"].as<int>()));
    }
  }
  return space_time;
}


BOOST_AUTO_TEST_CASE(TestDistributedMulticast_Model)
{
  using namespace analysis;

  std::string TEST_CASES_FILE = "./src/unit-test/multicast/test_cases.yaml";
  YAML::Node test_cases = YAML::LoadFile(TEST_CASES_FILE);
  
  for (auto test : test_cases) {
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
    isl::map dist_func = isl::map(GetIslCtx(), test["dist_func"].as<std::string>());
    auto multicast_model = DistributedMulticastModel(true);
    TransferInfo info = multicast_model.Apply(buf_id, fill, occ);
    isl_val *sum_extract = isl_pw_qpolynomial_eval(info.p_hops, isl_point_zero(isl_pw_qpolynomial_get_domain_space(info.p_hops)));
    long ret = isl_val_get_num_si(sum_extract);
    std::cout << "Number of hops: " << ret << std::endl;
  }
}