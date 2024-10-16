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

BOOST_AUTO_TEST_CASE(TestDistributedMulticast_Model)
{
  using namespace analysis;

  int M_int = 1024;
  int N_int = 1024;
  std::string M = std::to_string(M_int);
  std::string N = std::to_string(N_int);
  std::vector<int> D_vals({1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024});
  isl_ctx *p_ctx = isl_ctx_alloc();

  for (int D_int : D_vals) {
    std::cout << "D: " << D_int << std::endl;
    std::string D = std::to_string(D_int);
    // Defines the dst fill map as a string.
    std::string dst_fill =  "{dst[xd, yd] -> data[a, b] : b=yd and 0 <= xd < "+M+
                            " and 0 <= yd < "+N+" and 0 <= a < "+M+" and 0 <= b < "+N+" }";
    auto fill = Fill(
      {Spatial(0, 0), Spatial(1, 0)},
      isl::map(
        GetIslCtx(),
        dst_fill
      )
    );
    // Defines the src occupancy map as a string.
    std::string src_occupancy = "{src[xs, ys] -> data[a, b] : ("+D+"*xs)%"+M+" <= a <= ("+
                                D+"*xs+"+D+"-1)%"+M+" and b=ys and 0 <= xs < "+M+
                                " and 0 <= ys < "+N+" and 0 <= a < "+M+" and 0 <= b < "+N+" }";
    auto occ = Occupancy(fill.dim_in_tags, 
      isl::map(
        GetIslCtx(),
        src_occupancy
    ));

    auto multicast_model = DistributedMulticastModel(true);
    auto info = multicast_model.Apply(0, fill, occ);

    // BOOST_CHECK(info.fulfilled_fill.map.is_equal(
    //   isl::map(
    //     GetIslCtx(),
    //     nullptr
    //   )
    // ));

    // BOOST_CHECK(info.parent_reads.map.is_equal(
    //   isl::map(
    //     GetIslCtx(),
    //     nullptr
    //   )
    // ));

    // BOOST_CHECK(info.compat_access_stats.size() == 1);
    // for (const auto& [multicast_scatter, stats] : info.compat_access_stats)
    // {
    //   auto [multicast, scatter] = multicast_scatter;

    //   BOOST_CHECK(multicast == 1);
    //   BOOST_CHECK(scatter == 1);
    //   BOOST_CHECK(stats.accesses == 40);
    //   BOOST_TEST(stats.hops == 5.2, boost::test_tools::tolerance(0.001));
    // }
  }
}