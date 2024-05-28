#include <boost/test/unit_test.hpp>

#include "binding/binding.hpp"

BOOST_AUTO_TEST_CASE(TestBinding)
{
  auto spec = BindingSpec{
    .virtual_buffer = 0,
    .physical_buffer = 0,
    .descriptors = {}
  };

  spec.descriptors.push_back(
    SpatialDuplicate{.num_duplication=2, .spatial_dimension=0}
  );

  std::cout << spec << std::endl;
}