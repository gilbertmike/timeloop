// represents YAML maps
#include <map>
// for type-safe unions of YAML types
#include <variant>
// for YAML arrays
#include <vector>
// for std errors from variants
#include <stdexcept>
#include <memory>
#include <iostream>

#include "compound-config/compound-config-struct.hpp"

namespace structured_config {
CCRet CCRet::operator [](int idx) const {
  if (idx != 0) return CCRet();
  else return CCRet();
  // assert(isList());

  // return At((YAMLVector::size_type) idx);
}
}