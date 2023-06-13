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
  assert(isList());

  if (idx < Size()) return At((YAMLVector::size_type) idx);
  else return CCRet();
}
}