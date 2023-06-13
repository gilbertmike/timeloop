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
CCRet& CCRet::operator [](int idx) {
  assert(isList());

  return At(idx);
}

bool CCRet::exists(std::string name) const {
  return std::get<YAMLMap>(data_).count(name) != 0;
}
}