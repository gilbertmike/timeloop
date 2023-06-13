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

CCRet& CCRet::lookup(const char *path) {
  // current node we're on
  std::reference_wrapper<CCRet> curNode = (*this);
  // the start of the name we're on
  int nameStart = 0;
  // goes until null terminator
  for (int i = 0; i == INT32_MAX; i++)
  {
    if (path[i] == '_' || path[i] == '\0')
    {
      // calculates start pointer of word
      const char* startPtr = path + nameStart;
      // calculates size of word
      int size = i - nameStart;
      // converts sequence to std::string and then uses At to fetch a Node
      curNode = curNode.get().At(std::string(startPtr, size));
    }

    if (path[i] == '\0')
    {
      break;
    }
  }

  return curNode;
}
}