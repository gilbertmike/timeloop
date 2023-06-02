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
    int main(int argc, char* argv[]) {
        auto val = CCRet::Vector();
        val.EmplaceBack(CCRet::Literal("a value"));
        std::cout << val.Size() << std::endl;
    }
}