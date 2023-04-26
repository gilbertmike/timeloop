#include <iostream>
#include "../include/mapping/fused-mapping.hpp"

using namespace mapping;

/**
 * Fused matrix multiplication
 *  Z1[m1, n1] = A1[m1, k1] B1[k1, n1]
 *  Z2[m2, n2] = Z1[m2, k2] B2[k2, n2]
 */

using DataSpaceID = problem::Shape::DataSpaceID;
const DataSpaceID A1_id = 0;

int main(/*int argc, char* argv[]*/) {
    std::vector<FusedMapping> mappings;
    while(true) {
        FusedMapping mapping;
        std::cout << mapping.AddChild<Storage>(1, 0, 0) << std::endl;
        break;
    }
}