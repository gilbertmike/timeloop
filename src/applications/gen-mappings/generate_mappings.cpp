#include <iostream>
#include "../include/mapping/fused-mapping.hpp"
#include <algorithm>

using namespace mapping;

/**
 * Fused matrix multiplication
 *  Z1[m1, n1] = A1[m1, k1] B1[k1, n1]
 *  Z2[m2, n2] = Z1[m2, k2] B2[k2, n2]
 */

using DataSpaceID = problem::Shape::DataSpaceID;
const DataSpaceID A1_id = 0;

std::vector<std::vector<int> > getForLoopOrderings(std::vector<int> dims) {
    std::vector<std::vector<int> > perms;
    do {
        std::vector<int> this_dim = dims;
        perms.push_back(this_dim);
        for(auto dim : dims) {
            std::cout << dim << ' ';
        }
        std::cout << std::endl;
    } while(std::next_permutation(dims.begin(), dims.end()));
    return perms;
}

// TODO: 5/3 embed workload parsing
int main(/*int argc, char* argv[]*/) {
    std::vector<FusedMapping> mappings;
    int m2_id = 0, n2_id = 1, k2_id  = 2; // all permutations
    std::vector<int> dims = {m2_id, n2_id, k2_id};
    std::vector<std::vector<int> > orderings;
    do {
        // std::cout << "Storage DRAM ID: " << mapping.AddChild<Storage>(1, 0, 0) << std::endl; // DRAM
        // mapping.AddChild<For>(2, "it1", 0); // 0 is dim_id
        // mapping.AddChild<Storage>(3, 0, 0); // storage for Z1
        
        std::vector<int> this_dim = dims;
        orderings.push_back(this_dim);
        for(auto dim : dims) {
            std::cout << dim << ' ';
        }
        std::cout << std::endl;
        for(size_t i = 0; i < dims.size(); i++) {
            FusedMapping mapping;
            for(size_t j = 0; j < dims.size(); j++) {
                if(j == i) {
                    mapping.AddChild<Storage>(j+1, 0, 0);
                    continue;
                }
                int p_id = (int) (j < i ? j : j + 1);
                mapping.AddChild<For>(p_id, "it" + (char)('0' + p_id), dims[j]);
            }
            mappings.push_back(mapping);
        }
    } while(std::next_permutation(dims.begin(), dims.end()));
    std::cout << mappings.size() << std::endl;
}