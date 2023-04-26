#include <iostream>
using namespace std;

/**
 * Fused matrix multiplication
 *  Z1[m1, n1] = A1[m1, k1] B1[k1, n1]
 *  Z2[m2, n2] = Z1[m2, k2] B2[k2, n2]
 *  m1 = m2, n1 = k2
 */

int num_tiles(int dim) {
    int num_tiles = 0;
    for(int i = 1; i <= dim; i++) {
        if(dim % i == 0) num_tiles++;
    }
    return num_tiles;
}
int num_mappings(int m1, int k1, int n1, int n2) {
    // storage can be between DRAM and pipeline
    // ordering of the loops
    return num_tiles(m1) * num_tiles(n1) * num_tiles(n2);  // tile only on consumer: m2 = m1, k2 = n1, n2
}

int main() {
    cout << num_mappings(4, 3, 4, 4) << endl;
    cout << num_mappings(8, 8, 4, 4) << endl;
    cout << num_mappings(100000, 100000, 100000, 10000) << endl;
    return 0;
}