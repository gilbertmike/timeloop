#include "distributed-buffers.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#pragma O2
/// NOTES FOR NON-TREE MULTICAST SCENARIO
// - Load balancing issues for multiple minimally distant sources.
// - Compose minimal distances with the other set to remove non-minimal pairs then
// move on with the rest of the algorithm.
__isl_give isl_map *identify_mesh_casts( 
    __isl_take isl_map *src_occupancy, 
    __isl_take isl_map *dst_fill, 
    __isl_take isl_map *dist_func
) {
    /* Makes [[dst -> data] -> dst] -> [data] */
    isl_set *wrapped_dst_fill = isl_map_wrap(dst_fill);
    isl_map *wrapped_fill_identity =
        isl_map_identity(isl_space_map_from_set(isl_set_get_space(
            wrapped_dst_fill
        )));
    DUMP(wrapped_fill_identity);
    wrapped_fill_identity = isl_map_intersect_domain(
        wrapped_fill_identity,
        wrapped_dst_fill
    );
    DUMP(wrapped_fill_identity);
    /* Makes [dst -> data] -> [dst -> data] */
    isl_map *uncurried_fill_identity = isl_map_uncurry(wrapped_fill_identity);
    DUMP(uncurried_fill_identity);
    /* Inverts src_occupancy such that data implies source.
    * i.e. {[xs, ys] -> [d0, d1]} becomes {[d0, d1] -> [xs, ys]} */
    isl_map *src_occupancy_inverted = isl_map_reverse(src_occupancy);

    isl_map *dst_to_data_to_dst_TO_src = isl_map_apply_range(
        uncurried_fill_identity,
        src_occupancy_inverted
    );
    DUMP(dst_to_data_to_dst_TO_src);

    isl_map *dst_to_data_TO_dst_to_src =
        isl_map_curry(dst_to_data_to_dst_TO_src);
    DUMP(dst_to_data_TO_dst_to_src);

    // Calculates the distance of all the dst-src pairs with matching data.
    DUMP(dist_func);
    isl_map *distances_map = isl_map_apply_range(
        isl_map_copy(dst_to_data_TO_dst_to_src), isl_map_copy(dist_func)
    );
    DUMP(distances_map);
    isl_map *dst_to_data_TO_dst_to_src_TO2_dst_to_src = isl_map_range_map(dst_to_data_TO_dst_to_src);
    isl_map *dst_to_data_TO_dst_to_src_TO2_dist = isl_map_apply_range(dst_to_data_TO_dst_to_src_TO2_dst_to_src, dist_func);
    DUMP(dst_to_data_TO_dst_to_src_TO2_dist);

    // Gets the minimal distance pairs.
    isl_map *lexmin_distances = isl_map_lexmin(distances_map);
    isl_map *assoc_dist_with_src = isl_map_apply_range(lexmin_distances, isl_map_reverse(
        dst_to_data_TO_dst_to_src_TO2_dist
    ));
    // Isolates the relevant minimal pairs.
    isl_map *minimal_pairs = isl_set_unwrap(isl_map_range(assoc_dist_with_src));
    DUMP(minimal_pairs);
    // Isolates the multicast networks.
    isl_map *multicast_networks = isl_map_curry(minimal_pairs);
    multicast_networks = isl_set_unwrap(isl_map_range(multicast_networks));
    DUMP(multicast_networks);
    multicast_networks = isl_map_uncurry(multicast_networks);
    DUMP(multicast_networks);
    multicast_networks = isl_map_lexmin(multicast_networks);
    DUMP(multicast_networks);
    multicast_networks = isl_map_curry(multicast_networks);
    DUMP(multicast_networks);

    return multicast_networks;
}   

__isl_give isl_map *identify_mesh_casts(
    isl_ctx *const p_ctx,
    const std::string& src_occupancy, 
    const std::string& dst_fill, 
    const std::string& dist_func
) {
    // Reads the string representations of the maps into isl objects.
    isl_map *p_src_occupancy = isl_map_read_from_str(
        p_ctx,
        src_occupancy.c_str()
    );
    isl_map *p_dst_fill = isl_map_read_from_str(
        p_ctx,
        dst_fill.c_str()
    );
    isl_map *p_dist_func = isl_map_read_from_str(
        p_ctx,
        dist_func.c_str()
    );

    // Calls the isl version of analyze_latency.
    isl_map *ret = identify_mesh_casts(
        p_src_occupancy,
        p_dst_fill,
        p_dist_func
    );

    return ret;
}

long cost_mesh_cast(
    __isl_take isl_map *mesh_cast_networks,
    __isl_take isl_map *dist_func
) {
    DUMP(mesh_cast_networks);
    DUMP(dist_func);
    
    /**
     * Makes mesh_cash_networks from [a, b] -> [[xd, yd] -> [xs -> ys]] to 
     * [[a, b] -> [xs, ys]] -> [xd, yd]
     */
    mesh_cast_networks = isl_map_range_reverse(mesh_cast_networks);
    DUMP(mesh_cast_networks);
    // Uncurrys the mesh_cast_networks to [[a, b] -> [xs, ys]] -> [xd, yd]
    mesh_cast_networks = isl_map_uncurry(mesh_cast_networks);
    DUMP(mesh_cast_networks);
    
    // Projects away the xd dimension from mesh_cast_networks.
    isl_map *multicast_simplification = isl_map_project_out(mesh_cast_networks, isl_dim_out, 1, 1);
    DUMP(multicast_simplification);
    // Finds max(yd) - min(yd) for each [a, b] -> [xs, ys].
    isl_map *multicast_max = isl_map_lexmax(isl_map_copy(multicast_simplification));
    DUMP(multicast_max);
    isl_map *multicast_min = isl_map_lexmin(multicast_simplification);
    DUMP(multicast_min);
    // Subtracts the max from the min to get the range.
    isl_map *multicast_min_neg = isl_map_neg(multicast_min);
    DUMP(multicast_min_neg);
    isl_map *multi_cast_cost = isl_map_sum(multicast_max, multicast_min);
    DUMP(multi_cast_cost);

    // Converts to a qpolynomial for addition over range.
    isl_multi_pw_aff *dirty_distances_aff =isl_multi_pw_aff_from_pw_multi_aff(
        isl_pw_multi_aff_from_map(multi_cast_cost)
    );
    DUMP(dirty_distances_aff);
    assert(isl_multi_pw_aff_size(dirty_distances_aff) == 1);
    isl_pw_aff *distances_aff = isl_multi_pw_aff_get_at(dirty_distances_aff, 0);
    DUMP(distances_aff);
    isl_multi_pw_aff_free(dirty_distances_aff);
    auto *dirty_distances_fold = isl_pw_qpolynomial_from_pw_aff(distances_aff);
    DUMP(dirty_distances_fold);

    // Does the addition over range.
    isl_pw_qpolynomial *sum = isl_pw_qpolynomial_sum(isl_pw_qpolynomial_sum(dirty_distances_fold));
    // Grabs the return value as an isl_val.
    isl_val *sum_extract = isl_pw_qpolynomial_eval(sum, isl_point_zero(isl_pw_qpolynomial_get_domain_space(sum)));
    long ret = isl_val_get_num_si(sum_extract);

    return ret;
}

/**
 * @return The cost per datum of each network.
 */
long cost_mesh_cast(
    isl_ctx *const p_ctx,
    const std::string& mesh_cast_networks,
    const std::string& dist_func
) {

    // Reads the string representations of the maps into isl objects.
    isl_map *p_mesh_cast_networks = isl_map_read_from_str(
        p_ctx,
        mesh_cast_networks.c_str()
    );
    isl_map *p_dist_func = isl_map_read_from_str(
        p_ctx,
        dist_func.c_str()
    );

    // Calls the isl version of analyze_latency.
    auto ret = cost_mesh_cast(
        p_mesh_cast_networks,
        p_dist_func
    );

    return ret;
}

int main(int argc, char* argv[])
{
    int M_int = 1024;
    int N_int = 1024;
    std::string M = std::to_string(M_int);
    std::string N = std::to_string(N_int);
    std::vector<int> D_vals({1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024});
    clock_t start, end;
    double cpu_time_used;
    isl_ctx *p_ctx = isl_ctx_alloc();

    for (int D_int : D_vals) {
        start = clock();
        std::string D = std::to_string(D_int);
        // Defines the src occupancy map as a string.
        std::string src_occupancy = "{src[xs, ys] -> data[a, b] : ("+D+"*xs)%"+M+" <= a <= ("+
                                    D+"*xs+"+D+"-1)%"+M+" and b=ys and 0 <= xs < "+M+
                                    " and 0 <= ys < "+N+" and 0 <= a < "+M+" and 0 <= b < "+N+" }";
        // Defines the dst fill map as a string.
        std::string dst_fill =  "{dst[xd, yd] -> data[a, b] : b=yd and 0 <= xd < "+M+
                                " and 0 <= yd < "+N+" and 0 <= a < "+M+" and 0 <= b < "+N+" }";

        // Defines the distance function string.
        std::string dist_func_str = R"DIST({
            [dst[xd, yd] -> src[xs, ys]] -> dist[(xd - xs) + (yd - ys)] : 
                xd >= xs and yd >= ys;
            [dst[xd, yd] -> src[xs, ys]] -> dist[-(xd - xs) + -(yd - ys)] : 
                xd < xs and yd < ys;
            [dst[xd, yd] -> src[xs, ys]] -> dist[-(xd - xs) + (yd - ys)] : 
                xd < xs and yd >= ys;
            [dst[xd, yd] -> src[xs, ys]] -> dist[(xd - xs) + -(yd - ys)] : 
                xd >= xs and yd < ys
            })DIST";
    
        isl_map *mcs = identify_mesh_casts(p_ctx, src_occupancy, dst_fill, dist_func_str);
        DUMP(mcs);
        long res = cost_mesh_cast(p_ctx, isl_map_to_str(mcs), dist_func_str);
        std::cout << res << std::endl;
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        // std::cout << "Time: " << cpu_time_used << std::endl;
    }
}