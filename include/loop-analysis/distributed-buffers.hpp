#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

// Includes ISL affine list/piecewise functions.
#include <isl/aff.h>
#include <isl/ilp.h>
// Includes isl qpolynomials.
#include <isl/polynomial.h>
// Includes ISL maps/binary relations.
#include <isl/map.h>
// Includes ISL ids and dspaces.
#include <isl/id.h>
#include <isl/space.h>
// Imports ISL sets.
#include <isl/set.h>
// Imports ISL val.
#include <isl/val.h>
#include <barvinok/isl.h>
#include <barvinok/polylib.h>

__isl_give isl_pw_qpolynomial* gather_pw_qpolynomial_from_fold(__isl_take isl_pw_qpolynomial_fold* pwqpf);

long analyze_jumps(isl_map *p_src_occupancy, isl_map *p_dst_fill, isl_pw_aff *dist_func);
long analyze_jumps(const std::string& src_occupancy, const std::string& dst_fill, const std::string& dist_func);
long analyze_latency(isl_map *p_src_occupancy, isl_map *p_dst_fill, isl_pw_aff *dist_func);
long analyze_latency(const std::string& src_occupancy, const std::string& dst_fill, const std::string& dist_func);
std::string nd_manhattan_metric(std::vector<std::string> src_dims, std::vector<std::string> dst_dims);
std::string n_long_ring_metric(long n);

// Defines debug variables from environment variables.
#include <string.h>
bool islIntermediates = (getenv("ISL_INTERMEDIATES") != NULL) &&
                        (strcmp(getenv("ISL_INTERMEDIATES"), "0") != 0);

// Defines a function to programatically generate an n-dimensional Manhattan distance function.
std::string nd_manhattan_metric(std::vector<std::string> src_dims, std::vector<std::string> dst_dims);
// Defines a function to programatically generate an n-circumference ring distance function.
std::string n_long_ring_metric(long ring_circumference);

// Defines debug dump function.
void dump(const std::string& str, isl_map *map)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_map_dump(map);
    }
}

void dump(const std::string& str, isl_pw_aff *pw_aff)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_pw_aff_dump(pw_aff);
    }
}

void dump(const std::string& str, isl_multi_pw_aff *multi_pw_aff)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_multi_pw_aff_dump(multi_pw_aff);
    }
}

void dump(const std::string& str, isl_multi_val *multi_val)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_multi_val_dump(multi_val);
    }
}

void dump(const std::string& str, isl_set *set)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_set_dump(set);
    }
}

void dump(const std::string& str, isl_pw_qpolynomial* pwqp) {
    if (islIntermediates) {
        std::cout << str << std::endl;
        isl_pw_qpolynomial_dump(pwqp);
    }
}

#define DUMP(varname) dump(#varname, varname)

/// @brief Strings representing the src and dst datum holds/requests in ISL.
struct binding_struct
{
    const std::string srcs;
    const std::string dsts;
};
typedef std::shared_ptr<binding_struct> binding;
/// @brief Represents the cost of abstracting a binding to a higher level.
struct abstraction_cost_struct
{
    const long fold_cost;
    const long multicast_cost;
};
typedef std::shared_ptr<abstraction_cost_struct> abstraction_cost;
/// @brief Struct with the abstracted binding + the cost of abstracting it.
struct abstracted_binding_struct
{
    const binding abstraction;
    const abstraction_cost cost;
};
typedef std::unique_ptr<abstracted_binding_struct> abstracted_binding;
/** 
 * @brief Defines the struct that comprises the result of folding and the unique
 * ptr to it that represents what is returned by fold.
 */
struct fold_struct
{
    const long cost;
    const std::string folded_repr;
};
typedef std::unique_ptr<fold_struct> fold_result;
/// @brief Defines the struct characterizing the collapsing behavior of a layer.
struct collapse_struct
{
    const std::string src_collapser;
    const std::string dst_collapser;
};
typedef std::shared_ptr<collapse_struct> collapse;