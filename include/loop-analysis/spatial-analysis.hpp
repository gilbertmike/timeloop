#pragma once
#include "loop-analysis/isl-ir.hpp"
#include "isl/polynomial.h"

/******************************************************************************
 * Interface
 *****************************************************************************/

namespace analysis
{

struct LinkTransferInfo
{
  Transfers link_transfer;
  Fill unfulfilled_fill;
};

struct LinkTransferModel
{
  virtual LinkTransferInfo Apply(const Fill&, const Occupancy&) const = 0;
};

struct MulticastInfo
{
  isl::map reads;
  isl_pw_qpolynomial* p_hops;

  /***************** Compatibility with Timeloop v2.0 ************************/
  struct AccessStats
  {
    double accesses;
    double hops;
  };
  std::map<std::pair<uint64_t, uint64_t>, AccessStats> compat_access_stats;
  /***************************************************************************/

  ~MulticastInfo();
};

struct MulticastModel
{
  virtual MulticastInfo Apply(const Fill&) const = 0;
};

struct SpatialReuseInfo
{
  LinkTransferInfo link_transfer_info;
  MulticastInfo multicast_info;
};

struct SpatialReuseAnalysisInput
{
  const LogicalBuffer& buf;
  const Fill& children_fill;
  const Occupancy& children_occupancy;

  SpatialReuseAnalysisInput(const LogicalBuffer&,
                            const Fill&,
                            const Occupancy&);
};

struct SpatialReuseModels
{
  const LinkTransferModel& link_transfer_model;
  const MulticastModel& multicast_model;

  SpatialReuseModels(const LinkTransferModel&, const MulticastModel&);
};

SpatialReuseInfo SpatialReuseAnalysis(SpatialReuseAnalysisInput input,
                                      SpatialReuseModels models);

/******************************************************************************
 * Model Classes
 *****************************************************************************/

/**
 * @brief A link transfer model for 2-dimensional mesh interconnect.
 */
class SimpleLinkTransferModel : public LinkTransferModel
{
 public:
  SimpleLinkTransferModel();

  LinkTransferInfo
  Apply(const Fill& fills, const Occupancy& occupancies) const;

 private:
  isl::map connectivity_;
};

/**
 * @brief A multicast model for 2-dimensional array.
 *
 * @note Differs from Timeloop's original multicast model in terms of partial
 *   tile overlap multicast. This model assumes partial overlaps can benefit
 *   from multicasting. The original model does the opposite.
 * @note Does not directly model distributed multicast. Uses the same methods
 *   as the original multicast model (i.e., assumes parents have equally
 *   distributed tiles)
 */
class SimpleMulticastModel : public MulticastModel
{
 public:
  SimpleMulticastModel();

  MulticastInfo Apply(const Fill& fills) const;

 private:
  isl::map connectivity_;
};
}