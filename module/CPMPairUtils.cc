#include "CPMPairUtils.h"

#include "CPMHelixPoCA.h"
#include "CPMLocalLinePoCA.h"

#include <cmath>

bool cpmPairHasGoodPt(
    const int charge,
    const double pt,
    const double min_pt)
{
  return charge != 0 &&
         std::isfinite(pt) &&
         pt > 0.0 &&
         pt >= min_pt;
}

double cpmPairWeightFromPt(const double pt_a, const double pt_b)
{
  return 1.0 / (pt_a * pt_b);
}

double cpmWrapDeltaPhi(double value)
{
  constexpr double pi = 3.141592653589793238462643383279502884;
  while (value > pi)
  {
    value -= 2.0 * pi;
  }
  while (value <= -pi)
  {
    value += 2.0 * pi;
  }
  return value;
}

CPMPairInput makeCPMPairInput(const TrackStateRecord& record)
{
  return {
      record.track.charge,
      record.track.pt,
      record.state.position,
      record.state.momentum,
      record.cluster.cluster_minus_voxel_center};
}

CPMPairResult computeCPMPair(
    const CPMPairInput& input_a,
    const CPMPairInput& input_b,
    const Vector3& voxel_center,
    const CPMPairOptions& options)
{
  CPMPairResult result;

  if (!cpmPairHasGoodPt(input_a.charge, input_a.pt, options.min_pt) ||
      !cpmPairHasGoodPt(input_b.charge, input_b.pt, options.min_pt))
  {
    result.status = CPMPairStatus::PtRejected;
    return result;
  }

  result.pair_weight = cpmPairWeightFromPt(input_a.pt, input_b.pt);
  if (!(std::isfinite(result.pair_weight) && result.pair_weight > 0.0))
  {
    result.status = CPMPairStatus::InvalidWeight;
    return result;
  }

  const Vector3 point_a =
      input_a.state_position - input_a.cluster_minus_voxel_center;
  const Vector3 point_b =
      input_b.state_position - input_b.cluster_minus_voxel_center;

  bool poca_valid = false;
  if (options.solver == CPMPairSolver::Helix)
  {
    HelixPoCAOptions helix_options;
    helix_options.magnetic_field_z = options.magnetic_field_z;
    const auto poca = computeHelixPoCA(
        {point_a, input_a.state_momentum, input_a.charge},
        {point_b, input_b.state_momentum, input_b.charge},
        helix_options);
    poca_valid = poca.valid;
    result.used_line_fallback = poca.used_line_fallback;
    result.s = poca.s;
    result.t = poca.t;
    result.dca = poca.dca;
    result.point_a = poca.point_a;
    result.point_b = poca.point_b;
    result.midpoint = poca.midpoint;
  }
  else
  {
    LocalLinePoCAOptions line_options;
    line_options.min_sin_angle = options.min_sin_angle;
    const auto poca = computeLocalLinePoCA(
        point_a,
        input_a.state_momentum,
        point_b,
        input_b.state_momentum,
        line_options);
    poca_valid = poca.valid;
    result.s = poca.s;
    result.t = poca.t;
    result.dca = poca.dca;
    result.point_a = poca.point_a;
    result.point_b = poca.point_b;
    result.midpoint = poca.midpoint;
  }

  if (!poca_valid)
  {
    result.status = CPMPairStatus::InvalidPoCA;
    return result;
  }

  if (!(result.dca <= options.max_pair_dca))
  {
    result.status = CPMPairStatus::DcaRejected;
    return result;
  }

  result.delta_x = voxel_center.x - result.midpoint.x;
  result.delta_y = voxel_center.y - result.midpoint.y;
  result.delta_z = voxel_center.z - result.midpoint.z;

  const double voxel_r = std::hypot(voxel_center.x, voxel_center.y);
  const double midpoint_r = std::hypot(result.midpoint.x, result.midpoint.y);
  const double voxel_phi = std::atan2(voxel_center.y, voxel_center.x);
  const double midpoint_phi = std::atan2(result.midpoint.y, result.midpoint.x);
  result.delta_r = voxel_r - midpoint_r;
  result.delta_phi = cpmWrapDeltaPhi(voxel_phi - midpoint_phi);
  result.delta_rphi = voxel_r * result.delta_phi;
  result.status = CPMPairStatus::Accepted;
  return result;
}
