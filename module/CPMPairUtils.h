#ifndef CPM_CPMPAIRUTILS_H
#define CPM_CPMPAIRUTILS_H

#include "CPMRecord.h"
#include "CPMTypes.h"

#include <limits>

enum class CPMPairSolver
{
  Line,
  Helix
};

enum class CPMPairStatus
{
  Accepted,
  PtRejected,
  InvalidWeight,
  InvalidPoCA,
  DcaRejected
};

struct CPMPairInput
{
  int charge = 0;
  double pt = std::numeric_limits<double>::quiet_NaN();
  Vector3 state_position;
  Vector3 state_momentum;
  Vector3 cluster_minus_voxel_center;
};

struct CPMPairOptions
{
  CPMPairSolver solver = CPMPairSolver::Helix;
  double min_pt = 0.5;
  double max_pair_dca = 2.0;
  double magnetic_field_z = 1.4;
  double min_sin_angle = 1.0e-4;
};

struct CPMPairResult
{
  CPMPairStatus status = CPMPairStatus::InvalidPoCA;
  bool used_line_fallback = false;
  double pair_weight = std::numeric_limits<double>::quiet_NaN();
  double dca = std::numeric_limits<double>::quiet_NaN();
  double s = std::numeric_limits<double>::quiet_NaN();
  double t = std::numeric_limits<double>::quiet_NaN();
  Vector3 point_a;
  Vector3 point_b;
  Vector3 midpoint;
  double delta_x = std::numeric_limits<double>::quiet_NaN();
  double delta_y = std::numeric_limits<double>::quiet_NaN();
  double delta_z = std::numeric_limits<double>::quiet_NaN();
  double delta_r = std::numeric_limits<double>::quiet_NaN();
  double delta_phi = std::numeric_limits<double>::quiet_NaN();
  double delta_rphi = std::numeric_limits<double>::quiet_NaN();

  [[nodiscard]] bool accepted() const { return status == CPMPairStatus::Accepted; }
};

bool cpmPairHasGoodPt(int charge, double pt, double min_pt);
double cpmPairWeightFromPt(double pt_a, double pt_b);
double cpmWrapDeltaPhi(double value);

CPMPairInput makeCPMPairInput(const TrackStateRecord& record);

CPMPairResult computeCPMPair(
    const CPMPairInput& input_a,
    const CPMPairInput& input_b,
    const Vector3& voxel_center,
    const CPMPairOptions& options = {});

#endif
