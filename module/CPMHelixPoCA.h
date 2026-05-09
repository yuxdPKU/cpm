#ifndef CPM_CPMHELIXPOCA_H
#define CPM_CPMHELIXPOCA_H

#include "CPMTypes.h"

#include <limits>

namespace cpm
{
  struct HelixState
  {
    Vector3 position;
    Vector3 momentum;
    int charge = 0;
  };

  struct HelixPoCAOptions
  {
    double magnetic_field_z = 1.4;
    double min_pt = 1.0e-6;
    double min_momentum = 1.0e-9;
    double min_hessian_determinant = 1.0e-18;
    double gradient_tolerance = 1.0e-9;
    double step_tolerance = 1.0e-9;
    double max_abs_path = 100.0;
    double max_step = 10.0;
    unsigned int max_iterations = 50;
    bool allow_line_fallback = true;
  };

  struct HelixEvaluation
  {
    Vector3 position;
    Vector3 tangent;
    Vector3 curvature;
    bool valid = false;
  };

  struct HelixPoCAResult
  {
    bool valid = false;
    bool converged = false;
    bool used_line_fallback = false;
    unsigned int iterations = 0;
    double s = 0.0;
    double t = 0.0;
    double dca = std::numeric_limits<double>::quiet_NaN();
    double gradient_norm = std::numeric_limits<double>::quiet_NaN();
    Vector3 point_a;
    Vector3 point_b;
    Vector3 midpoint;
  };

  HelixEvaluation evaluateHelix(
      const HelixState& state,
      double path,
      const HelixPoCAOptions& options = {});

  HelixPoCAResult computeHelixPoCA(
      const HelixState& state_a,
      const HelixState& state_b,
      const HelixPoCAOptions& options = {});
}

#endif
