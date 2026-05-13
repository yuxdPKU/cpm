#include "CPMHelixPoCA.h"

#include "CPMLocalLinePoCA.h"

#include <cmath>

namespace
{
  constexpr double kCurvaturePerGeVPerTesla = 0.003;

  bool finite_vector(const Vector3& value)
  {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
  }

  double transverse_momentum(const Vector3& momentum)
  {
    return std::hypot(momentum.x, momentum.y);
  }

  double signed_omega(
      const HelixState& state,
      const HelixPoCAOptions& options,
      const double momentum_norm)
  {
    return -static_cast<double>(state.charge) *
           kCurvaturePerGeVPerTesla *
           options.magnetic_field_z / momentum_norm;
  }

  double clamp_step(const double step, const double max_step)
  {
    if (!(max_step > 0.0) || std::fabs(step) <= max_step)
    {
      return step;
    }
    return std::copysign(max_step, step);
  }

  bool within_path_window(
      const double s,
      const double t,
      const HelixPoCAOptions& options)
  {
    if (!(options.max_abs_path > 0.0))
    {
      return true;
    }
    return std::fabs(s) <= options.max_abs_path &&
           std::fabs(t) <= options.max_abs_path;
  }

  HelixPoCAResult line_fallback(
      const HelixState& state_a,
      const HelixState& state_b)
  {
    HelixPoCAResult result;
    const auto line_result = computeLocalLinePoCA(
        state_a.position,
        state_a.momentum,
        state_b.position,
        state_b.momentum);

    if (!line_result.valid)
    {
      return result;
    }

    result.valid = true;
    result.converged = true;
    result.used_line_fallback = true;
    result.s = line_result.s;
    result.t = line_result.t;
    result.dca = line_result.dca;
    result.gradient_norm = 0.0;
    result.point_a = line_result.point_a;
    result.point_b = line_result.point_b;
    result.midpoint = line_result.midpoint;
    return result;
  }
}

HelixEvaluation evaluateHelix(
    const HelixState& state,
    const double path,
    const HelixPoCAOptions& options)
{
  HelixEvaluation result;
  if (!finite_vector(state.position) ||
      !finite_vector(state.momentum) ||
      !std::isfinite(path) ||
      !std::isfinite(options.magnetic_field_z))
  {
    return result;
  }

  const double momentum_norm = state.momentum.norm();
  const double pt = transverse_momentum(state.momentum);
  if (!(momentum_norm > options.min_momentum) || !(pt > options.min_pt))
  {
    return result;
  }

  const Vector3 tangent0 = state.momentum / momentum_norm;
  const double transverse_speed = pt / momentum_norm;
  const double theta0 = std::atan2(tangent0.y, tangent0.x);
  const double omega = signed_omega(state, options, momentum_norm);

  if (state.charge == 0 || std::fabs(omega) < 1.0e-18)
  {
    result.position = state.position + tangent0 * path;
    result.tangent = tangent0;
    result.curvature = {0.0, 0.0, 0.0};
    result.valid = true;
    return result;
  }

  const double theta = theta0 + omega * path;
  const double sin_theta = std::sin(theta);
  const double cos_theta = std::cos(theta);
  const double sin_theta0 = std::sin(theta0);
  const double cos_theta0 = std::cos(theta0);

  result.position = {
      state.position.x + transverse_speed / omega * (sin_theta - sin_theta0),
      state.position.y - transverse_speed / omega * (cos_theta - cos_theta0),
      state.position.z + tangent0.z * path};
  result.tangent = {
      transverse_speed * cos_theta,
      transverse_speed * sin_theta,
      tangent0.z};
  result.curvature = {
      -transverse_speed * omega * sin_theta,
      transverse_speed * omega * cos_theta,
      0.0};
  result.valid = finite_vector(result.position) &&
                 finite_vector(result.tangent) &&
                 finite_vector(result.curvature);
  return result;
}

HelixPoCAResult computeHelixPoCA(
    const HelixState& state_a,
    const HelixState& state_b,
    const HelixPoCAOptions& options)
{
  HelixPoCAResult result;

  const auto initial = computeLocalLinePoCA(
      state_a.position,
      state_a.momentum,
      state_b.position,
      state_b.momentum);
  double s = initial.valid ? initial.s : 0.0;
  double t = initial.valid ? initial.t : 0.0;

  if (!within_path_window(s, t, options))
  {
    s = 0.0;
    t = 0.0;
  }

  for (unsigned int iteration = 0; iteration < options.max_iterations; ++iteration)
  {
    const auto eval_a = evaluateHelix(state_a, s, options);
    const auto eval_b = evaluateHelix(state_b, t, options);
    if (!eval_a.valid || !eval_b.valid)
    {
      if (options.allow_line_fallback)
      {
        return line_fallback(state_a, state_b);
      }
      return result;
    }

    const Vector3 delta = eval_a.position - eval_b.position;
    const double gradient_s = dot(delta, eval_a.tangent);
    const double gradient_t = -dot(delta, eval_b.tangent);
    result.gradient_norm = std::hypot(gradient_s, gradient_t);

    result.iterations = iteration + 1;
    result.s = s;
    result.t = t;
    result.point_a = eval_a.position;
    result.point_b = eval_b.position;
    result.midpoint = (result.point_a + result.point_b) * 0.5;
    result.dca = delta.norm();
    result.valid = std::isfinite(result.dca);

    if (result.gradient_norm <= options.gradient_tolerance)
    {
      result.converged = result.valid;
      return result;
    }

    const double hss = dot(eval_a.tangent, eval_a.tangent) + dot(delta, eval_a.curvature);
    const double htt = dot(eval_b.tangent, eval_b.tangent) - dot(delta, eval_b.curvature);
    const double hst = -dot(eval_a.tangent, eval_b.tangent);
    const double determinant = hss * htt - hst * hst;

    if (std::fabs(determinant) <= options.min_hessian_determinant)
    {
      if (options.allow_line_fallback)
      {
        return line_fallback(state_a, state_b);
      }
      return result;
    }

    double delta_s = (-gradient_s * htt + hst * gradient_t) / determinant;
    double delta_t = (hst * gradient_s - hss * gradient_t) / determinant;
    delta_s = clamp_step(delta_s, options.max_step);
    delta_t = clamp_step(delta_t, options.max_step);

    if (std::hypot(delta_s, delta_t) <= options.step_tolerance)
    {
      result.converged = result.valid;
      return result;
    }

    s += delta_s;
    t += delta_t;
    if (!within_path_window(s, t, options))
    {
      return result;
    }
  }

  return result;
}
