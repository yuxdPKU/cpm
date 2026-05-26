#include "CPMReconstructionHelper.h"

#include <tpccalib/TpcSpaceChargeReconstructionHelper.h>

#include <TFile.h>

#include <cmath>
#include <fstream>

std::vector<std::string> CPMReconstructionHelper::read_file_list(
    const std::string& input_list)
{
  std::ifstream input(input_list);
  std::vector<std::string> files;
  std::string line;
  while (std::getline(input, line))
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    files.push_back(line);
  }
  return files;
}

CPMReconstructionHelper::LocalLinePoCAResult
CPMReconstructionHelper::compute_local_line_poca(
    const TVector3& point_a,
    const TVector3& direction_a,
    const TVector3& point_b,
    const TVector3& direction_b,
    const LocalLinePoCAOptions& options)
{
  LocalLinePoCAResult result;

  const TVector3 u = direction_a.Mag() > 0.0 ? direction_a.Unit() : TVector3{};
  const TVector3 v = direction_b.Mag() > 0.0 ? direction_b.Unit() : TVector3{};
  const TVector3 w0 = point_a - point_b;

  const double a = u.Dot(u);
  const double b = u.Dot(v);
  const double c = v.Dot(v);
  const double d = u.Dot(w0);
  const double e = v.Dot(w0);
  const double denom = a * c - b * b;

  if (!(a > 0.0 && c > 0.0) ||
      denom <= options.min_sin_angle * options.min_sin_angle)
  {
    return result;
  }

  result.s = (b * e - c * d) / denom;
  result.t = (a * e - b * d) / denom;
  result.point_a = point_a + u * result.s;
  result.point_b = point_b + v * result.t;
  result.midpoint = (result.point_a + result.point_b) * 0.5;
  result.dca = (result.point_a - result.point_b).Mag();
  result.valid = std::isfinite(result.dca);

  return result;
}

CPMReconstructionHelper::LocalLinePoCAResult
CPMReconstructionHelper::compute_voxel_center_poca(
    const TrackStateRecord& record_a,
    const TrackStateRecord& record_b,
    const LocalLinePoCAOptions& options)
{
  const TVector3 point_a =
      record_a.state.position - record_a.cluster.cluster_minus_voxel_center;
  const TVector3 point_b =
      record_b.state.position - record_b.cluster.cluster_minus_voxel_center;

  return compute_local_line_poca(
      point_a,
      record_a.state.momentum,
      point_b,
      record_b.state.momentum,
      options);
}

namespace
{
  constexpr double kCurvaturePerGeVPerTesla = 0.003;

  bool finite_vector(const TVector3& value)
  {
    return std::isfinite(value.X()) &&
           std::isfinite(value.Y()) &&
           std::isfinite(value.Z());
  }

  double transverse_momentum(const TVector3& momentum)
  {
    return std::hypot(momentum.X(), momentum.Y());
  }

  double signed_omega(
      const CPMReconstructionHelper::HelixState& state,
      const CPMReconstructionHelper::HelixPoCAOptions& options,
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
      const CPMReconstructionHelper::HelixPoCAOptions& options)
  {
    if (!(options.max_abs_path > 0.0))
    {
      return true;
    }
    return std::fabs(s) <= options.max_abs_path &&
           std::fabs(t) <= options.max_abs_path;
  }

  CPMReconstructionHelper::HelixPoCAResult line_fallback(
      const CPMReconstructionHelper::HelixState& state_a,
      const CPMReconstructionHelper::HelixState& state_b)
  {
    CPMReconstructionHelper::HelixPoCAResult result;
    const auto line_result = CPMReconstructionHelper::compute_local_line_poca(
        state_a.position,
        state_a.momentum,
        state_b.position,
        state_b.momentum,
        CPMReconstructionHelper::LocalLinePoCAOptions{});

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
}  // namespace

CPMReconstructionHelper::HelixEvaluation
CPMReconstructionHelper::evaluate_helix(
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

  const double momentum_norm = state.momentum.Mag();
  const double pt = transverse_momentum(state.momentum);
  if (!(momentum_norm > options.min_momentum) || !(pt > options.min_pt))
  {
    return result;
  }

  const TVector3 tangent0 = state.momentum * (1.0 / momentum_norm);
  const double transverse_speed = pt / momentum_norm;
  const double theta0 = std::atan2(tangent0.Y(), tangent0.X());
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
      state.position.X() + transverse_speed / omega * (sin_theta - sin_theta0),
      state.position.Y() - transverse_speed / omega * (cos_theta - cos_theta0),
      state.position.Z() + tangent0.Z() * path};
  result.tangent = {
      transverse_speed * cos_theta,
      transverse_speed * sin_theta,
      tangent0.Z()};
  result.curvature = {
      -transverse_speed * omega * sin_theta,
      transverse_speed * omega * cos_theta,
      0.0};
  result.valid = finite_vector(result.position) &&
                 finite_vector(result.tangent) &&
                 finite_vector(result.curvature);
  return result;
}

CPMReconstructionHelper::HelixPoCAResult
CPMReconstructionHelper::compute_helix_poca(
    const HelixState& state_a,
    const HelixState& state_b,
    const HelixPoCAOptions& options)
{
  HelixPoCAResult result;

  const auto initial = compute_local_line_poca(
      state_a.position,
      state_a.momentum,
      state_b.position,
      state_b.momentum,
      LocalLinePoCAOptions{});
  double s = initial.valid ? initial.s : 0.0;
  double t = initial.valid ? initial.t : 0.0;

  if (!within_path_window(s, t, options))
  {
    s = 0.0;
    t = 0.0;
  }

  for (unsigned int iteration = 0; iteration < options.max_iterations; ++iteration)
  {
    const auto eval_a = evaluate_helix(state_a, s, options);
    const auto eval_b = evaluate_helix(state_b, t, options);
    if (!eval_a.valid || !eval_b.valid)
    {
      if (options.allow_line_fallback)
      {
        return line_fallback(state_a, state_b);
      }
      return result;
    }

    const TVector3 delta = eval_a.position - eval_b.position;
    const double gradient_s = delta.Dot(eval_a.tangent);
    const double gradient_t = -delta.Dot(eval_b.tangent);
    result.gradient_norm = std::hypot(gradient_s, gradient_t);

    result.iterations = iteration + 1;
    result.s = s;
    result.t = t;
    result.point_a = eval_a.position;
    result.point_b = eval_b.position;
    result.midpoint = (result.point_a + result.point_b) * 0.5;
    result.dca = delta.Mag();
    result.valid = std::isfinite(result.dca);

    if (result.gradient_norm <= options.gradient_tolerance)
    {
      result.converged = result.valid;
      return result;
    }

    const double hss =
        eval_a.tangent.Dot(eval_a.tangent) + delta.Dot(eval_a.curvature);
    const double htt =
        eval_b.tangent.Dot(eval_b.tangent) - delta.Dot(eval_b.curvature);
    const double hst = -eval_a.tangent.Dot(eval_b.tangent);
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

bool CPMReconstructionHelper::pair_has_good_pt(
    const int charge,
    const double pt,
    const double min_pt)
{
  return charge != 0 &&
         std::isfinite(pt) &&
         pt > 0.0 &&
         pt >= min_pt;
}

double CPMReconstructionHelper::pair_weight_from_pt(
    const double pt_a,
    const double pt_b)
{
  return 1.0 / (pt_a * pt_b);
}

double CPMReconstructionHelper::wrap_delta_phi(double value)
{
  while (value > Pi)
  {
    value -= 2.0 * Pi;
  }
  while (value <= -Pi)
  {
    value += 2.0 * Pi;
  }
  return value;
}

CPMReconstructionHelper::PairInput
CPMReconstructionHelper::make_pair_input(const TrackStateRecord& record)
{
  return {
      record.track.charge,
      record.track.pt,
      record.state.position,
      record.state.momentum,
      record.cluster.cluster_minus_voxel_center};
}

CPMReconstructionHelper::PairResult
CPMReconstructionHelper::compute_pair(
    const PairInput& input_a,
    const PairInput& input_b,
    const TVector3& voxel_center,
    const PairOptions& options)
{
  PairResult result;

  if (!pair_has_good_pt(input_a.charge, input_a.pt, options.min_pt) ||
      !pair_has_good_pt(input_b.charge, input_b.pt, options.min_pt))
  {
    result.status = PairStatus::PtRejected;
    return result;
  }

  result.pair_weight = pair_weight_from_pt(input_a.pt, input_b.pt);
  if (!(std::isfinite(result.pair_weight) && result.pair_weight > 0.0))
  {
    result.status = PairStatus::InvalidWeight;
    return result;
  }

  const TVector3 point_a =
      input_a.state_position - input_a.cluster_minus_voxel_center;
  const TVector3 point_b =
      input_b.state_position - input_b.cluster_minus_voxel_center;

  bool poca_valid = false;
  if (options.solver == PairSolver::Helix)
  {
    HelixPoCAOptions helix_options;
    helix_options.magnetic_field_z = options.magnetic_field_z;
    const auto poca = compute_helix_poca(
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
    const auto poca = compute_local_line_poca(
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
    result.status = PairStatus::InvalidPoCA;
    return result;
  }

  if (!(result.dca <= options.max_pair_dca))
  {
    result.status = PairStatus::DcaRejected;
    return result;
  }

  result.delta_x = voxel_center.X() - result.midpoint.X();
  result.delta_y = voxel_center.Y() - result.midpoint.Y();
  result.delta_z = voxel_center.Z() - result.midpoint.Z();

  const double voxel_r = std::hypot(voxel_center.X(), voxel_center.Y());
  const double midpoint_r = std::hypot(result.midpoint.X(), result.midpoint.Y());
  const double voxel_phi = std::atan2(voxel_center.Y(), voxel_center.X());
  const double midpoint_phi = std::atan2(result.midpoint.Y(), result.midpoint.X());
  result.delta_r = voxel_r - midpoint_r;
  result.delta_phi = wrap_delta_phi(voxel_phi - midpoint_phi);
  result.delta_rphi = voxel_r * result.delta_phi;
  result.status = PairStatus::Accepted;
  return result;
}

CPMReconstructionHelper::EventTrackKey
CPMReconstructionHelper::event_track_key(const TrackStateRecord& record)
{
  return {
      record.event_ref.cluster_source,
      record.event_ref.track_source,
      record.event_ref.run,
      record.event_ref.segment,
      record.event_ref.sync_event,
      record.event_ref.event_sequence,
      record.event_ref.stream_event_ordinal,
      record.track_ref.track_id};
}

namespace
{
  void hash_mix_byte(std::uint64_t& hash, const unsigned char value)
  {
    constexpr std::uint64_t fnv_prime = 1099511628211ULL;
    hash ^= value;
    hash *= fnv_prime;
  }

  void hash_mix_uint64(std::uint64_t& hash, std::uint64_t value)
  {
    for (int i = 0; i < 8; ++i)
    {
      hash_mix_byte(hash, static_cast<unsigned char>(value & 0xffU));
      value >>= 8U;
    }
  }

  void hash_mix_string(std::uint64_t& hash, const std::string& value)
  {
    for (const auto character : value)
    {
      hash_mix_byte(hash, static_cast<unsigned char>(character));
    }
    hash_mix_byte(hash, 0xffU);
  }

  bool same_grid_value(const double lhs, const double rhs)
  {
    return std::fabs(lhs - rhs) < 1.0e-9;
  }
}  // namespace

std::uint64_t CPMReconstructionHelper::record_selection_hash(
    const TrackStateRecord& record)
{
  std::uint64_t hash = 14695981039346656037ULL;
  hash_mix_string(hash, record.event_ref.cluster_source);
  hash_mix_string(hash, record.event_ref.track_source);
  hash_mix_uint64(hash, static_cast<std::uint64_t>(record.event_ref.run));
  hash_mix_uint64(hash, static_cast<std::uint64_t>(record.event_ref.segment));
  hash_mix_uint64(hash, static_cast<std::uint64_t>(record.event_ref.sync_event));
  hash_mix_uint64(hash, static_cast<std::uint64_t>(record.event_ref.event_sequence));
  hash_mix_uint64(hash, record.event_ref.stream_event_ordinal);
  hash_mix_uint64(hash, record.track_ref.track_id);
  hash_mix_uint64(hash, record.cluster_ref.cluskey);
  return hash;
}

double CPMReconstructionHelper::offset_magnitude2(const TrackStateRecord& record)
{
  const auto& offset = record.cluster.cluster_minus_voxel_center;
  return offset.X() * offset.X() +
         offset.Y() * offset.Y() +
         offset.Z() * offset.Z();
}

bool CPMReconstructionHelper::is_closer_to_voxel_center(
    const TrackStateRecord& candidate,
    const TrackStateRecord& current_best)
{
  const double candidate_distance2 = offset_magnitude2(candidate);
  const double current_distance2 = offset_magnitude2(current_best);
  if (std::isfinite(candidate_distance2) &&
      std::isfinite(current_distance2) &&
      candidate_distance2 != current_distance2)
  {
    return candidate_distance2 < current_distance2;
  }
  if (std::isfinite(candidate_distance2) != std::isfinite(current_distance2))
  {
    return std::isfinite(candidate_distance2);
  }
  return candidate.cluster_ref.cluskey < current_best.cluster_ref.cluskey;
}

bool CPMReconstructionHelper::same_grid(
    const CPMVoxelContainer::Grid& lhs,
    const CPMVoxelContainer::Grid& rhs)
{
  return lhs.valid() &&
         rhs.valid() &&
         lhs.phi_bins == rhs.phi_bins &&
         lhs.r_bins == rhs.r_bins &&
         lhs.z_bins == rhs.z_bins &&
         same_grid_value(lhs.phi_min, rhs.phi_min) &&
         same_grid_value(lhs.phi_max, rhs.phi_max) &&
         same_grid_value(lhs.r_min, rhs.r_min) &&
         same_grid_value(lhs.r_max, rhs.r_max) &&
         same_grid_value(lhs.z_min, rhs.z_min) &&
         same_grid_value(lhs.z_max, rhs.z_max);
}

CPMReconstructionHelper::PairSolver
CPMReconstructionHelper::solver_from_string(
    const std::string& name,
    bool& ok)
{
  if (name == "helix")
  {
    ok = true;
    return PairSolver::Helix;
  }
  if (name == "line")
  {
    ok = true;
    return PairSolver::Line;
  }
  ok = false;
  return PairSolver::Helix;
}

const char* CPMReconstructionHelper::solver_name(const PairSolver solver)
{
  return solver == PairSolver::Line ? "line" : "helix";
}

CPMReconstructionHelper::HistogramPair
CPMReconstructionHelper::make_guarded_histograms(
    TH3* source,
    const TString& name)
{
  const auto result = TpcSpaceChargeReconstructionHelper::split(source);
  std::unique_ptr<TH3> hneg(std::get<0>(result));
  std::unique_ptr<TH3> hpos(std::get<1>(result));

  return std::make_tuple(
      std::unique_ptr<TH3>(TpcSpaceChargeReconstructionHelper::add_guarding_bins(
          hneg.get(),
          name + "_negz")),
      std::unique_ptr<TH3>(TpcSpaceChargeReconstructionHelper::add_guarding_bins(
          hpos.get(),
          name + "_posz")));
}

void CPMReconstructionHelper::write_guarded_histograms(
    TFile& output,
    TH3* source,
    const TString& name)
{
  auto finished = make_guarded_histograms(source, name);
  output.cd();
  std::get<0>(finished)->Write();
  std::get<1>(finished)->Write();
}
