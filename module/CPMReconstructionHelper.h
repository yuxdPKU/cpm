#ifndef CPM_CPMRECONSTRUCTIONHELPER_H
#define CPM_CPMRECONSTRUCTIONHELPER_H

#include "CPMVoxelContainer.h"

#include <TH3.h>
#include <TString.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

class TFile;

class CPMReconstructionHelper
{
 public:
  static constexpr double Pi = 3.14159265358979323846;

  struct EventTrackKey
  {
    std::string cluster_source;
    std::string track_source;
    int run = -1;
    int segment = -1;
    int sync_event = -1;
    int event_sequence = -1;
    std::uint64_t stream_event_ordinal = 0;
    unsigned int track_id = 0;

    friend bool operator==(const EventTrackKey& lhs, const EventTrackKey& rhs)
    {
      return std::tie(lhs.cluster_source, lhs.track_source, lhs.run, lhs.segment,
                      lhs.sync_event, lhs.event_sequence,
                      lhs.stream_event_ordinal, lhs.track_id) ==
             std::tie(rhs.cluster_source, rhs.track_source, rhs.run, rhs.segment,
                      rhs.sync_event, rhs.event_sequence,
                      rhs.stream_event_ordinal, rhs.track_id);
    }

    friend bool operator<(const EventTrackKey& lhs, const EventTrackKey& rhs)
    {
      return std::tie(lhs.cluster_source, lhs.track_source, lhs.run, lhs.segment,
                      lhs.sync_event, lhs.event_sequence,
                      lhs.stream_event_ordinal, lhs.track_id) <
             std::tie(rhs.cluster_source, rhs.track_source, rhs.run, rhs.segment,
                      rhs.sync_event, rhs.event_sequence,
                      rhs.stream_event_ordinal, rhs.track_id);
    }
  };

  using HistogramPair = std::tuple<std::unique_ptr<TH3>, std::unique_ptr<TH3>>;

  enum class PairSolver
  {
    Line,
    Helix
  };

  enum class PairStatus
  {
    Accepted,
    PtRejected,
    InvalidWeight,
    InvalidPoCA,
    DcaRejected
  };

  struct LocalLinePoCAOptions
  {
    double min_sin_angle = 1.0e-4;
  };

  struct LocalLinePoCAResult
  {
    bool valid = false;
    double s = 0.0;
    double t = 0.0;
    double dca = std::numeric_limits<double>::quiet_NaN();
    TVector3 point_a;
    TVector3 point_b;
    TVector3 midpoint;
  };

  struct HelixState
  {
    TVector3 position;
    TVector3 momentum;
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
    TVector3 position;
    TVector3 tangent;
    TVector3 curvature;
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
    TVector3 point_a;
    TVector3 point_b;
    TVector3 midpoint;
  };

  struct PairInput
  {
    int charge = 0;
    double pt = std::numeric_limits<double>::quiet_NaN();
    TVector3 state_position;
    TVector3 state_momentum;
    TVector3 cluster_minus_voxel_center;
  };

  struct PairOptions
  {
    PairSolver solver = PairSolver::Helix;
    double min_pt = 0.5;
    double max_pair_dca = 2.0;
    double magnetic_field_z = 1.4;
    double min_sin_angle = 1.0e-4;
  };

  struct PairResult
  {
    PairStatus status = PairStatus::InvalidPoCA;
    bool used_line_fallback = false;
    double pair_weight = std::numeric_limits<double>::quiet_NaN();
    double dca = std::numeric_limits<double>::quiet_NaN();
    double s = std::numeric_limits<double>::quiet_NaN();
    double t = std::numeric_limits<double>::quiet_NaN();
    TVector3 point_a;
    TVector3 point_b;
    TVector3 midpoint;
    double delta_x = std::numeric_limits<double>::quiet_NaN();
    double delta_y = std::numeric_limits<double>::quiet_NaN();
    double delta_z = std::numeric_limits<double>::quiet_NaN();
    double delta_r = std::numeric_limits<double>::quiet_NaN();
    double delta_phi = std::numeric_limits<double>::quiet_NaN();
    double delta_rphi = std::numeric_limits<double>::quiet_NaN();

    [[nodiscard]] bool accepted() const { return status == PairStatus::Accepted; }
  };

  struct CorrectionAccumulator
  {
    unsigned long long entries = 0;
    double sum_pair_weight = 0.0;
    double sum_pair_weight2 = 0.0;
    double sum_delta_r = 0.0;
    double sum_delta_r2 = 0.0;
    double sum_delta_rphi = 0.0;
    double sum_delta_rphi2 = 0.0;
    double sum_delta_phi = 0.0;
    double sum_delta_phi2 = 0.0;
    double sum_delta_z = 0.0;
    double sum_delta_z2 = 0.0;
    double sum_weighted_delta_r = 0.0;
    double sum_weighted_delta_r2 = 0.0;
    double sum_weighted_delta_rphi = 0.0;
    double sum_weighted_delta_rphi2 = 0.0;
    double sum_weighted_delta_phi = 0.0;
    double sum_weighted_delta_phi2 = 0.0;
    double sum_weighted_delta_z = 0.0;
    double sum_weighted_delta_z2 = 0.0;
    double sum_dca = 0.0;
    double sum_dca2 = 0.0;
    double sum_voxel_x = 0.0;
    double sum_voxel_y = 0.0;
    double sum_voxel_z = 0.0;

    void add(
        double delta_r,
        double delta_rphi,
        double delta_phi,
        double delta_z,
        double dca,
        double weight,
        double voxel_x,
        double voxel_y,
        double voxel_z)
    {
      ++entries;
      sum_pair_weight += weight;
      sum_pair_weight2 += weight * weight;
      sum_delta_r += delta_r;
      sum_delta_r2 += delta_r * delta_r;
      sum_delta_rphi += delta_rphi;
      sum_delta_rphi2 += delta_rphi * delta_rphi;
      sum_delta_phi += delta_phi;
      sum_delta_phi2 += delta_phi * delta_phi;
      sum_delta_z += delta_z;
      sum_delta_z2 += delta_z * delta_z;
      sum_weighted_delta_r += weight * delta_r;
      sum_weighted_delta_r2 += weight * delta_r * delta_r;
      sum_weighted_delta_rphi += weight * delta_rphi;
      sum_weighted_delta_rphi2 += weight * delta_rphi * delta_rphi;
      sum_weighted_delta_phi += weight * delta_phi;
      sum_weighted_delta_phi2 += weight * delta_phi * delta_phi;
      sum_weighted_delta_z += weight * delta_z;
      sum_weighted_delta_z2 += weight * delta_z * delta_z;
      sum_dca += dca;
      sum_dca2 += dca * dca;
      sum_voxel_x += voxel_x;
      sum_voxel_y += voxel_y;
      sum_voxel_z += voxel_z;
    }

    void add_sums(
        unsigned long long added_entries,
        double added_sum_pair_weight,
        double added_sum_pair_weight2,
        double added_sum_delta_r,
        double added_sum_delta_r2,
        double added_sum_delta_rphi,
        double added_sum_delta_rphi2,
        double added_sum_delta_phi,
        double added_sum_delta_phi2,
        double added_sum_delta_z,
        double added_sum_delta_z2,
        double added_sum_dca,
        double added_sum_dca2,
        double voxel_x,
        double voxel_y,
        double voxel_z)
    {
      entries += added_entries;
      sum_pair_weight += added_sum_pair_weight;
      sum_pair_weight2 += added_sum_pair_weight2;
      sum_weighted_delta_r += added_sum_delta_r;
      sum_weighted_delta_r2 += added_sum_delta_r2;
      sum_weighted_delta_rphi += added_sum_delta_rphi;
      sum_weighted_delta_rphi2 += added_sum_delta_rphi2;
      sum_weighted_delta_phi += added_sum_delta_phi;
      sum_weighted_delta_phi2 += added_sum_delta_phi2;
      sum_weighted_delta_z += added_sum_delta_z;
      sum_weighted_delta_z2 += added_sum_delta_z2;
      sum_dca += added_sum_dca;
      sum_dca2 += added_sum_dca2;
      sum_voxel_x += voxel_x * static_cast<double>(added_entries);
      sum_voxel_y += voxel_y * static_cast<double>(added_entries);
      sum_voxel_z += voxel_z * static_cast<double>(added_entries);
    }
  };

  static double correction_mean(
      double sum,
      unsigned long long entries)
  {
    return entries > 0 ? sum / static_cast<double>(entries) : std::numeric_limits<double>::quiet_NaN();
  }

  static double correction_weighted_mean(
      double sum_weighted,
      double sum_weight)
  {
    return sum_weight > 0.0 ? sum_weighted / sum_weight : std::numeric_limits<double>::quiet_NaN();
  }

  static double correction_rms(
      double sum,
      double sum2,
      unsigned long long entries)
  {
    if (entries == 0)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double average = correction_mean(sum, entries);
    const double variance = sum2 / static_cast<double>(entries) - average * average;
    return std::sqrt(std::max(0.0, variance));
  }

  static double correction_weighted_rms(
      double sum_weighted,
      double sum_weighted2,
      double sum_weight)
  {
    if (sum_weight <= 0.0)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double average = correction_weighted_mean(sum_weighted, sum_weight);
    const double variance = sum_weighted2 / sum_weight - average * average;
    return std::sqrt(std::max(0.0, variance));
  }

  static double correction_effective_entries(
      double sum_weight,
      double sum_weight2)
  {
    return sum_weight2 > 0.0 ?
        sum_weight * sum_weight / sum_weight2 :
        std::numeric_limits<double>::quiet_NaN();
  }

  static std::vector<std::string> read_file_list(const std::string& input_list);

  static LocalLinePoCAResult compute_local_line_poca(
      const TVector3& point_a,
      const TVector3& direction_a,
      const TVector3& point_b,
      const TVector3& direction_b,
      const LocalLinePoCAOptions& options);
  static LocalLinePoCAResult compute_voxel_center_poca(
      const TrackStateRecord& record_a,
      const TrackStateRecord& record_b,
      const LocalLinePoCAOptions& options);

  static HelixEvaluation evaluate_helix(
      const HelixState& state,
      double path,
      const HelixPoCAOptions& options);
  static HelixPoCAResult compute_helix_poca(
      const HelixState& state_a,
      const HelixState& state_b,
      const HelixPoCAOptions& options);

  static bool pair_has_good_pt(int charge, double pt, double min_pt);
  static double pair_weight_from_pt(double pt_a, double pt_b);
  static double wrap_delta_phi(double value);
  static PairInput make_pair_input(const TrackStateRecord& record);
  static PairResult compute_pair(
      const PairInput& input_a,
      const PairInput& input_b,
      const TVector3& voxel_center,
      const PairOptions& options);

  static EventTrackKey event_track_key(const TrackStateRecord& record);
  static std::uint64_t record_selection_hash(const TrackStateRecord& record);
  static double offset_magnitude2(const TrackStateRecord& record);
  static bool is_closer_to_voxel_center(
      const TrackStateRecord& candidate,
      const TrackStateRecord& current_best);

  static bool same_grid(
      const CPMVoxelContainer::Grid& lhs,
      const CPMVoxelContainer::Grid& rhs);

  static PairSolver solver_from_string(
      const std::string& name,
      bool& ok);
  static const char* solver_name(PairSolver solver);

  static HistogramPair make_guarded_histograms(
      TH3* source,
      const TString& name);
  static void write_guarded_histograms(
      TFile& output,
      TH3* source,
      const TString& name);
};

#endif
