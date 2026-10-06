#include "CPMReconstructionHelper.h"
#include "CPMVoxelContainerv1.h"
#include "CPMAverageCorrectionReconstruction.h"
#include <TFile.h>
#include <TTree.h>
#include <cstdio>

#include <cassert>
#include <cmath>

namespace
{
  bool near(double lhs, double rhs, double tolerance = 1.0e-9)
  {
    return std::fabs(lhs - rhs) < tolerance;
  }
}

int main()
{
  {
    const auto result = CPMReconstructionHelper::compute_local_line_poca(
        {-1.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, -1.0, 0.0},
        {0.0, 1.0, 0.0},
        CPMReconstructionHelper::LocalLinePoCAOptions{});

    assert(result.valid);
    assert(near(result.midpoint.X(), 0.0));
    assert(near(result.midpoint.Y(), 0.0));
    assert(near(result.midpoint.Z(), 0.0));
    assert(near(result.dca, 0.0));
  }

  {
    const auto result = CPMReconstructionHelper::compute_local_line_poca(
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {1.0, 0.0, 0.0},
        CPMReconstructionHelper::LocalLinePoCAOptions{});

    assert(!result.valid);
  }

  {
    TrackStateRecord first;
    first.voxel = {1, 2, 3};
    first.state.position = {-1.0, 0.2, 0.0};
    first.state.momentum = {1.0, 0.0, 0.0};
    first.cluster.cluster_minus_voxel_center = {0.0, 0.2, 0.0};

    TrackStateRecord second;
    second.voxel = {1, 2, 3};
    second.state.position = {0.1, -1.0, 0.0};
    second.state.momentum = {0.0, 1.0, 0.0};
    second.cluster.cluster_minus_voxel_center = {0.1, 0.0, 0.0};

    CPMVoxelContainerv1 container;
    container.set_grid(36, 16, 80, 20.0, 78.0, -105.5, 105.5);
    container.add(first);
    container.add(second);

    assert(container.voxel_count() == 1);
    assert(container.record_count() == 2);

    const auto records = container.find({1, 2, 3});
    assert(records != nullptr);
    assert(records->size() == 2);
    assert(container.record_count({1, 2, 3}) == 2);
    assert(container.find_by_index(1, 2, 3) == records);

    const auto result = CPMReconstructionHelper::compute_voxel_center_poca(
        records->at(0),
        records->at(1),
        CPMReconstructionHelper::LocalLinePoCAOptions{});
    assert(result.valid);
    assert(near(result.midpoint.X(), 0.0));
    assert(near(result.midpoint.Y(), 0.0));
    assert(near(result.midpoint.Z(), 0.0));
  }

  {
    TrackStateRecord first;
    first.voxel = {0, 0, 0};
    first.event_ref.run = 1;
    first.event_ref.event_sequence = 5;
    first.track_ref.track_id = 2;
    first.cluster_ref.cluskey = 20;

    TrackStateRecord second;
    second.voxel = {0, 0, 0};
    second.event_ref.run = 1;
    second.event_ref.event_sequence = 3;
    second.track_ref.track_id = 1;
    second.cluster_ref.cluskey = 10;

    CPMVoxelContainerv1 lhs;
    lhs.set_grid(4, 4, 4, 0.0, 4.0, -2.0, 2.0);
    lhs.add(first);

    CPMVoxelContainerv1 rhs;
    rhs.add(second);
    lhs.add(rhs);

    const auto records = lhs.find_by_position(0.1, 0.5, -1.5);
    assert(records != nullptr);
    assert(records->size() == 2);
    assert(records->at(0).event_ref.event_sequence == 3);
    assert(records->at(1).event_ref.event_sequence == 5);
  }

  {
    CPMReconstructionHelper::HelixPoCAOptions options;
    options.magnetic_field_z = 0.0;
    const auto result = CPMReconstructionHelper::compute_helix_poca(
        {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, 1},
        {{0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, -1},
        options);

    assert(result.valid);
    assert(result.converged);
    assert(near(result.midpoint.X(), 0.0));
    assert(near(result.midpoint.Y(), 0.0));
    assert(near(result.midpoint.Z(), 0.0));
    assert(near(result.dca, 0.0));
  }

  {
    const auto result = CPMReconstructionHelper::compute_helix_poca(
        {{0.0, 0.0, 0.0}, {1.0, 0.2, 0.1}, 1},
        {{0.0, 0.0, 0.0}, {-0.2, 1.0, 0.1}, -1},
        CPMReconstructionHelper::HelixPoCAOptions{});

    assert(result.valid);
    assert(result.converged);
    assert(near(result.midpoint.X(), 0.0, 1.0e-8));
    assert(near(result.midpoint.Y(), 0.0, 1.0e-8));
    assert(near(result.midpoint.Z(), 0.0, 1.0e-8));
    assert(near(result.dca, 0.0, 1.0e-8));
  }

  {
    CPMReconstructionHelper::HelixPoCAOptions options;
    options.magnetic_field_z = 0.0;
    const CPMReconstructionHelper::HelixState state{{1.0, 2.0, 3.0}, {0.0, 3.0, 4.0}, 1};
    const auto eval = CPMReconstructionHelper::evaluate_helix(state, 5.0, options);

    assert(eval.valid);
    assert(near(eval.position.X(), 1.0));
    assert(near(eval.position.Y(), 5.0));
    assert(near(eval.position.Z(), 7.0));
    assert(near(eval.tangent.X(), 0.0));
    assert(near(eval.tangent.Y(), 0.6));
    assert(near(eval.tangent.Z(), 0.8));
  }

  {
    CPMReconstructionHelper::PairOptions options;
    options.solver = CPMReconstructionHelper::PairSolver::Line;
    options.min_pt = 0.5;
    const CPMReconstructionHelper::PairInput first{
        1,
        2.0,
        {-1.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 0.0, 0.0}};
    const CPMReconstructionHelper::PairInput second{
        -1,
        4.0,
        {0.0, -1.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0}};
    const auto result = CPMReconstructionHelper::compute_pair(first, second, {0.0, 0.0, 1.0}, options);

    assert(result.accepted());
    assert(near(result.pair_weight, 0.125));
    assert(near(result.delta_z, 1.0));

    CPMReconstructionHelper::CorrectionAccumulator accumulator;
    accumulator.add(
        result.delta_r,
        result.delta_rphi,
        result.delta_phi,
        result.delta_z,
        result.dca,
        result.pair_weight,
        0.0,
        0.0,
        1.0);
    assert(accumulator.entries == 1);
    assert(near(CPMReconstructionHelper::correction_weighted_mean(accumulator.sum_weighted_delta_z, accumulator.sum_pair_weight), 1.0));
  }

  // Identical translated origins with nonparallel tangents: zero residual stays zero.
  {
    const auto result = CPMReconstructionHelper::compute_pair(
        {1, 1.0, {60, 0, 0}, {1, .2, .1}, {}},
        {-1, 1.0, {60, 0, 0}, {1, -.2, -.1}, {}},
        {60, 0, 0}, CPMReconstructionHelper::PairOptions{});
    assert(result.accepted());
    assert(near(result.delta_r, 0.0));
  }
  // Curved tracks sampled away from a known local crossing recover that crossing.
  {
    const TVector3 crossing{60.1, .05, .02};
    CPMReconstructionHelper::HelixPoCAOptions options;
    const auto a = CPMReconstructionHelper::evaluate_helix(
        {crossing, TVector3{1, .3, .1}.Unit(), 1}, -.5, options);
    const auto b = CPMReconstructionHelper::evaluate_helix(
        {crossing, TVector3{.2, 1, -.1}.Unit(), -1}, .7, options);
    assert(a.valid && b.valid);
    const auto result = CPMReconstructionHelper::compute_pair(
        {1, 1.0, a.position, a.tangent, {}},
        {-1, 1.0, b.position, b.tangent, {}},
        {60, 0, 0}, CPMReconstructionHelper::PairOptions{});
    assert(result.accepted());
    assert(result.dca < 1e-3);
    assert((result.midpoint - crossing).Mag() < 1e-3);
  }
  // Parallel zero-DCA tracks do not determine an isolated crossing.
  {
    const auto result = CPMReconstructionHelper::compute_helix_poca(
        {{60, 0, 0}, {1, 0, 0}, 1},
        {{60, 0, 0}, {1, 0, 0}, -1}, {});
    assert(!result.valid && !result.converged);
  }
  // A zero-DCA intersection 60 cm away must fail in both solver modes.
  for (const auto solver : {CPMReconstructionHelper::PairSolver::Line,
                            CPMReconstructionHelper::PairSolver::Helix})
  {
    CPMReconstructionHelper::PairOptions options;
    options.solver = solver;
    options.magnetic_field_z = 0.0;
    options.allow_line_fallback = true;
    const auto result = CPMReconstructionHelper::compute_pair(
        {1, 1.0, {60, 0, 0}, {1, 0, 0}, {}},
        {-1, 1.0, {60, .6, 0}, {1, .01, 0}, {}}, {60, 0, 0}, options);
    assert(!result.accepted());
  }
  // Local path limits alone do not catch bad input origins outside the voxel.
  {
    CPMReconstructionHelper::PairOptions options;
    options.magnetic_field_z = 0.0;
    const auto result = CPMReconstructionHelper::compute_pair(
        {1, 1.0, {0, 0, 0}, {1, 0, 0}, {}},
        {-1, 1.0, {0, 0, 0}, {0, 1, 0}, {}}, {60, 0, 0}, options);
    assert(!result.accepted());
  }
  // Exhaustion, stagnation, and nonfinite inputs cannot become valid results.
  {
    CPMReconstructionHelper::HelixPoCAOptions options;
    options.magnetic_field_z = 0.0;
    options.max_iterations = 1;
    auto result = CPMReconstructionHelper::compute_helix_poca(
        {{-1, 0, 0}, {1, 0, 0}, 1}, {{0, -1, 0}, {0, 1, 0}, -1}, options);
    assert(!result.valid && !result.converged);
    options.max_iterations = 50;
    options.step_tolerance = 10;
    result = CPMReconstructionHelper::compute_helix_poca(
        {{-1, 0, 0}, {1, 0, 0}, 1}, {{0, -1, 0}, {0, 1, 0}, -1}, options);
    assert(!result.valid && !result.converged);
    options.step_tolerance = 1e-9;
    options.max_abs_path = std::numeric_limits<double>::infinity();
    result = CPMReconstructionHelper::compute_helix_poca(
        {{-1, 0, 0}, {1, 0, 0}, 1}, {{0, -1, 0}, {0, 1, 0}, -1}, options);
    assert(!result.valid);
  }
  // Raising the angular threshold is honored by helix and bounded fallback.
  {
    CPMReconstructionHelper::HelixPoCAOptions options;
    options.magnetic_field_z = 0.0;
    options.min_sin_angle = .1;
    options.allow_line_fallback = true;
    const auto result = CPMReconstructionHelper::compute_helix_poca(
        {{0, 0, 0}, {1, 0, 0}, 1}, {{0, 0, 0}, {1, .01, 0}, -1}, options);
    assert(!result.valid);
  }

  // Real ROOT summary round-trip must preserve cuts and reject mixed partials.
  {
    CPMVoxelContainerv1 empty;
    empty.set_grid(4, 4, 4, 20, 78, -100, 100);
    CPMAverageCorrectionReconstruction source;
    source.set_use_pair_weights(false);
    source.set_max_pair_dca(.2);
    source.set_max_abs_path(3);
    source.set_max_midpoint_distance(4);
    assert(source.add(empty));
    assert(source.finalize_average_corrections());
    assert(source.save_average_corrections("test_cpm_partial.root"));
    CPMAverageCorrectionReconstruction merged;
    merged.set_use_pair_weights(false);
    assert(merged.add_accumulators_from_file("test_cpm_partial.root"));
    assert(merged.add_accumulators_from_file("test_cpm_partial.root"));
    assert(merged.finalize_average_corrections());
    assert(merged.save_average_corrections("test_cpm_merged.root"));
    {
      TFile file("test_cpm_merged.root", "READ");
      auto* summary = dynamic_cast<TTree*>(file.Get("cpm_b3_summary"));
      assert(summary);
      double dca = 0, path = 0, distance = 0;
      unsigned int version = 0;
      summary->SetBranchAddress("max_pair_dca", &dca);
      summary->SetBranchAddress("max_abs_path", &path);
      summary->SetBranchAddress("max_midpoint_distance", &distance);
      summary->SetBranchAddress("solver_protection_version", &version);
      summary->GetEntry(0);
      assert(near(dca, .2) && near(path, 3) && near(distance, 4) && version == 1);
    }
    source.set_max_abs_path(2);
    assert(source.save_average_corrections("test_cpm_mismatch.root"));
    assert(!merged.add_accumulators_from_file("test_cpm_mismatch.root"));
    CPMAverageCorrectionReconstruction weighted;
    assert(!weighted.add_accumulators_from_file("test_cpm_partial.root"));
    std::remove("test_cpm_partial.root");
    std::remove("test_cpm_merged.root");
    std::remove("test_cpm_mismatch.root");
  }
  return 0;
}
