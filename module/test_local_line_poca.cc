#include "CPMReconstructionHelper.h"
#include "CPMVoxelContainerv1.h"

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

  return 0;
}
