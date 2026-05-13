#include "CPMHelixPoCA.h"
#include "CPMLocalLinePoCA.h"
#include "CPMVoxelContainer.h"

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
    const auto result = computeLocalLinePoCA(
        {-1.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, -1.0, 0.0},
        {0.0, 1.0, 0.0});

    assert(result.valid);
    assert(near(result.midpoint.x, 0.0));
    assert(near(result.midpoint.y, 0.0));
    assert(near(result.midpoint.z, 0.0));
    assert(near(result.dca, 0.0));
  }

  {
    const auto result = computeLocalLinePoCA(
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {1.0, 0.0, 0.0});

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

    VoxelContainer container;
    container.add(first);
    container.add(second);

    assert(container.voxel_count() == 1);
    assert(container.record_count() == 2);

    const auto records = container.find({1, 2, 3});
    assert(records != nullptr);
    assert(records->size() == 2);

    const auto result = computeVoxelCenterPoCA(records->at(0), records->at(1));
    assert(result.valid);
    assert(near(result.midpoint.x, 0.0));
    assert(near(result.midpoint.y, 0.0));
    assert(near(result.midpoint.z, 0.0));
  }

  {
    HelixPoCAOptions options;
    options.magnetic_field_z = 0.0;
    const auto result = computeHelixPoCA(
        {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, 1},
        {{0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, -1},
        options);

    assert(result.valid);
    assert(result.converged);
    assert(near(result.midpoint.x, 0.0));
    assert(near(result.midpoint.y, 0.0));
    assert(near(result.midpoint.z, 0.0));
    assert(near(result.dca, 0.0));
  }

  {
    const auto result = computeHelixPoCA(
        {{0.0, 0.0, 0.0}, {1.0, 0.2, 0.1}, 1},
        {{0.0, 0.0, 0.0}, {-0.2, 1.0, 0.1}, -1});

    assert(result.valid);
    assert(result.converged);
    assert(near(result.midpoint.x, 0.0, 1.0e-8));
    assert(near(result.midpoint.y, 0.0, 1.0e-8));
    assert(near(result.midpoint.z, 0.0, 1.0e-8));
    assert(near(result.dca, 0.0, 1.0e-8));
  }

  {
    HelixPoCAOptions options;
    options.magnetic_field_z = 0.0;
    const HelixState state{{1.0, 2.0, 3.0}, {0.0, 3.0, 4.0}, 1};
    const auto eval = evaluateHelix(state, 5.0, options);

    assert(eval.valid);
    assert(near(eval.position.x, 1.0));
    assert(near(eval.position.y, 5.0));
    assert(near(eval.position.z, 7.0));
    assert(near(eval.tangent.x, 0.0));
    assert(near(eval.tangent.y, 0.6));
    assert(near(eval.tangent.z, 0.8));
  }

  return 0;
}
