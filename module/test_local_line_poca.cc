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
    const auto result = cpm::computeLocalLinePoCA(
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
    const auto result = cpm::computeLocalLinePoCA(
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {1.0, 0.0, 0.0});

    assert(!result.valid);
  }

  {
    cpm::TrackStateRecord first;
    first.voxel = {1, 2, 3};
    first.state.position = {-1.0, 0.2, 0.0};
    first.state.momentum = {1.0, 0.0, 0.0};
    first.cluster.cluster_minus_voxel_center = {0.0, 0.2, 0.0};

    cpm::TrackStateRecord second;
    second.voxel = {1, 2, 3};
    second.state.position = {0.1, -1.0, 0.0};
    second.state.momentum = {0.0, 1.0, 0.0};
    second.cluster.cluster_minus_voxel_center = {0.1, 0.0, 0.0};

    cpm::VoxelContainer container;
    container.add(first);
    container.add(second);

    assert(container.voxel_count() == 1);
    assert(container.record_count() == 2);

    const auto records = container.find({1, 2, 3});
    assert(records != nullptr);
    assert(records->size() == 2);

    const auto result = cpm::computeVoxelCenterPoCA(records->at(0), records->at(1));
    assert(result.valid);
    assert(near(result.midpoint.x, 0.0));
    assert(near(result.midpoint.y, 0.0));
    assert(near(result.midpoint.z, 0.0));
  }

  return 0;
}
