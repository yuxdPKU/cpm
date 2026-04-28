#include "CPMLocalLinePoCA.h"

#include <cmath>

namespace cpm
{
  LocalLinePoCAResult computeLocalLinePoCA(
      const Vector3& point_a,
      const Vector3& direction_a,
      const Vector3& point_b,
      const Vector3& direction_b,
      const LocalLinePoCAOptions& options)
  {
    LocalLinePoCAResult result;

    const Vector3 u = unit(direction_a);
    const Vector3 v = unit(direction_b);
    const Vector3 w0 = point_a - point_b;

    const double a = dot(u, u);
    const double b = dot(u, v);
    const double c = dot(v, v);
    const double d = dot(u, w0);
    const double e = dot(v, w0);
    const double denom = a * c - b * b;

    if (!(a > 0.0 && c > 0.0) || denom <= options.min_sin_angle * options.min_sin_angle)
    {
      return result;
    }

    result.s = (b * e - c * d) / denom;
    result.t = (a * e - b * d) / denom;
    result.point_a = point_a + u * result.s;
    result.point_b = point_b + v * result.t;
    result.midpoint = (result.point_a + result.point_b) * 0.5;
    result.dca = (result.point_a - result.point_b).norm();
    result.valid = std::isfinite(result.dca);

    return result;
  }

  LocalLinePoCAResult computeVoxelCenterPoCA(
      const TrackStateRecord& record_a,
      const TrackStateRecord& record_b,
      const LocalLinePoCAOptions& options)
  {
    const Vector3 point_a = record_a.state.position - record_a.cluster.cluster_minus_voxel_center;
    const Vector3 point_b = record_b.state.position - record_b.cluster.cluster_minus_voxel_center;

    return computeLocalLinePoCA(
        point_a,
        record_a.state.momentum,
        point_b,
        record_b.state.momentum,
        options);
  }
}
