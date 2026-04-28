#ifndef CPM_CPMLOCALLINEPOCA_H
#define CPM_CPMLOCALLINEPOCA_H

#include "CPMRecord.h"
#include "CPMTypes.h"

namespace cpm
{
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
    Vector3 point_a;
    Vector3 point_b;
    Vector3 midpoint;
  };

  LocalLinePoCAResult computeLocalLinePoCA(
      const Vector3& point_a,
      const Vector3& direction_a,
      const Vector3& point_b,
      const Vector3& direction_b,
      const LocalLinePoCAOptions& options = {});

  LocalLinePoCAResult computeVoxelCenterPoCA(
      const TrackStateRecord& record_a,
      const TrackStateRecord& record_b,
      const LocalLinePoCAOptions& options = {});
}

#endif
