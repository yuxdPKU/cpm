#ifndef CPM_CPMTYPES_H
#define CPM_CPMTYPES_H

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace cpm
{
  using ClusterKey = std::uint64_t;
  using HitSetKey = std::uint64_t;
  using SubSurfKey = std::uint32_t;
  using TrackId = std::uint32_t;

  inline constexpr ClusterKey InvalidClusterKey = std::numeric_limits<ClusterKey>::max();
  inline constexpr HitSetKey InvalidHitSetKey = std::numeric_limits<HitSetKey>::max();
  inline constexpr SubSurfKey InvalidSubSurfKey = std::numeric_limits<SubSurfKey>::max();
  inline constexpr TrackId InvalidTrackId = std::numeric_limits<TrackId>::max();

  struct Vector3
  {
    double x = std::numeric_limits<double>::quiet_NaN();
    double y = std::numeric_limits<double>::quiet_NaN();
    double z = std::numeric_limits<double>::quiet_NaN();

    [[nodiscard]] double norm2() const { return x * x + y * y + z * z; }
    [[nodiscard]] double norm() const { return std::sqrt(norm2()); }
  };

  inline Vector3 operator+(const Vector3& lhs, const Vector3& rhs)
  {
    return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
  }

  inline Vector3 operator-(const Vector3& lhs, const Vector3& rhs)
  {
    return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
  }

  inline Vector3 operator*(const Vector3& lhs, const double scale)
  {
    return {lhs.x * scale, lhs.y * scale, lhs.z * scale};
  }

  inline Vector3 operator*(const double scale, const Vector3& rhs)
  {
    return rhs * scale;
  }

  inline Vector3 operator/(const Vector3& lhs, const double scale)
  {
    return {lhs.x / scale, lhs.y / scale, lhs.z / scale};
  }

  inline double dot(const Vector3& lhs, const Vector3& rhs)
  {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
  }

  inline Vector3 unit(const Vector3& value)
  {
    const double length = value.norm();
    return (length > 0.0) ? value / length : Vector3{};
  }

  using Matrix6 = std::array<double, 36>;
}

#endif
