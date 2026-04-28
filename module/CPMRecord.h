#ifndef CPM_CPMRECORD_H
#define CPM_CPMRECORD_H

#include "CPMTypes.h"

#include <array>
#include <cstdint>
#include <limits>
#include <string>

namespace cpm
{
  struct EventReference
  {
    std::string cluster_source;
    std::string track_source;
    std::int32_t run = -1;
    std::int32_t segment = -1;
    std::int32_t sync_event = -1;
    std::int32_t event_sequence = -1;
    std::uint64_t stream_event_ordinal = 0;
  };

  struct TrackReference
  {
    std::string track_map_name = "SvtxSiliconMMTrackMap";
    TrackId track_id = InvalidTrackId;
    ClusterKey state_cluskey = InvalidClusterKey;
  };

  struct TrackSummary
  {
    std::int32_t charge = 0;
    float pt = std::numeric_limits<float>::quiet_NaN();
    float quality = std::numeric_limits<float>::quiet_NaN();
    std::uint16_t n_mvtx = 0;
    std::uint16_t n_intt = 0;
    std::uint16_t n_tpc = 0;
    std::uint16_t n_tpot = 0;
    std::uint16_t n_mvtx_states = 0;
    std::uint16_t n_intt_states = 0;
    std::uint16_t n_tpc_states = 0;
    std::uint16_t n_tpot_states = 0;
  };

  struct ClusterReference
  {
    ClusterKey cluskey = InvalidClusterKey;
    HitSetKey hitsetkey = InvalidHitSetKey;
    SubSurfKey subsurfkey = InvalidSubSurfKey;
    std::uint16_t layer = 0;
    std::uint16_t side = 0;
  };

  struct VoxelId
  {
    std::int32_t iphi = -1;
    std::int32_t ir = -1;
    std::int32_t iz = -1;

    [[nodiscard]] bool valid() const
    {
      return iphi >= 0 && ir >= 0 && iz >= 0;
    }

    friend bool operator==(const VoxelId& lhs, const VoxelId& rhs)
    {
      return lhs.iphi == rhs.iphi && lhs.ir == rhs.ir && lhs.iz == rhs.iz;
    }

    friend bool operator<(const VoxelId& lhs, const VoxelId& rhs)
    {
      if (lhs.iphi != rhs.iphi) return lhs.iphi < rhs.iphi;
      if (lhs.ir != rhs.ir) return lhs.ir < rhs.ir;
      return lhs.iz < rhs.iz;
    }
  };

  struct ClusterSnapshot
  {
    Vector3 corrected_position;
    Vector3 voxel_center;
    Vector3 cluster_minus_voxel_center;
  };

  struct TrackStateSnapshot
  {
    double pathlength = std::numeric_limits<double>::quiet_NaN();
    double local_x = std::numeric_limits<double>::quiet_NaN();
    double local_y = std::numeric_limits<double>::quiet_NaN();
    Vector3 position;
    Vector3 momentum;
    Matrix6 covariance{};
  };

  struct SurfaceSnapshot
  {
    std::uint64_t geometry_id = 0;
    std::uint64_t volume = 0;
    std::uint64_t layer = 0;
    std::uint64_t sensitive = 0;
    std::uint64_t approach = 0;
    std::uint64_t boundary = 0;
    Vector3 center;
    Vector3 normal;
    Vector3 local_x_axis;
    Vector3 local_y_axis;
  };

  struct SelectionFlags
  {
    bool has_crossing = false;
    bool passes_cm = false;
    bool passes_tpot = false;
    bool passes_track_quality = false;
    bool passes_geometry = false;
  };

  struct TrackStateRecord
  {
    EventReference event_ref;
    TrackReference track_ref;
    TrackSummary track;
    ClusterReference cluster_ref;
    VoxelId voxel;
    ClusterSnapshot cluster;
    TrackStateSnapshot state;
    SurfaceSnapshot surface;
    SelectionFlags selection;
  };
}

#endif
