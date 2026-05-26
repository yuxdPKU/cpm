#include "CPMVoxelContainerv1.h"

#include <iostream>
#include <utility>

void CPMVoxelContainerv1::identify(std::ostream& out) const
{
  out << "CPMVoxelContainerv1 with " << record_count()
      << " records in " << voxel_count() << " voxels" << std::endl;
}

void CPMVoxelContainerv1::Reset()
{
  m_records.clear();
}

int CPMVoxelContainerv1::isValid() const
{
  return m_grid.valid() ? 1 : 0;
}

void CPMVoxelContainerv1::set_grid(
    const int phi_bins,
    const int r_bins,
    const int z_bins,
    const double r_min,
    const double r_max,
    const double z_min,
    const double z_max,
    const double phi_min,
    const double phi_max)
{
  m_grid.phi_bins = phi_bins;
  m_grid.r_bins = r_bins;
  m_grid.z_bins = z_bins;
  m_grid.phi_min = phi_min;
  m_grid.phi_max = phi_max;
  m_grid.r_min = r_min;
  m_grid.r_max = r_max;
  m_grid.z_min = z_min;
  m_grid.z_max = z_max;
}

void CPMVoxelContainerv1::add(TrackStateRecord record)
{
  m_records[record.voxel].push_back(std::move(record));
}

bool CPMVoxelContainerv1::add(const CPMVoxelContainer& other)
{
  const auto& other_grid = other.grid();
  if (!m_grid.valid() && other_grid.valid())
  {
    m_grid = other_grid;
  }
  else if (m_grid.valid() && other_grid.valid() &&
           (m_grid.phi_bins != other_grid.phi_bins ||
            m_grid.r_bins != other_grid.r_bins ||
            m_grid.z_bins != other_grid.z_bins ||
            m_grid.phi_min != other_grid.phi_min ||
            m_grid.phi_max != other_grid.phi_max ||
            m_grid.r_min != other_grid.r_min ||
            m_grid.r_max != other_grid.r_max ||
            m_grid.z_min != other_grid.z_min ||
            m_grid.z_max != other_grid.z_max))
  {
    std::cout << "CPMVoxelContainerv1::add - inconsistent grid" << std::endl;
    return false;
  }

  for (const auto& [voxel, records] : other.records())
  {
    auto& out = m_records[voxel];
    out.insert(out.end(), records.begin(), records.end());
  }
  sort_records();
  return true;
}

const std::vector<TrackStateRecord>* CPMVoxelContainerv1::find(const VoxelId& voxel) const
{
  const auto iter = m_records.find(voxel);
  return iter == m_records.end() ? nullptr : &iter->second;
}

const std::vector<TrackStateRecord>* CPMVoxelContainerv1::find_by_index(
    const int iphi,
    const int ir,
    const int iz) const
{
  return find({iphi, ir, iz});
}

const std::vector<TrackStateRecord>* CPMVoxelContainerv1::find_by_position(
    double phi,
    const double radius,
    const double z) const
{
  VoxelId voxel;
  return get_voxel_id(phi, radius, z, voxel) ? find(voxel) : nullptr;
}

const std::vector<TrackStateRecord>* CPMVoxelContainerv1::find_by_cartesian(
    const double x,
    const double y,
    const double z) const
{
  double phi = std::atan2(y, x);
  if (phi < 0.0)
  {
    phi += m_grid.phi_max - m_grid.phi_min;
  }
  const double radius = std::sqrt(x * x + y * y);
  return find_by_position(phi, radius, z);
}

bool CPMVoxelContainerv1::get_voxel_id(
    double phi,
    const double radius,
    const double z,
    VoxelId& voxel) const
{
  if (!m_grid.valid())
  {
    return false;
  }

  const double phi_width = m_grid.phi_max - m_grid.phi_min;
  while (phi < m_grid.phi_min)
  {
    phi += phi_width;
  }
  while (phi >= m_grid.phi_max)
  {
    phi -= phi_width;
  }

  if (phi < m_grid.phi_min || phi >= m_grid.phi_max ||
      radius < m_grid.r_min || radius >= m_grid.r_max ||
      z < m_grid.z_min || z >= m_grid.z_max)
  {
    return false;
  }

  voxel.iphi = static_cast<int>(
      m_grid.phi_bins * (phi - m_grid.phi_min) / phi_width);
  voxel.ir = static_cast<int>(
      m_grid.r_bins * (radius - m_grid.r_min) / (m_grid.r_max - m_grid.r_min));
  voxel.iz = static_cast<int>(
      m_grid.z_bins * (z - m_grid.z_min) / (m_grid.z_max - m_grid.z_min));
  return voxel.valid();
}

std::size_t CPMVoxelContainerv1::record_count(const VoxelId& voxel) const
{
  const auto* records = find(voxel);
  return records ? records->size() : 0;
}

std::size_t CPMVoxelContainerv1::record_count() const
{
  std::size_t out = 0;
  for (const auto& [voxel, records] : m_records)
  {
    (void) voxel;
    out += records.size();
  }
  return out;
}

void CPMVoxelContainerv1::sort_records()
{
  for (auto& [voxel, records] : m_records)
  {
    (void) voxel;
    stable_sort_records(records);
  }
}

bool CPMVoxelContainerv1::has_event_order(const TrackStateRecord& record)
{
  return record.event_ref.run >= 0 &&
         (record.event_ref.event_sequence >= 0 ||
          record.event_ref.sync_event >= 0);
}

std::tuple<int, int, int, unsigned long long, int, unsigned int, TrkrDefs::cluskey>
CPMVoxelContainerv1::event_order_key(const TrackStateRecord& record)
{
  const int event_sequence =
      record.event_ref.event_sequence >= 0 ?
      record.event_ref.event_sequence :
      std::numeric_limits<int>::max();
  const int sync_event =
      record.event_ref.sync_event >= 0 ?
      record.event_ref.sync_event :
      std::numeric_limits<int>::max();
  return {
      record.event_ref.run,
      event_sequence,
      sync_event,
      record.event_ref.stream_event_ordinal,
      record.event_ref.segment,
      record.track_ref.track_id,
      record.cluster_ref.cluskey};
}

void CPMVoxelContainerv1::stable_sort_records(std::vector<TrackStateRecord>& records)
{
  std::stable_sort(
      records.begin(),
      records.end(),
      [](const TrackStateRecord& lhs, const TrackStateRecord& rhs)
      {
        const bool lhs_has_order = has_event_order(lhs);
        const bool rhs_has_order = has_event_order(rhs);
        if (lhs_has_order != rhs_has_order)
        {
          return lhs_has_order;
        }
        if (!lhs_has_order)
        {
          return false;
        }
        return event_order_key(lhs) < event_order_key(rhs);
      });
}
