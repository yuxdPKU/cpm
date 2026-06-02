#ifndef CPM_CPMVOXELCONTAINERV1_H
#define CPM_CPMVOXELCONTAINERV1_H

#include "CPMVoxelContainer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

class CPMVoxelContainerv1 : public CPMVoxelContainer
{
 public:
  ~CPMVoxelContainerv1() override = default;

  void identify(std::ostream& out = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  void set_grid(
      int phi_bins,
      int r_bins,
      int z_bins,
      double r_min,
      double r_max,
      double z_min,
      double z_max,
      double phi_min = 0.0,
      double phi_max = 6.28318530717958647692) override;

  [[nodiscard]] const Grid& grid() const override { return m_grid; }

  void add(TrackStateRecord record) override;
  bool add(const CPMVoxelContainer& other) override;

  [[nodiscard]] const std::vector<TrackStateRecord>* find(const VoxelId& voxel) const override;
  [[nodiscard]] const std::vector<TrackStateRecord>* find_by_index(int iphi, int ir, int iz) const override;
  [[nodiscard]] const std::vector<TrackStateRecord>* find_by_position(double phi, double radius, double z) const override;
  [[nodiscard]] const std::vector<TrackStateRecord>* find_by_cartesian(double x, double y, double z) const override;
  [[nodiscard]] bool get_voxel_id(double phi, double radius, double z, VoxelId& voxel) const override;

  [[nodiscard]] std::size_t record_count(const VoxelId& voxel) const override;
  [[nodiscard]] std::size_t voxel_count() const override { return m_records.size(); }
  [[nodiscard]] std::size_t record_count() const override;
  [[nodiscard]] bool empty() const override { return m_records.empty(); }
  [[nodiscard]] const std::map<VoxelId, std::vector<TrackStateRecord>>& records() const override { return m_records; }

  [[nodiscard]] std::map<VoxelId, std::vector<TrackStateRecord>>::const_iterator begin() const { return m_records.begin(); }
  [[nodiscard]] std::map<VoxelId, std::vector<TrackStateRecord>>::const_iterator end() const { return m_records.end(); }

  std::vector<TrackStateRecord> take_records(const VoxelId& voxel);
  void sort_records() override;
  void clear() override { m_records.clear(); }

 private:
  static bool has_event_order(const TrackStateRecord& record);
  static std::tuple<int, int, int, unsigned long long, int, unsigned int, TrkrDefs::cluskey>
  event_order_key(const TrackStateRecord& record);
  static void stable_sort_records(std::vector<TrackStateRecord>& records);

  Grid m_grid;
  std::map<VoxelId, std::vector<TrackStateRecord>> m_records;

  ClassDefOverride(CPMVoxelContainerv1, 2)
};

#endif
