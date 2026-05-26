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
  using const_iterator = Map::const_iterator;

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

  void add(Record record) override;
  bool add(const CPMVoxelContainer& other) override;

  [[nodiscard]] const RecordVector* find(const VoxelId& voxel) const override;
  [[nodiscard]] const RecordVector* find_by_index(int iphi, int ir, int iz) const override;
  [[nodiscard]] const RecordVector* find_by_position(double phi, double radius, double z) const override;
  [[nodiscard]] const RecordVector* find_by_cartesian(double x, double y, double z) const override;
  [[nodiscard]] bool get_voxel_id(double phi, double radius, double z, VoxelId& voxel) const override;

  [[nodiscard]] std::size_t record_count(const VoxelId& voxel) const override;
  [[nodiscard]] std::size_t voxel_count() const override { return m_records.size(); }
  [[nodiscard]] std::size_t record_count() const override;
  [[nodiscard]] bool empty() const override { return m_records.empty(); }
  [[nodiscard]] const Map& records() const override { return m_records; }

  [[nodiscard]] const_iterator begin() const { return m_records.begin(); }
  [[nodiscard]] const_iterator end() const { return m_records.end(); }

  void sort_records() override;
  void clear() override { m_records.clear(); }

 private:
  static bool has_event_order(const Record& record);
  static std::tuple<int, int, int, unsigned long long, int, TrackId, ClusterKey>
  event_order_key(const Record& record);
  static void stable_sort_records(RecordVector& records);

  Grid m_grid;
  Map m_records;

  ClassDefOverride(CPMVoxelContainerv1, 1)
};

#endif
