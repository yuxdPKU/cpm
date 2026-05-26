#ifndef CPM_CPMVOXELCONTAINER_H
#define CPM_CPMVOXELCONTAINER_H

#include "CPMRecord.h"

#include <iostream>

#if __has_include(<phool/PHObject.h>)
#include <phool/PHObject.h>
#elif __has_include(<TObject.h>)
#include <TObject.h>
#define CPM_FALLBACK_PHOBJECT_WITH_ROOT 1
class PHObject : public TObject
{
 public:
  ~PHObject() override = default;
  virtual void identify(std::ostream& = std::cout) const {}
  virtual void Reset() {}
  virtual int isValid() const { return 0; }

 private:
  ClassDefOverride(PHObject, 0)
};
#else
#ifndef ClassDefOverride
#define ClassDefOverride(name, version)
#endif
class PHObject
{
 public:
  virtual ~PHObject() = default;
  virtual void identify(std::ostream& = std::cout) const {}
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
};
#endif

#include <cstddef>
#include <map>
#include <vector>

class CPMVoxelContainer : public PHObject
{
 public:
  using Record = TrackStateRecord;
  using RecordVector = std::vector<Record>;
  using Map = std::map<VoxelId, RecordVector>;
  using const_iterator = Map::const_iterator;

  struct Grid
  {
    int phi_bins = 0;
    int r_bins = 0;
    int z_bins = 0;
    double phi_min = 0.0;
    double phi_max = 6.28318530717958647692;
    double r_min = 0.0;
    double r_max = 0.0;
    double z_min = 0.0;
    double z_max = 0.0;

    [[nodiscard]] bool valid() const
    {
      return phi_bins > 0 && r_bins > 0 && z_bins > 0 &&
             phi_min < phi_max && r_min < r_max && z_min < z_max;
    }
  };

  ~CPMVoxelContainer() override = default;

  void identify(std::ostream& out = std::cout) const override
  {
    out << "CPMVoxelContainer base class" << std::endl;
  }

  void Reset() override {}
  int isValid() const override { return 0; }

  virtual void set_grid(
      int /*phi_bins*/,
      int /*r_bins*/,
      int /*z_bins*/,
      double /*r_min*/,
      double /*r_max*/,
      double /*z_min*/,
      double /*z_max*/,
      double /*phi_min*/ = 0.0,
      double /*phi_max*/ = 6.28318530717958647692)
  {
  }

  [[nodiscard]] virtual const Grid& grid() const
  {
    static const Grid empty_grid;
    return empty_grid;
  }

  virtual void add(Record /*record*/) {}
  virtual bool add(const CPMVoxelContainer& /*other*/) { return false; }

  [[nodiscard]] virtual const RecordVector* find(const VoxelId& /*voxel*/) const
  {
    return nullptr;
  }

  [[nodiscard]] virtual const RecordVector* find_by_index(
      int iphi,
      int ir,
      int iz) const
  {
    return find({iphi, ir, iz});
  }

  [[nodiscard]] virtual const RecordVector* find_by_position(
      double /*phi*/,
      double /*radius*/,
      double /*z*/) const
  {
    return nullptr;
  }

  [[nodiscard]] virtual const RecordVector* find_by_cartesian(
      double /*x*/,
      double /*y*/,
      double /*z*/) const
  {
    return nullptr;
  }

  [[nodiscard]] virtual bool get_voxel_id(
      double /*phi*/,
      double /*radius*/,
      double /*z*/,
      VoxelId& /*voxel*/) const
  {
    return false;
  }

  [[nodiscard]] virtual std::size_t record_count(const VoxelId& voxel) const
  {
    const auto* records = find(voxel);
    return records ? records->size() : 0;
  }

  [[nodiscard]] virtual std::size_t voxel_count() const { return 0; }
  [[nodiscard]] virtual std::size_t record_count() const { return 0; }
  [[nodiscard]] virtual bool empty() const { return true; }
  [[nodiscard]] virtual const Map& records() const
  {
    static const Map empty_map;
    return empty_map;
  }

  virtual void sort_records() {}
  virtual void clear() {}

 private:
  ClassDefOverride(CPMVoxelContainer, 1)
};

#endif
