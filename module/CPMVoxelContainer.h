#ifndef CPM_CPMVOXELCONTAINER_H
#define CPM_CPMVOXELCONTAINER_H

#include "CPMRecord.h"

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace cpm
{
  class VoxelContainer
  {
   public:
    using Record = TrackStateRecord;
    using RecordVector = std::vector<Record>;
    using Map = std::map<VoxelId, RecordVector>;
    using const_iterator = Map::const_iterator;

    void add(Record record)
    {
      m_records[record.voxel].push_back(std::move(record));
    }

    [[nodiscard]] const RecordVector* find(const VoxelId& voxel) const
    {
      const auto iter = m_records.find(voxel);
      return iter == m_records.end() ? nullptr : &iter->second;
    }

    [[nodiscard]] std::size_t voxel_count() const { return m_records.size(); }

    [[nodiscard]] std::size_t record_count() const
    {
      std::size_t out = 0;
      for (const auto& [voxel, records] : m_records)
      {
        (void) voxel;
        out += records.size();
      }
      return out;
    }

    [[nodiscard]] bool empty() const { return m_records.empty(); }

    [[nodiscard]] const_iterator begin() const { return m_records.begin(); }
    [[nodiscard]] const_iterator end() const { return m_records.end(); }

    void clear() { m_records.clear(); }

   private:
    Map m_records;
  };
}

#endif
