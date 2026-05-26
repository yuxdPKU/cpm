#ifndef CPM_JOBB_CPMVOXELCONTAINERREADER_H
#define CPM_JOBB_CPMVOXELCONTAINERREADER_H

#include <CPMVoxelContainer.h>

#include <TFile.h>
#include <TVector3.h>

#include <iostream>
#include <limits>
#include <string>
#include <vector>

class CPMVoxelContainerReader
{
 public:
  static constexpr const char* ObjectName = "CPMVoxelContainer";

  struct Record
  {
    Long64_t entry = -1;
    std::string cluster_source;
    std::string track_source;
    int run = -1;
    int segment = -1;
    int sync_event = -1;
    int event_sequence = -1;
    unsigned long long stream_event_ordinal = 0;
    unsigned int track_id = 0;
    unsigned long long cluskey = 0;
    unsigned long long hitsetkey = 0;
    unsigned int subsurfkey = 0;
    int charge = 0;
    double pt = std::numeric_limits<double>::quiet_NaN();
    int iphi = -1;
    int ir = -1;
    int iz = -1;
    TVector3 voxel_center;
    TVector3 offset;
    TVector3 state_position;
    TVector3 state_momentum;
  };

  static Record make_record(
      const TrackStateRecord& input,
      const VoxelId& voxel,
      const Long64_t entry)
  {
    Record record;
    record.entry = entry;
    record.cluster_source = input.event_ref.cluster_source;
    record.track_source = input.event_ref.track_source;
    record.run = input.event_ref.run;
    record.segment = input.event_ref.segment;
    record.sync_event = input.event_ref.sync_event;
    record.event_sequence = input.event_ref.event_sequence;
    record.stream_event_ordinal = input.event_ref.stream_event_ordinal;
    record.track_id = input.track_ref.track_id;
    record.cluskey = input.cluster_ref.cluskey;
    record.hitsetkey = input.cluster_ref.hitsetkey;
    record.subsurfkey = input.cluster_ref.subsurfkey;
    record.charge = input.track.charge;
    record.pt = input.track.pt;
    record.iphi = voxel.iphi;
    record.ir = voxel.ir;
    record.iz = voxel.iz;
    record.voxel_center = input.cluster.voxel_center;
    record.offset = input.cluster.cluster_minus_voxel_center;
    record.state_position = input.state.position;
    record.state_momentum = input.state.momentum;
    return record;
  }

  template <class Callback>
  static Long64_t read_records(
      const CPMVoxelContainer& container,
      Callback callback,
      const Long64_t first_entry = 0)
  {
    Long64_t records_read = 0;
    for (const auto& [voxel, records] : container.records())
    {
      for (const auto& record : records)
      {
        callback(make_record(record, voxel, first_entry + records_read));
        ++records_read;
      }
    }
    return records_read;
  }

  template <class Callback>
  static Long64_t read_records(
      const std::string& filename,
      Callback callback,
      const Long64_t first_entry = 0)
  {
    TFile input(filename.c_str(), "READ");
    if (input.IsZombie())
    {
      std::cerr << "CPMVoxelContainerReader::read_records - failed to open "
                << filename << std::endl;
      return -1;
    }

    auto* container = dynamic_cast<CPMVoxelContainer*>(input.Get(ObjectName));
    if (!container)
    {
      std::cerr << "CPMVoxelContainerReader::read_records - missing "
                << ObjectName << " in " << filename << std::endl;
      return -1;
    }

    return read_records(*container, callback, first_entry);
  }

  template <class Callback>
  static Long64_t read_records(
      const std::vector<std::string>& input_files,
      Callback callback)
  {
    Long64_t records_read = 0;
    for (const auto& file : input_files)
    {
      const auto file_records = read_records(file, callback, records_read);
      if (file_records < 0)
      {
        return -1;
      }
      records_read += file_records;
    }
    return records_read;
  }
};

#endif
