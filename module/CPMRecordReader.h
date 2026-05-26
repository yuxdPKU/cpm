#ifndef CPM_CPMRECORDREADER_H
#define CPM_CPMRECORDREADER_H

#include <CPMVoxelContainer.h>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

class CPMRecordReader
{
 public:
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
    Vector3 voxel_center;
    Vector3 offset;
    Vector3 state_position;
    Vector3 state_momentum;
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
    record.voxel_center = {
        input.cluster.voxel_center.x,
        input.cluster.voxel_center.y,
        input.cluster.voxel_center.z};
    record.offset = {
        input.cluster.cluster_minus_voxel_center.x,
        input.cluster.cluster_minus_voxel_center.y,
        input.cluster.cluster_minus_voxel_center.z};
    record.state_position = {
        input.state.position.x,
        input.state.position.y,
        input.state.position.z};
    record.state_momentum = {
        input.state.momentum.x,
        input.state.momentum.y,
        input.state.momentum.z};
    return record;
  }

  static bool has_branch(TChain& chain, const char* name)
  {
    if (chain.GetEntries() <= 0)
    {
      return chain.GetBranch(name) != nullptr;
    }
    chain.LoadTree(0);
    return chain.GetBranch(name) != nullptr;
  }

  template <class T>
  static T vector_value(const std::vector<T>* values, const std::size_t index, const T& fallback)
  {
    return values && index < values->size() ? values->at(index) : fallback;
  }

  template <class Callback>
  static Long64_t read_flat_records(TChain& chain, Callback callback)
  {
    std::string* cluster_source = nullptr;
    std::string* track_source = nullptr;
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
    float pt = std::numeric_limits<float>::quiet_NaN();
    int iphi = -1;
    int ir = -1;
    int iz = -1;
    double voxel_x = std::numeric_limits<double>::quiet_NaN();
    double voxel_y = std::numeric_limits<double>::quiet_NaN();
    double voxel_z = std::numeric_limits<double>::quiet_NaN();
    double offset_x = std::numeric_limits<double>::quiet_NaN();
    double offset_y = std::numeric_limits<double>::quiet_NaN();
    double offset_z = std::numeric_limits<double>::quiet_NaN();
    double state_x = std::numeric_limits<double>::quiet_NaN();
    double state_y = std::numeric_limits<double>::quiet_NaN();
    double state_z = std::numeric_limits<double>::quiet_NaN();
    double state_px = std::numeric_limits<double>::quiet_NaN();
    double state_py = std::numeric_limits<double>::quiet_NaN();
    double state_pz = std::numeric_limits<double>::quiet_NaN();

    chain.SetBranchAddress("cluster_source", &cluster_source);
    chain.SetBranchAddress("track_source", &track_source);
    chain.SetBranchAddress("run", &run);
    chain.SetBranchAddress("segment", &segment);
    chain.SetBranchAddress("sync_event", &sync_event);
    chain.SetBranchAddress("event_sequence", &event_sequence);
    chain.SetBranchAddress("stream_event_ordinal", &stream_event_ordinal);
    chain.SetBranchAddress("track_id", &track_id);
    chain.SetBranchAddress("cluskey", &cluskey);
    if (has_branch(chain, "hitsetkey"))
    {
      chain.SetBranchAddress("hitsetkey", &hitsetkey);
    }
    if (has_branch(chain, "subsurfkey"))
    {
      chain.SetBranchAddress("subsurfkey", &subsurfkey);
    }
    chain.SetBranchAddress("charge", &charge);
    chain.SetBranchAddress("pt", &pt);
    chain.SetBranchAddress("iphi", &iphi);
    chain.SetBranchAddress("ir", &ir);
    chain.SetBranchAddress("iz", &iz);
    chain.SetBranchAddress("voxel_x", &voxel_x);
    chain.SetBranchAddress("voxel_y", &voxel_y);
    chain.SetBranchAddress("voxel_z", &voxel_z);
    chain.SetBranchAddress("offset_x", &offset_x);
    chain.SetBranchAddress("offset_y", &offset_y);
    chain.SetBranchAddress("offset_z", &offset_z);
    chain.SetBranchAddress("state_x", &state_x);
    chain.SetBranchAddress("state_y", &state_y);
    chain.SetBranchAddress("state_z", &state_z);
    chain.SetBranchAddress("state_px", &state_px);
    chain.SetBranchAddress("state_py", &state_py);
    chain.SetBranchAddress("state_pz", &state_pz);

    const auto entries = chain.GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry)
    {
      chain.GetEntry(entry);

      Record record;
      record.entry = entry;
      record.cluster_source = cluster_source ? *cluster_source : "";
      record.track_source = track_source ? *track_source : "";
      record.run = run;
      record.segment = segment;
      record.sync_event = sync_event;
      record.event_sequence = event_sequence;
      record.stream_event_ordinal = stream_event_ordinal;
      record.track_id = track_id;
      record.cluskey = cluskey;
      record.hitsetkey = hitsetkey;
      record.subsurfkey = subsurfkey;
      record.charge = charge;
      record.pt = pt;
      record.iphi = iphi;
      record.ir = ir;
      record.iz = iz;
      record.voxel_center = {voxel_x, voxel_y, voxel_z};
      record.offset = {offset_x, offset_y, offset_z};
      record.state_position = {state_x, state_y, state_z};
      record.state_momentum = {state_px, state_py, state_pz};
      callback(record);
    }
    return entries;
  }

  template <class Callback>
  static Long64_t read_voxel_records(TChain& chain, Callback callback)
  {
    int iphi = -1;
    int ir = -1;
    int iz = -1;
    unsigned int entry_count = 0;
    double voxel_x = std::numeric_limits<double>::quiet_NaN();
    double voxel_y = std::numeric_limits<double>::quiet_NaN();
    double voxel_z = std::numeric_limits<double>::quiet_NaN();

    std::vector<std::string>* cluster_source = nullptr;
    std::vector<std::string>* track_source = nullptr;
    std::vector<int>* run = nullptr;
    std::vector<int>* segment = nullptr;
    std::vector<int>* sync_event = nullptr;
    std::vector<int>* event_sequence = nullptr;
    std::vector<unsigned long long>* stream_event_ordinal = nullptr;
    std::vector<unsigned int>* track_id = nullptr;
    std::vector<unsigned long long>* cluskey = nullptr;
    std::vector<unsigned long long>* hitsetkey = nullptr;
    std::vector<unsigned int>* subsurfkey = nullptr;
    std::vector<int>* charge = nullptr;
    std::vector<float>* pt = nullptr;
    std::vector<double>* offset_x = nullptr;
    std::vector<double>* offset_y = nullptr;
    std::vector<double>* offset_z = nullptr;
    std::vector<double>* state_x = nullptr;
    std::vector<double>* state_y = nullptr;
    std::vector<double>* state_z = nullptr;
    std::vector<double>* state_px = nullptr;
    std::vector<double>* state_py = nullptr;
    std::vector<double>* state_pz = nullptr;

    chain.SetBranchAddress("iphi", &iphi);
    chain.SetBranchAddress("ir", &ir);
    chain.SetBranchAddress("iz", &iz);
    chain.SetBranchAddress("entry_count", &entry_count);
    chain.SetBranchAddress("voxel_x", &voxel_x);
    chain.SetBranchAddress("voxel_y", &voxel_y);
    chain.SetBranchAddress("voxel_z", &voxel_z);
    chain.SetBranchAddress("cluster_source", &cluster_source);
    chain.SetBranchAddress("track_source", &track_source);
    chain.SetBranchAddress("run", &run);
    chain.SetBranchAddress("segment", &segment);
    chain.SetBranchAddress("sync_event", &sync_event);
    chain.SetBranchAddress("event_sequence", &event_sequence);
    chain.SetBranchAddress("stream_event_ordinal", &stream_event_ordinal);
    chain.SetBranchAddress("track_id", &track_id);
    chain.SetBranchAddress("cluskey", &cluskey);
    chain.SetBranchAddress("hitsetkey", &hitsetkey);
    chain.SetBranchAddress("subsurfkey", &subsurfkey);
    chain.SetBranchAddress("charge", &charge);
    chain.SetBranchAddress("pt", &pt);
    chain.SetBranchAddress("offset_x", &offset_x);
    chain.SetBranchAddress("offset_y", &offset_y);
    chain.SetBranchAddress("offset_z", &offset_z);
    chain.SetBranchAddress("state_x", &state_x);
    chain.SetBranchAddress("state_y", &state_y);
    chain.SetBranchAddress("state_z", &state_z);
    chain.SetBranchAddress("state_px", &state_px);
    chain.SetBranchAddress("state_py", &state_py);
    chain.SetBranchAddress("state_pz", &state_pz);

    Long64_t records_read = 0;
    const auto entries = chain.GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry)
    {
      chain.GetEntry(entry);
      const std::size_t records =
          track_id ? track_id->size() : static_cast<std::size_t>(entry_count);
      for (std::size_t index = 0; index < records; ++index)
      {
        Record record;
        record.entry = records_read;
        record.cluster_source = vector_value(cluster_source, index, std::string());
        record.track_source = vector_value(track_source, index, std::string());
        record.run = vector_value(run, index, -1);
        record.segment = vector_value(segment, index, -1);
        record.sync_event = vector_value(sync_event, index, -1);
        record.event_sequence = vector_value(event_sequence, index, -1);
        record.stream_event_ordinal =
            vector_value(stream_event_ordinal, index, 0ULL);
        record.track_id = vector_value(track_id, index, 0U);
        record.cluskey = vector_value(cluskey, index, 0ULL);
        record.hitsetkey = vector_value(hitsetkey, index, 0ULL);
        record.subsurfkey = vector_value(subsurfkey, index, 0U);
        record.charge = vector_value(charge, index, 0);
        record.pt = vector_value(pt, index, std::numeric_limits<float>::quiet_NaN());
        record.iphi = iphi;
        record.ir = ir;
        record.iz = iz;
        record.voxel_center = {voxel_x, voxel_y, voxel_z};
        record.offset = {
            vector_value(offset_x, index, std::numeric_limits<double>::quiet_NaN()),
            vector_value(offset_y, index, std::numeric_limits<double>::quiet_NaN()),
            vector_value(offset_z, index, std::numeric_limits<double>::quiet_NaN())};
        record.state_position = {
            vector_value(state_x, index, std::numeric_limits<double>::quiet_NaN()),
            vector_value(state_y, index, std::numeric_limits<double>::quiet_NaN()),
            vector_value(state_z, index, std::numeric_limits<double>::quiet_NaN())};
        record.state_momentum = {
            vector_value(state_px, index, std::numeric_limits<double>::quiet_NaN()),
            vector_value(state_py, index, std::numeric_limits<double>::quiet_NaN()),
            vector_value(state_pz, index, std::numeric_limits<double>::quiet_NaN())};
        callback(record);
        ++records_read;
      }
    }
    return records_read;
  }

  template <class Callback>
  static Long64_t read_records(TChain& chain, Callback callback)
  {
    return has_branch(chain, "entry_count") ?
           read_voxel_records(chain, callback) :
           read_flat_records(chain, callback);
  }

  template <class Callback>
  static Long64_t read_object_records(const std::string& filename, Callback callback)
  {
    TFile input(filename.c_str(), "READ");
    if (input.IsZombie())
    {
      return -1;
    }

    auto* container = dynamic_cast<CPMVoxelContainer*>(input.Get("cpm_records"));
    if (!container)
    {
      container = dynamic_cast<CPMVoxelContainer*>(input.Get("CPMVoxelContainer"));
    }
    if (!container)
    {
      return -1;
    }

    Long64_t records_read = 0;
    for (const auto& [voxel, records] : container->records())
    {
      for (const auto& record : records)
      {
        callback(make_record(record, voxel, records_read));
        ++records_read;
      }
    }
    return records_read;
  }

  template <class Callback>
  static Long64_t read_records(const std::vector<std::string>& input_files, Callback callback)
  {
    Long64_t records_read = 0;
    std::vector<std::string> flat_tree_files;
    for (const auto& file : input_files)
    {
      const auto object_records = read_object_records(file, callback);
      if (object_records >= 0)
      {
        records_read += object_records;
      }
      else
      {
        flat_tree_files.push_back(file);
      }
    }

    if (!flat_tree_files.empty())
    {
      TChain chain("cpm_records");
      for (const auto& file : flat_tree_files)
      {
        chain.Add(file.c_str());
      }
      records_read += read_records(chain, callback);
    }
    return records_read;
  }
};

#endif
