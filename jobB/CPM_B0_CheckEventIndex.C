/*
 * CPM Job B0 event-index QA.
 *
 * This macro checks the index produced by CPM_B0_BuildEventIndex.C before a
 * later mini-DST rehydration pass uses it for sequential event readback.
 *
 * Example:
 *
 *   root -l -b -q 'macro/CPM_B0_CheckEventIndex.C("CPM_B0_event_index.root")'
 */

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <map>
#include <set>
#include <string>

bool CPM_B0_CheckEventIndex(
    const std::string& index_file = "CPM_B0_event_index.root")
{
  TFile input(index_file.c_str(), "READ");
  if (input.IsZombie())
  {
    std::cerr << "CPM_B0_CheckEventIndex - cannot open " << index_file << std::endl;
    return false;
  }

  auto* events = dynamic_cast<TTree*>(input.Get("cpm_event_requests"));
  auto* objects = dynamic_cast<TTree*>(input.Get("cpm_object_requests"));
  if (!events || !objects)
  {
    std::cerr << "CPM_B0_CheckEventIndex - missing cpm_event_requests or cpm_object_requests" << std::endl;
    return false;
  }

  unsigned int event_index = 0;
  std::string* track_source = nullptr;
  std::string* cluster_source = nullptr;
  int run = -1;
  int segment = -1;
  int sync_event = -1;
  int event_sequence = -1;
  unsigned long long stream_event_ordinal = 0;
  unsigned int object_count = 0;

  events->SetBranchAddress("event_index", &event_index);
  events->SetBranchAddress("track_source", &track_source);
  events->SetBranchAddress("cluster_source", &cluster_source);
  events->SetBranchAddress("run", &run);
  events->SetBranchAddress("segment", &segment);
  events->SetBranchAddress("sync_event", &sync_event);
  events->SetBranchAddress("event_sequence", &event_sequence);
  events->SetBranchAddress("stream_event_ordinal", &stream_event_ordinal);
  events->SetBranchAddress("object_count", &object_count);

  std::set<unsigned int> event_indices;
  std::map<unsigned int, unsigned int> declared_counts;
  std::map<std::string, unsigned int> events_by_track_source;
  unsigned int empty_track_source_events = 0;
  unsigned int empty_cluster_source_events = 0;
  bool ok = true;

  const auto n_events = events->GetEntries();
  for (Long64_t entry = 0; entry < n_events; ++entry)
  {
    events->GetEntry(entry);

    if (!event_indices.insert(event_index).second)
    {
      std::cerr << "CPM_B0_CheckEventIndex - duplicate event_index " << event_index << std::endl;
      ok = false;
    }

    declared_counts[event_index] = object_count;
    const std::string track_source_value = track_source ? *track_source : "";
    const std::string cluster_source_value = cluster_source ? *cluster_source : "";
    events_by_track_source[track_source_value] += 1;
    if (track_source_value.empty())
    {
      ++empty_track_source_events;
    }
    if (cluster_source_value.empty())
    {
      ++empty_cluster_source_events;
    }

    (void) run;
    (void) segment;
    (void) sync_event;
    (void) event_sequence;
    (void) stream_event_ordinal;
  }

  for (unsigned int expected = 0; expected < static_cast<unsigned int>(n_events); ++expected)
  {
    if (event_indices.find(expected) == event_indices.end())
    {
      std::cerr << "CPM_B0_CheckEventIndex - missing sequential event_index " << expected << std::endl;
      ok = false;
    }
  }

  unsigned int object_event_index = 0;
  unsigned int track_id = 0;
  unsigned long long cluskey = 0;
  unsigned long long hitsetkey = 0;
  unsigned int subsurfkey = 0;
  int iphi = -1;
  int ir = -1;
  int iz = -1;

  objects->SetBranchAddress("event_index", &object_event_index);
  objects->SetBranchAddress("track_id", &track_id);
  objects->SetBranchAddress("cluskey", &cluskey);
  objects->SetBranchAddress("hitsetkey", &hitsetkey);
  objects->SetBranchAddress("subsurfkey", &subsurfkey);
  objects->SetBranchAddress("iphi", &iphi);
  objects->SetBranchAddress("ir", &ir);
  objects->SetBranchAddress("iz", &iz);

  std::map<unsigned int, unsigned int> actual_counts;
  unsigned int invalid_object_event = 0;
  unsigned int invalid_voxel = 0;
  unsigned int invalid_keys = 0;

  const auto n_objects = objects->GetEntries();
  for (Long64_t entry = 0; entry < n_objects; ++entry)
  {
    objects->GetEntry(entry);

    if (event_indices.find(object_event_index) == event_indices.end())
    {
      ++invalid_object_event;
      ok = false;
    }
    actual_counts[object_event_index] += 1;

    if (iphi < 0 || ir < 0 || iz < 0)
    {
      ++invalid_voxel;
      ok = false;
    }
    if (cluskey == 0 || hitsetkey == 0)
    {
      ++invalid_keys;
      ok = false;
    }

    (void) track_id;
    (void) subsurfkey;
  }

  for (const auto& [idx, declared] : declared_counts)
  {
    const auto actual_it = actual_counts.find(idx);
    const unsigned int actual = actual_it == actual_counts.end() ? 0 : actual_it->second;
    if (declared != actual)
    {
      std::cerr << "CPM_B0_CheckEventIndex - event_index " << idx
                << " declares " << declared << " objects but has " << actual << std::endl;
      ok = false;
    }
  }

  std::cout << "CPM_B0_CheckEventIndex - file: " << index_file << std::endl;
  std::cout << "CPM_B0_CheckEventIndex - events: " << n_events << std::endl;
  std::cout << "CPM_B0_CheckEventIndex - objects: " << n_objects << std::endl;
  std::cout << "CPM_B0_CheckEventIndex - empty track_source events: " << empty_track_source_events << std::endl;
  std::cout << "CPM_B0_CheckEventIndex - empty cluster_source events: " << empty_cluster_source_events << std::endl;
  std::cout << "CPM_B0_CheckEventIndex - invalid object event links: " << invalid_object_event << std::endl;
  std::cout << "CPM_B0_CheckEventIndex - invalid voxel indices: " << invalid_voxel << std::endl;
  std::cout << "CPM_B0_CheckEventIndex - invalid keys: " << invalid_keys << std::endl;
  std::cout << "CPM_B0_CheckEventIndex - events by track_source:" << std::endl;
  for (const auto& [source, count] : events_by_track_source)
  {
    std::cout << "  " << (source.empty() ? "<empty>" : source) << ": " << count << std::endl;
  }
  std::cout << "CPM_B0_CheckEventIndex - status: " << (ok ? "OK" : "FAILED") << std::endl;

  return ok;
}
