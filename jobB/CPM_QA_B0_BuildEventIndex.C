/*
 * CPM QA Job B0 skeleton.
 *
 * This macro reads one or more Job A CPM ROOT files, scans CPMVoxelContainer,
 * groups all referenced objects by event, and writes an event-ordered request
 * index. The next B0 step can use this index to read each mini-DST event once
 * and validate/rehydrate the requested
 * SvtxTrack/SvtxTrackState/TrkrCluster objects.
 *
 * Example:
 *
 *   root -l -b -q 'jobB/CPM_QA_B0_BuildEventIndex.C("jobA_CPMVoxelContainer.root",
 *                                                  "CPM_QA_B0_event_index.root")'
 */

#include "CPMVoxelContainerReader.h"

#include <TFile.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

R__LOAD_LIBRARY(libcpm.so)

namespace CPMB0
{
  struct EventKey
  {
    std::string track_source;
    std::string cluster_source;
    int run = -1;
    int segment = -1;
    int sync_event = -1;
    int event_sequence = -1;
    unsigned long long stream_event_ordinal = 0;

    bool operator<(const EventKey& rhs) const
    {
      return std::tie(track_source, cluster_source, run, segment,
                      sync_event, event_sequence, stream_event_ordinal) <
             std::tie(rhs.track_source, rhs.cluster_source, rhs.run, rhs.segment,
                      rhs.sync_event, rhs.event_sequence, rhs.stream_event_ordinal);
    }
  };

  struct ObjectRequest
  {
    unsigned int track_id = 0;
    unsigned long long cluskey = 0;
    unsigned long long hitsetkey = 0;
    unsigned int subsurfkey = 0;
    int iphi = -1;
    int ir = -1;
    int iz = -1;

    bool operator<(const ObjectRequest& rhs) const
    {
      return std::tie(track_id, cluskey, hitsetkey, subsurfkey, iphi, ir, iz) <
             std::tie(rhs.track_id, rhs.cluskey, rhs.hitsetkey,
                      rhs.subsurfkey, rhs.iphi, rhs.ir, rhs.iz);
    }
  };
}

void CPM_QA_B0_BuildEventIndex(
    const std::vector<std::string>& input_files,
    const std::string& output_file = "CPM_QA_B0_event_index.root")
{
  std::map<CPMB0::EventKey, std::set<CPMB0::ObjectRequest>> event_requests;

  const auto records_read = CPMVoxelContainerReader::read_records(
      input_files,
      [&](const CPMVoxelContainerReader::Record& input_record)
      {
    CPMB0::EventKey event_key;
    event_key.track_source = input_record.track_source;
    event_key.cluster_source = input_record.cluster_source;
    event_key.run = input_record.run;
    event_key.segment = input_record.segment;
    event_key.sync_event = input_record.sync_event;
    event_key.event_sequence = input_record.event_sequence;
    event_key.stream_event_ordinal = input_record.stream_event_ordinal;

    CPMB0::ObjectRequest request;
    request.track_id = input_record.track_id;
    request.cluskey = input_record.cluskey;
    request.hitsetkey = input_record.hitsetkey;
    request.subsurfkey = input_record.subsurfkey;
    request.iphi = input_record.iphi;
    request.ir = input_record.ir;
    request.iz = input_record.iz;
    event_requests[event_key].insert(request);
      });
  if (records_read < 0)
  {
    return;
  }

  TFile output(output_file.c_str(), "RECREATE");

  TTree events("cpm_event_requests", "Unique CPM events requested for B0 rehydration");
  unsigned int event_index = 0;
  std::string out_track_source;
  std::string out_cluster_source;
  int out_run = -1;
  int out_segment = -1;
  int out_sync_event = -1;
  int out_event_sequence = -1;
  unsigned long long out_stream_event_ordinal = 0;
  unsigned int object_count = 0;

  events.Branch("event_index", &event_index);
  events.Branch("track_source", &out_track_source);
  events.Branch("cluster_source", &out_cluster_source);
  events.Branch("run", &out_run);
  events.Branch("segment", &out_segment);
  events.Branch("sync_event", &out_sync_event);
  events.Branch("event_sequence", &out_event_sequence);
  events.Branch("stream_event_ordinal", &out_stream_event_ordinal);
  events.Branch("object_count", &object_count);

  unsigned int next_event_index = 0;
  std::map<CPMB0::EventKey, unsigned int> event_indices;
  for (const auto& [key, requests] : event_requests)
  {
    event_index = next_event_index++;
    event_indices.emplace(key, event_index);
    out_track_source = key.track_source;
    out_cluster_source = key.cluster_source;
    out_run = key.run;
    out_segment = key.segment;
    out_sync_event = key.sync_event;
    out_event_sequence = key.event_sequence;
    out_stream_event_ordinal = key.stream_event_ordinal;
    object_count = requests.size();
    events.Fill();
  }

  TTree objects("cpm_object_requests", "Track/state/cluster objects requested by CPM B0");
  unsigned int object_event_index = 0;
  unsigned int object_track_id = 0;
  unsigned long long object_cluskey = 0;
  unsigned long long object_hitsetkey = 0;
  unsigned int object_subsurfkey = 0;
  int object_iphi = -1;
  int object_ir = -1;
  int object_iz = -1;

  objects.Branch("event_index", &object_event_index);
  objects.Branch("track_id", &object_track_id);
  objects.Branch("cluskey", &object_cluskey);
  objects.Branch("hitsetkey", &object_hitsetkey);
  objects.Branch("subsurfkey", &object_subsurfkey);
  objects.Branch("iphi", &object_iphi);
  objects.Branch("ir", &object_ir);
  objects.Branch("iz", &object_iz);

  for (const auto& [key, requests] : event_requests)
  {
    object_event_index = event_indices[key];
    for (const auto& request : requests)
    {
      object_track_id = request.track_id;
      object_cluskey = request.cluskey;
      object_hitsetkey = request.hitsetkey;
      object_subsurfkey = request.subsurfkey;
      object_iphi = request.iphi;
      object_ir = request.ir;
      object_iz = request.iz;
      objects.Fill();
    }
  }

  events.Write();
  objects.Write();
  output.Close();

  std::cout << "CPM_QA_B0_BuildEventIndex - input records: " << records_read << std::endl;
  std::cout << "CPM_QA_B0_BuildEventIndex - unique events: " << event_requests.size() << std::endl;
  unsigned long long n_object_requests = 0;
  for (const auto& [key, requests] : event_requests)
  {
    (void) key;
    n_object_requests += requests.size();
  }
  std::cout << "CPM_QA_B0_BuildEventIndex - unique object requests: " << n_object_requests << std::endl;
  std::cout << "CPM_QA_B0_BuildEventIndex - output: " << output_file << std::endl;
}

void CPM_QA_B0_BuildEventIndexFromList(
    const std::string& input_list,
    const std::string& output_file = "CPM_QA_B0_event_index.root")
{
  std::ifstream input(input_list);
  if (!input.good())
  {
    std::cerr << "CPM_QA_B0_BuildEventIndexFromList - cannot open " << input_list << std::endl;
    return;
  }

  std::vector<std::string> input_files;
  std::string line;
  while (std::getline(input, line))
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    input_files.push_back(line);
  }

  std::cout << "CPM_QA_B0_BuildEventIndexFromList - input list: " << input_list << std::endl;
  std::cout << "CPM_QA_B0_BuildEventIndexFromList - files: " << input_files.size() << std::endl;
  CPM_QA_B0_BuildEventIndex(input_files, output_file);
}

void CPM_QA_B0_BuildEventIndex(
    const std::string& input_file_or_list,
    const std::string& output_file,
    const bool input_is_list)
{
  if (input_is_list)
  {
    CPM_QA_B0_BuildEventIndexFromList(input_file_or_list, output_file);
  }
  else
  {
    CPM_QA_B0_BuildEventIndex(std::vector<std::string>{input_file_or_list}, output_file);
  }
}

void CPM_QA_B0_BuildEventIndex(
    const std::string& input_file,
    const std::string& output_file = "CPM_QA_B0_event_index.root")
{
  CPM_QA_B0_BuildEventIndex(std::vector<std::string>{input_file}, output_file);
}
