/*
 * CPM Job B1 local line-line PoCA prototype.
 *
 * This macro reads Job A cpm_records, groups ACTS-ready state snapshots by
 * voxel, forms track-state pairs inside each voxel, and computes the first CPM
 * local line-line point of closest approach. It does not require seed objects
 * or TRKR_CLUSTER in the CPM mini-DST.
 */

#include <CPMLocalLinePoCA.h>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <vector>

R__LOAD_LIBRARY(libcpm.so)

namespace CPMB1
{
  struct VoxelKey
  {
    int iphi = -1;
    int ir = -1;
    int iz = -1;

    bool operator<(const VoxelKey& rhs) const
    {
      return std::tie(iphi, ir, iz) < std::tie(rhs.iphi, rhs.ir, rhs.iz);
    }
  };

  struct EventTrackKey
  {
    std::string cluster_source;
    std::string track_source;
    int run = -1;
    int segment = -1;
    int sync_event = -1;
    int event_sequence = -1;
    unsigned long long stream_event_ordinal = 0;
    unsigned int track_id = 0;

    bool operator==(const EventTrackKey& rhs) const
    {
      return std::tie(cluster_source, track_source, run, segment, sync_event,
                      event_sequence, stream_event_ordinal, track_id) ==
             std::tie(rhs.cluster_source, rhs.track_source, rhs.run, rhs.segment,
                      rhs.sync_event, rhs.event_sequence, rhs.stream_event_ordinal,
                      rhs.track_id);
    }
  };

  struct Record
  {
    Long64_t entry = -1;
    EventTrackKey event_track;
    unsigned long long cluskey = 0;
    VoxelKey voxel;
    cpm::Vector3 voxel_center;
    cpm::Vector3 offset;
    cpm::Vector3 state_position;
    cpm::Vector3 state_momentum;
  };

  double wrap_delta_phi(double value)
  {
    constexpr double pi = 3.141592653589793238462643383279502884;
    while (value > pi)
    {
      value -= 2.0 * pi;
    }
    while (value <= -pi)
    {
      value += 2.0 * pi;
    }
    return value;
  }

  std::vector<std::string> read_file_list(const std::string& input_list)
  {
    std::ifstream input(input_list);
    std::vector<std::string> files;
    std::string line;
    while (std::getline(input, line))
    {
      if (line.empty() || line[0] == '#')
      {
        continue;
      }
      files.push_back(line);
    }
    return files;
  }
}

void CPM_B1_LocalLinePoCA(
    const std::vector<std::string>& input_files,
    const std::string& output_file = "CPM_B1_local_line_poca.root",
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 500)
{
  TChain chain("cpm_records");
  for (const auto& file : input_files)
  {
    chain.Add(file.c_str());
  }

  std::string* cluster_source = nullptr;
  std::string* track_source = nullptr;
  int run = -1;
  int segment = -1;
  int sync_event = -1;
  int event_sequence = -1;
  unsigned long long stream_event_ordinal = 0;
  unsigned int track_id = 0;
  unsigned long long cluskey = 0;
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

  std::map<CPMB1::VoxelKey, std::vector<CPMB1::Record>> records_by_voxel;

  const auto entries = chain.GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    chain.GetEntry(entry);

    CPMB1::Record record;
    record.entry = entry;
    record.event_track.cluster_source = cluster_source ? *cluster_source : "";
    record.event_track.track_source = track_source ? *track_source : "";
    record.event_track.run = run;
    record.event_track.segment = segment;
    record.event_track.sync_event = sync_event;
    record.event_track.event_sequence = event_sequence;
    record.event_track.stream_event_ordinal = stream_event_ordinal;
    record.event_track.track_id = track_id;
    record.cluskey = cluskey;
    record.voxel = {iphi, ir, iz};
    record.voxel_center = {voxel_x, voxel_y, voxel_z};
    record.offset = {offset_x, offset_y, offset_z};
    record.state_position = {state_x, state_y, state_z};
    record.state_momentum = {state_px, state_py, state_pz};

    if (record.voxel.iphi < 0 || record.voxel.ir < 0 || record.voxel.iz < 0)
    {
      continue;
    }

    records_by_voxel[record.voxel].push_back(record);
  }

  TFile output(output_file.c_str(), "RECREATE");

  TTree pairs("cpm_poca_pairs", "CPM local line-line PoCA pairs");
  int out_iphi = -1;
  int out_ir = -1;
  int out_iz = -1;
  Long64_t entry_a = -1;
  Long64_t entry_b = -1;
  unsigned int track_id_a = 0;
  unsigned int track_id_b = 0;
  unsigned long long cluskey_a = 0;
  unsigned long long cluskey_b = 0;
  double dca = std::numeric_limits<double>::quiet_NaN();
  double s = std::numeric_limits<double>::quiet_NaN();
  double t = std::numeric_limits<double>::quiet_NaN();
  double point_ax = std::numeric_limits<double>::quiet_NaN();
  double point_ay = std::numeric_limits<double>::quiet_NaN();
  double point_az = std::numeric_limits<double>::quiet_NaN();
  double point_bx = std::numeric_limits<double>::quiet_NaN();
  double point_by = std::numeric_limits<double>::quiet_NaN();
  double point_bz = std::numeric_limits<double>::quiet_NaN();
  double midpoint_x = std::numeric_limits<double>::quiet_NaN();
  double midpoint_y = std::numeric_limits<double>::quiet_NaN();
  double midpoint_z = std::numeric_limits<double>::quiet_NaN();
  double voxel_center_x = std::numeric_limits<double>::quiet_NaN();
  double voxel_center_y = std::numeric_limits<double>::quiet_NaN();
  double voxel_center_z = std::numeric_limits<double>::quiet_NaN();
  double delta_x = std::numeric_limits<double>::quiet_NaN();
  double delta_y = std::numeric_limits<double>::quiet_NaN();
  double delta_z = std::numeric_limits<double>::quiet_NaN();
  double delta_r = std::numeric_limits<double>::quiet_NaN();
  double delta_phi = std::numeric_limits<double>::quiet_NaN();
  double delta_rphi = std::numeric_limits<double>::quiet_NaN();

  pairs.Branch("iphi", &out_iphi);
  pairs.Branch("ir", &out_ir);
  pairs.Branch("iz", &out_iz);
  pairs.Branch("entry_a", &entry_a);
  pairs.Branch("entry_b", &entry_b);
  pairs.Branch("track_id_a", &track_id_a);
  pairs.Branch("track_id_b", &track_id_b);
  pairs.Branch("cluskey_a", &cluskey_a);
  pairs.Branch("cluskey_b", &cluskey_b);
  pairs.Branch("dca", &dca);
  pairs.Branch("s", &s);
  pairs.Branch("t", &t);
  pairs.Branch("point_ax", &point_ax);
  pairs.Branch("point_ay", &point_ay);
  pairs.Branch("point_az", &point_az);
  pairs.Branch("point_bx", &point_bx);
  pairs.Branch("point_by", &point_by);
  pairs.Branch("point_bz", &point_bz);
  pairs.Branch("midpoint_x", &midpoint_x);
  pairs.Branch("midpoint_y", &midpoint_y);
  pairs.Branch("midpoint_z", &midpoint_z);
  pairs.Branch("voxel_center_x", &voxel_center_x);
  pairs.Branch("voxel_center_y", &voxel_center_y);
  pairs.Branch("voxel_center_z", &voxel_center_z);
  pairs.Branch("delta_x", &delta_x);
  pairs.Branch("delta_y", &delta_y);
  pairs.Branch("delta_z", &delta_z);
  pairs.Branch("delta_r", &delta_r);
  pairs.Branch("delta_phi", &delta_phi);
  pairs.Branch("delta_rphi", &delta_rphi);

  cpm::LocalLinePoCAOptions options;
  options.min_sin_angle = min_sin_angle;

  unsigned long long candidate_pairs = 0;
  unsigned long long accepted_pairs = 0;
  unsigned int processed_voxels = 0;
  unsigned int skipped_large_voxels = 0;

  for (const auto& [voxel, records] : records_by_voxel)
  {
    if (records.size() < 2)
    {
      continue;
    }
    if (max_records_per_voxel > 0 && records.size() > max_records_per_voxel)
    {
      ++skipped_large_voxels;
      continue;
    }

    ++processed_voxels;
    for (std::size_t ia = 0; ia < records.size(); ++ia)
    {
      for (std::size_t ib = ia + 1; ib < records.size(); ++ib)
      {
        const auto& a = records[ia];
        const auto& b = records[ib];
        if (a.event_track == b.event_track)
        {
          continue;
        }

        ++candidate_pairs;
        const cpm::Vector3 point_a = a.state_position - a.offset;
        const cpm::Vector3 point_b = b.state_position - b.offset;
        const auto result = cpm::computeLocalLinePoCA(
            point_a,
            a.state_momentum,
            point_b,
            b.state_momentum,
            options);

        if (!result.valid || !(result.dca <= max_pair_dca))
        {
          continue;
        }

        out_iphi = voxel.iphi;
        out_ir = voxel.ir;
        out_iz = voxel.iz;
        entry_a = a.entry;
        entry_b = b.entry;
        track_id_a = a.event_track.track_id;
        track_id_b = b.event_track.track_id;
        cluskey_a = a.cluskey;
        cluskey_b = b.cluskey;
        dca = result.dca;
        s = result.s;
        t = result.t;
        point_ax = result.point_a.x;
        point_ay = result.point_a.y;
        point_az = result.point_a.z;
        point_bx = result.point_b.x;
        point_by = result.point_b.y;
        point_bz = result.point_b.z;
        midpoint_x = result.midpoint.x;
        midpoint_y = result.midpoint.y;
        midpoint_z = result.midpoint.z;
        voxel_center_x = a.voxel_center.x;
        voxel_center_y = a.voxel_center.y;
        voxel_center_z = a.voxel_center.z;
        delta_x = voxel_center_x - midpoint_x;
        delta_y = voxel_center_y - midpoint_y;
        delta_z = voxel_center_z - midpoint_z;

        const double voxel_r = std::hypot(voxel_center_x, voxel_center_y);
        const double midpoint_r = std::hypot(midpoint_x, midpoint_y);
        const double voxel_phi = std::atan2(voxel_center_y, voxel_center_x);
        const double midpoint_phi = std::atan2(midpoint_y, midpoint_x);
        delta_r = voxel_r - midpoint_r;
        delta_phi = CPMB1::wrap_delta_phi(voxel_phi - midpoint_phi);
        delta_rphi = voxel_r * delta_phi;

        pairs.Fill();
        ++accepted_pairs;
      }
    }
  }

  TTree summary("cpm_b1_summary", "CPM B1 local line-line PoCA summary");
  unsigned long long input_records = entries;
  unsigned int input_files_count = input_files.size();
  unsigned int voxel_count = records_by_voxel.size();
  double summary_max_pair_dca = max_pair_dca;
  double summary_min_sin_angle = min_sin_angle;
  unsigned int summary_max_records_per_voxel = max_records_per_voxel;

  summary.Branch("input_records", &input_records);
  summary.Branch("input_files", &input_files_count);
  summary.Branch("voxel_count", &voxel_count);
  summary.Branch("processed_voxels", &processed_voxels);
  summary.Branch("skipped_large_voxels", &skipped_large_voxels);
  summary.Branch("candidate_pairs", &candidate_pairs);
  summary.Branch("accepted_pairs", &accepted_pairs);
  summary.Branch("max_pair_dca", &summary_max_pair_dca);
  summary.Branch("min_sin_angle", &summary_min_sin_angle);
  summary.Branch("max_records_per_voxel", &summary_max_records_per_voxel);
  summary.Fill();

  pairs.Write();
  summary.Write();
  output.Close();

  std::cout << "CPM_B1_LocalLinePoCA - input records: " << input_records << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - voxels: " << voxel_count << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - processed voxels: " << processed_voxels << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - skipped large voxels: " << skipped_large_voxels << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - candidate pairs: " << candidate_pairs << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - accepted pairs: " << accepted_pairs << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - output: " << output_file << std::endl;
}

void CPM_B1_LocalLinePoCA(
    const std::string& input_file,
    const std::string& output_file = "CPM_B1_local_line_poca.root",
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 500)
{
  CPM_B1_LocalLinePoCA(
      std::vector<std::string>{input_file},
      output_file,
      max_pair_dca,
      min_sin_angle,
      max_records_per_voxel);
}

void CPM_B1_LocalLinePoCA(
    const std::string& input_file_or_list,
    const std::string& output_file,
    const bool input_is_list,
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 500)
{
  const auto input_files = input_is_list ?
      CPMB1::read_file_list(input_file_or_list) :
      std::vector<std::string>{input_file_or_list};

  CPM_B1_LocalLinePoCA(
      input_files,
      output_file,
      max_pair_dca,
      min_sin_angle,
      max_records_per_voxel);
}
