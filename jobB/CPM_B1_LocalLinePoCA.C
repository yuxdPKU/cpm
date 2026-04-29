/*
 * CPM Job B1 local line-line PoCA prototype.
 *
 * This macro reads Job A cpm_records, groups ACTS-ready state snapshots by
 * voxel, forms same-charge track-state pairs inside each voxel, and computes
 * the first CPM local line-line point of closest approach. It does not require
 * seed objects or TRKR_CLUSTER in the CPM mini-DST.
 */

#include <CPMLocalLinePoCA.h>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
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
    int charge = 0;
    double pt = std::numeric_limits<double>::quiet_NaN();
    VoxelKey voxel;
    cpm::Vector3 voxel_center;
    cpm::Vector3 offset;
    cpm::Vector3 state_position;
    cpm::Vector3 state_momentum;
  };

  struct GridMetadata
  {
    int phi_bins = -1;
    int r_bins = -1;
    int z_bins = -1;
    bool valid = false;
  };

  struct VoxelSummary
  {
    unsigned long long unique_tracks = 0;
    unsigned long long unique_track_pairs = 0;
    unsigned long long selected_records = 0;
    unsigned long long selected_unique_tracks = 0;
    unsigned long long selected_record_pairs = 0;
    unsigned long long pt_rejected_records = 0;
    unsigned long long cap_dropped_records = 0;
    unsigned long long raw_pairs = 0;
    unsigned long long same_charge_pairs = 0;
    unsigned long long candidate_pairs = 0;
    unsigned long long accepted_pairs = 0;
  };

  using UniqueTrackId = std::tuple<
      std::string,
      std::string,
      int,
      int,
      int,
      int,
      unsigned long long,
      unsigned int>;

  unsigned long long pair_count(const std::size_t entries)
  {
    return entries < 2 ? 0 : static_cast<unsigned long long>(entries) * (entries - 1) / 2;
  }

  UniqueTrackId make_unique_track_id(const Record& record)
  {
    return std::make_tuple(
        record.event_track.cluster_source,
        record.event_track.track_source,
        record.event_track.run,
        record.event_track.segment,
        record.event_track.sync_event,
        record.event_track.event_sequence,
        record.event_track.stream_event_ordinal,
        record.event_track.track_id);
  }

  unsigned long long count_unique_tracks(const std::vector<Record>& records)
  {
    std::set<UniqueTrackId> unique_tracks;
    for (const auto& record : records)
    {
      unique_tracks.insert(make_unique_track_id(record));
    }
    return unique_tracks.size();
  }

  unsigned long long count_unique_tracks(const std::vector<const Record*>& records)
  {
    std::set<UniqueTrackId> unique_tracks;
    for (const auto* record : records)
    {
      unique_tracks.insert(make_unique_track_id(*record));
    }
    return unique_tracks.size();
  }

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

  bool has_good_curvature_proxy(const Record& record, const double min_pair_pt)
  {
    return record.charge != 0 &&
           std::isfinite(record.pt) &&
           record.pt > 0.0 &&
           record.pt >= min_pair_pt;
  }

  double pair_weight_from_curvature_proxy(const Record& a, const Record& b)
  {
    return 1.0 / (a.pt * b.pt);
  }

  GridMetadata load_grid_metadata(const std::vector<std::string>& input_files)
  {
    GridMetadata metadata;
    if (input_files.empty())
    {
      return metadata;
    }

    TFile input(input_files.front().c_str(), "READ");
    if (input.IsZombie())
    {
      return metadata;
    }

    auto* tree = dynamic_cast<TTree*>(input.Get("cpm_metadata"));
    if (!tree || tree->GetEntries() <= 0)
    {
      return metadata;
    }

    tree->SetBranchAddress("phi_bins", &metadata.phi_bins);
    tree->SetBranchAddress("r_bins", &metadata.r_bins);
    tree->SetBranchAddress("z_bins", &metadata.z_bins);
    tree->GetEntry(0);
    metadata.valid = metadata.phi_bins > 0 && metadata.r_bins > 0 && metadata.z_bins > 0;
    return metadata;
  }

  void print_voxel_summary(
      const VoxelKey& voxel,
      const GridMetadata& metadata,
      const std::size_t records,
      const unsigned int positive_records,
      const unsigned int negative_records,
      const VoxelSummary& summary,
      const std::string& status)
  {
    std::cout << "CPM_B1_LocalLinePoCA - voxel (iphi,ir,iz)=("
              << voxel.iphi << "," << voxel.ir << "," << voxel.iz << ")";
    if (metadata.valid)
    {
      std::cout << " / bins=(" << metadata.phi_bins << ","
                << metadata.r_bins << "," << metadata.z_bins << ")";
    }
    else
    {
      std::cout << " / bins=(unknown)";
    }
    std::cout << " records=" << records
              << " unique_tracks=" << summary.unique_tracks
              << " unique_track_pairs=" << summary.unique_track_pairs
              << " selected_records=" << summary.selected_records
              << " selected_unique_tracks=" << summary.selected_unique_tracks
              << " selected_charge_records(+,-)=("
              << positive_records << "," << negative_records << ")"
              << " pt_rejected_records=" << summary.pt_rejected_records
              << " cap_dropped_records=" << summary.cap_dropped_records
              << " raw_record_pairs=" << summary.raw_pairs
              << " selected_record_pairs=" << summary.selected_record_pairs
              << " same_charge_pairs=" << summary.same_charge_pairs
              << " candidate_pairs=" << summary.candidate_pairs
              << " accepted_pairs=" << summary.accepted_pairs
              << " status=" << status
              << std::endl;
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
    const unsigned int max_records_per_voxel = 500,
    const unsigned int min_records_per_charge = 2,
    const bool print_voxel_summaries = true,
    const double min_pair_pt = 0.0,
    const unsigned int max_pair_records_per_voxel = 0)
{
  TChain chain("cpm_records");
  for (const auto& file : input_files)
  {
    chain.Add(file.c_str());
  }
  const auto grid_metadata = CPMB1::load_grid_metadata(input_files);

  std::string* cluster_source = nullptr;
  std::string* track_source = nullptr;
  int run = -1;
  int segment = -1;
  int sync_event = -1;
  int event_sequence = -1;
  unsigned long long stream_event_ordinal = 0;
  unsigned int track_id = 0;
  unsigned long long cluskey = 0;
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
    record.charge = charge;
    record.pt = pt;
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
  int charge_a = 0;
  int charge_b = 0;
  double pt_a = std::numeric_limits<double>::quiet_NaN();
  double pt_b = std::numeric_limits<double>::quiet_NaN();
  double pair_weight = std::numeric_limits<double>::quiet_NaN();
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
  pairs.Branch("charge_a", &charge_a);
  pairs.Branch("charge_b", &charge_b);
  pairs.Branch("pt_a", &pt_a);
  pairs.Branch("pt_b", &pt_b);
  pairs.Branch("pair_weight", &pair_weight);
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

  TTree voxel_summaries("cpm_b1_voxel_summary", "CPM B1 per-voxel pair input summary");
  int summary_iphi = -1;
  int summary_ir = -1;
  int summary_iz = -1;
  unsigned long long summary_records = 0;
  unsigned long long summary_unique_tracks = 0;
  unsigned long long summary_unique_track_pairs = 0;
  unsigned long long summary_selected_records = 0;
  unsigned long long summary_selected_unique_tracks = 0;
  unsigned long long summary_pt_rejected_records = 0;
  unsigned long long summary_cap_dropped_records = 0;
  unsigned long long summary_raw_record_pairs = 0;
  unsigned long long summary_selected_record_pairs = 0;
  unsigned long long summary_same_charge_pairs = 0;
  unsigned long long summary_candidate_pairs = 0;
  unsigned long long summary_accepted_pairs = 0;
  unsigned int summary_positive_records = 0;
  unsigned int summary_negative_records = 0;
  std::string summary_status;

  voxel_summaries.Branch("iphi", &summary_iphi);
  voxel_summaries.Branch("ir", &summary_ir);
  voxel_summaries.Branch("iz", &summary_iz);
  voxel_summaries.Branch("records", &summary_records);
  voxel_summaries.Branch("unique_tracks", &summary_unique_tracks);
  voxel_summaries.Branch("unique_track_pairs", &summary_unique_track_pairs);
  voxel_summaries.Branch("selected_records", &summary_selected_records);
  voxel_summaries.Branch("selected_unique_tracks", &summary_selected_unique_tracks);
  voxel_summaries.Branch("pt_rejected_records", &summary_pt_rejected_records);
  voxel_summaries.Branch("cap_dropped_records", &summary_cap_dropped_records);
  voxel_summaries.Branch("raw_record_pairs", &summary_raw_record_pairs);
  voxel_summaries.Branch("selected_record_pairs", &summary_selected_record_pairs);
  voxel_summaries.Branch("same_charge_pairs", &summary_same_charge_pairs);
  voxel_summaries.Branch("candidate_pairs", &summary_candidate_pairs);
  voxel_summaries.Branch("accepted_pairs", &summary_accepted_pairs);
  voxel_summaries.Branch("positive_records", &summary_positive_records);
  voxel_summaries.Branch("negative_records", &summary_negative_records);
  voxel_summaries.Branch("status", &summary_status);

  auto fill_voxel_summary = [&](
      const CPMB1::VoxelKey& voxel,
      const std::size_t records,
      const unsigned int positive_records,
      const unsigned int negative_records,
      const CPMB1::VoxelSummary& voxel_summary,
      const std::string& status)
  {
    summary_iphi = voxel.iphi;
    summary_ir = voxel.ir;
    summary_iz = voxel.iz;
    summary_records = records;
    summary_unique_tracks = voxel_summary.unique_tracks;
    summary_unique_track_pairs = voxel_summary.unique_track_pairs;
    summary_selected_records = voxel_summary.selected_records;
    summary_selected_unique_tracks = voxel_summary.selected_unique_tracks;
    summary_pt_rejected_records = voxel_summary.pt_rejected_records;
    summary_cap_dropped_records = voxel_summary.cap_dropped_records;
    summary_raw_record_pairs = voxel_summary.raw_pairs;
    summary_selected_record_pairs = voxel_summary.selected_record_pairs;
    summary_same_charge_pairs = voxel_summary.same_charge_pairs;
    summary_candidate_pairs = voxel_summary.candidate_pairs;
    summary_accepted_pairs = voxel_summary.accepted_pairs;
    summary_positive_records = positive_records;
    summary_negative_records = negative_records;
    summary_status = status;
    voxel_summaries.Fill();
  };

  cpm::LocalLinePoCAOptions options;
  options.min_sin_angle = min_sin_angle;

  unsigned long long candidate_pairs = 0;
  unsigned long long accepted_pairs = 0;
  unsigned int processed_voxels = 0;
  unsigned int skipped_low_charge_voxels = 0;
  unsigned int skipped_large_voxels = 0;

  for (const auto& [voxel, records] : records_by_voxel)
  {
    CPMB1::VoxelSummary voxel_summary;
    std::vector<const CPMB1::Record*> selected_records;
    selected_records.reserve(records.size());

    for (const auto& record : records)
    {
      if (!CPMB1::has_good_curvature_proxy(record, min_pair_pt))
      {
        ++voxel_summary.pt_rejected_records;
        continue;
      }
      selected_records.push_back(&record);
    }

    std::stable_sort(
        selected_records.begin(),
        selected_records.end(),
        [](const CPMB1::Record* lhs, const CPMB1::Record* rhs)
        {
          if (lhs->pt != rhs->pt)
          {
            return lhs->pt > rhs->pt;
          }
          if (lhs->event_track.track_id != rhs->event_track.track_id)
          {
            return lhs->event_track.track_id < rhs->event_track.track_id;
          }
          return lhs->cluskey < rhs->cluskey;
        });

    if (max_pair_records_per_voxel > 0 &&
        selected_records.size() > max_pair_records_per_voxel)
    {
      voxel_summary.cap_dropped_records =
          selected_records.size() - max_pair_records_per_voxel;
      selected_records.resize(max_pair_records_per_voxel);
    }

    unsigned int positive_records = 0;
    unsigned int negative_records = 0;
    for (const auto* record : selected_records)
    {
      if (record->charge > 0)
      {
        ++positive_records;
      }
      else if (record->charge < 0)
      {
        ++negative_records;
      }
    }
    voxel_summary.unique_tracks = CPMB1::count_unique_tracks(records);
    voxel_summary.unique_track_pairs =
        CPMB1::pair_count(static_cast<std::size_t>(voxel_summary.unique_tracks));
    voxel_summary.selected_records = selected_records.size();
    voxel_summary.selected_unique_tracks = CPMB1::count_unique_tracks(selected_records);
    voxel_summary.raw_pairs = CPMB1::pair_count(records.size());
    voxel_summary.selected_record_pairs = CPMB1::pair_count(selected_records.size());
    voxel_summary.same_charge_pairs =
        CPMB1::pair_count(positive_records) + CPMB1::pair_count(negative_records);

    if (max_records_per_voxel > 0 && records.size() > max_records_per_voxel)
    {
      ++skipped_large_voxels;
      if (print_voxel_summaries)
      {
        CPMB1::print_voxel_summary(
            voxel,
            grid_metadata,
            records.size(),
            positive_records,
            negative_records,
            voxel_summary,
            "skipped_too_many_records");
      }
      fill_voxel_summary(
          voxel,
          records.size(),
          positive_records,
          negative_records,
          voxel_summary,
          "skipped_too_many_records");
      continue;
    }
    if (selected_records.size() < 2)
    {
      if (print_voxel_summaries)
      {
        CPMB1::print_voxel_summary(
            voxel,
            grid_metadata,
            records.size(),
            positive_records,
            negative_records,
            voxel_summary,
            "skipped_fewer_than_two_selected_records");
      }
      fill_voxel_summary(
          voxel,
          records.size(),
          positive_records,
          negative_records,
          voxel_summary,
          "skipped_fewer_than_two_selected_records");
      continue;
    }

    if (positive_records < min_records_per_charge &&
        negative_records < min_records_per_charge)
    {
      ++skipped_low_charge_voxels;
      if (print_voxel_summaries)
      {
        CPMB1::print_voxel_summary(
            voxel,
            grid_metadata,
            records.size(),
            positive_records,
            negative_records,
            voxel_summary,
            "skipped_low_charge_records");
      }
      fill_voxel_summary(
          voxel,
          records.size(),
          positive_records,
          negative_records,
          voxel_summary,
          "skipped_low_charge_records");
      continue;
    }

    ++processed_voxels;
    for (std::size_t ia = 0; ia < selected_records.size(); ++ia)
    {
      for (std::size_t ib = ia + 1; ib < selected_records.size(); ++ib)
      {
        const auto& a = *selected_records[ia];
        const auto& b = *selected_records[ib];
        if (a.event_track == b.event_track)
        {
          continue;
        }
        if (!CPMB1::has_good_curvature_proxy(a, min_pair_pt) ||
            !CPMB1::has_good_curvature_proxy(b, min_pair_pt) ||
            a.charge != b.charge)
        {
          continue;
        }
        if ((a.charge > 0 && positive_records < min_records_per_charge) ||
            (a.charge < 0 && negative_records < min_records_per_charge))
        {
          continue;
        }

        ++candidate_pairs;
        ++voxel_summary.candidate_pairs;
        pair_weight = CPMB1::pair_weight_from_curvature_proxy(a, b);
        if (!std::isfinite(pair_weight) || pair_weight <= 0.0)
        {
          continue;
        }

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
        charge_a = a.charge;
        charge_b = b.charge;
        pt_a = a.pt;
        pt_b = b.pt;
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
        ++voxel_summary.accepted_pairs;
      }
    }

    if (print_voxel_summaries)
    {
      CPMB1::print_voxel_summary(
          voxel,
          grid_metadata,
          records.size(),
          positive_records,
          negative_records,
          voxel_summary,
          "processed");
    }
    fill_voxel_summary(
        voxel,
        records.size(),
        positive_records,
        negative_records,
        voxel_summary,
        "processed");
  }

  TTree summary("cpm_b1_summary", "CPM B1 local line-line PoCA summary");
  unsigned long long input_records = entries;
  unsigned int input_files_count = input_files.size();
  unsigned int voxel_count = records_by_voxel.size();
  double summary_max_pair_dca = max_pair_dca;
  double summary_min_sin_angle = min_sin_angle;
  unsigned int summary_max_records_per_voxel = max_records_per_voxel;
  unsigned int summary_min_records_per_charge = min_records_per_charge;
  bool summary_print_voxel_summaries = print_voxel_summaries;
  double summary_min_pair_pt = min_pair_pt;
  unsigned int summary_max_pair_records_per_voxel = max_pair_records_per_voxel;
  int summary_phi_bins = grid_metadata.phi_bins;
  int summary_r_bins = grid_metadata.r_bins;
  int summary_z_bins = grid_metadata.z_bins;
  bool summary_grid_metadata_valid = grid_metadata.valid;

  summary.Branch("input_records", &input_records);
  summary.Branch("input_files", &input_files_count);
  summary.Branch("voxel_count", &voxel_count);
  summary.Branch("processed_voxels", &processed_voxels);
  summary.Branch("skipped_large_voxels", &skipped_large_voxels);
  summary.Branch("skipped_low_charge_voxels", &skipped_low_charge_voxels);
  summary.Branch("candidate_pairs", &candidate_pairs);
  summary.Branch("accepted_pairs", &accepted_pairs);
  summary.Branch("max_pair_dca", &summary_max_pair_dca);
  summary.Branch("min_sin_angle", &summary_min_sin_angle);
  summary.Branch("max_records_per_voxel", &summary_max_records_per_voxel);
  summary.Branch("min_records_per_charge", &summary_min_records_per_charge);
  summary.Branch("print_voxel_summaries", &summary_print_voxel_summaries);
  summary.Branch("min_pair_pt", &summary_min_pair_pt);
  summary.Branch("max_pair_records_per_voxel", &summary_max_pair_records_per_voxel);
  summary.Branch("phi_bins", &summary_phi_bins);
  summary.Branch("r_bins", &summary_r_bins);
  summary.Branch("z_bins", &summary_z_bins);
  summary.Branch("grid_metadata_valid", &summary_grid_metadata_valid);
  summary.Fill();

  pairs.Write();
  voxel_summaries.Write();
  summary.Write();
  output.Close();

  std::cout << "CPM_B1_LocalLinePoCA - input records: " << input_records << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - voxels: " << voxel_count << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - processed voxels: " << processed_voxels << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - skipped large voxels: " << skipped_large_voxels << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - skipped low-charge voxels: " << skipped_low_charge_voxels << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - candidate pairs: " << candidate_pairs << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - accepted pairs: " << accepted_pairs << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - min pair pt: " << min_pair_pt << std::endl;
  std::cout << "CPM_B1_LocalLinePoCA - max pair records per voxel: " << max_pair_records_per_voxel << std::endl;
  if (grid_metadata.valid)
  {
    std::cout << "CPM_B1_LocalLinePoCA - grid bins: ("
              << grid_metadata.phi_bins << ", "
              << grid_metadata.r_bins << ", "
              << grid_metadata.z_bins << ")" << std::endl;
  }
  std::cout << "CPM_B1_LocalLinePoCA - output: " << output_file << std::endl;
}

void CPM_B1_LocalLinePoCA(
    const std::string& input_file,
    const std::string& output_file = "CPM_B1_local_line_poca.root",
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 500,
    const unsigned int min_records_per_charge = 2,
    const bool print_voxel_summaries = true,
    const double min_pair_pt = 0.0,
    const unsigned int max_pair_records_per_voxel = 0)
{
  CPM_B1_LocalLinePoCA(
      std::vector<std::string>{input_file},
      output_file,
      max_pair_dca,
      min_sin_angle,
      max_records_per_voxel,
      min_records_per_charge,
      print_voxel_summaries,
      min_pair_pt,
      max_pair_records_per_voxel);
}

void CPM_B1_LocalLinePoCA(
    const std::string& input_file_or_list,
    const std::string& output_file,
    const bool input_is_list,
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 500,
    const unsigned int min_records_per_charge = 2,
    const bool print_voxel_summaries = true,
    const double min_pair_pt = 0.0,
    const unsigned int max_pair_records_per_voxel = 0)
{
  const auto input_files = input_is_list ?
      CPMB1::read_file_list(input_file_or_list) :
      std::vector<std::string>{input_file_or_list};

  CPM_B1_LocalLinePoCA(
      input_files,
      output_file,
      max_pair_dca,
      min_sin_angle,
      max_records_per_voxel,
      min_records_per_charge,
      print_voxel_summaries,
      min_pair_pt,
      max_pair_records_per_voxel);
}
