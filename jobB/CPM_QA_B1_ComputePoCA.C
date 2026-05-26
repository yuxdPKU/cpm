/*
 * CPM QA Job B1 crossing-point PoCA prototype.
 *
 * This macro reads Job A cpm_records, groups ACTS-ready state snapshots by
 * voxel, forms opposite-charge track-state pairs inside each voxel, and computes
 * CPM crossing-point estimates. The default solver is the ideal-helix PoCA;
 * the v1 local line-line PoCA solver can be enabled for comparison.
 * It does not require seed objects or TRKR_CLUSTER in the CPM mini-DST.
 */

#include <CPMRecordReader.h>

#include <CPMCorrectionAccumulator.h>
#include <CPMPairUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
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
    Vector3 voxel_center;
    Vector3 offset;
    Vector3 state_position;
    Vector3 state_momentum;
  };

  struct GridMetadata
  {
    int phi_bins = -1;
    int r_bins = -1;
    int z_bins = -1;
    double r_min = std::numeric_limits<double>::quiet_NaN();
    double r_max = std::numeric_limits<double>::quiet_NaN();
    double z_min = std::numeric_limits<double>::quiet_NaN();
    double z_max = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
    bool consistent = true;
    std::string source_file;
    std::string message;
  };

  struct VoxelSummary
  {
    unsigned long long unique_tracks = 0;
    unsigned long long unique_track_pairs = 0;
    unsigned long long selected_records = 0;
    unsigned long long selected_unique_tracks = 0;
    unsigned long long selected_record_pairs = 0;
    unsigned long long pt_rejected_records = 0;
    unsigned long long duplicate_dropped_records = 0;
    unsigned long long cap_dropped_records = 0;
    unsigned long long batch_count = 0;
    unsigned long long batched_positive_records = 0;
    unsigned long long batched_negative_records = 0;
    unsigned long long batched_record_pairs = 0;
    unsigned long long batched_opposite_charge_pairs = 0;
    unsigned long long raw_pairs = 0;
    unsigned long long same_charge_pairs = 0;
    unsigned long long opposite_charge_pairs = 0;
    unsigned long long candidate_pairs = 0;
    unsigned long long accepted_pairs = 0;
  };

  struct BatchAccumulator : public CPMCorrectionAccumulator
  {
    unsigned long long same_event_track_pairs = 0;
    unsigned long long pt_rejected_pairs = 0;
    unsigned long long candidate_pairs = 0;
    unsigned long long invalid_weight_pairs = 0;
    unsigned long long invalid_poca_pairs = 0;
    unsigned long long dca_rejected_pairs = 0;
    unsigned long long accepted_pairs = 0;

    void add(
        const double accepted_delta_r,
        const double accepted_delta_rphi,
        const double accepted_delta_phi,
        const double accepted_delta_z,
        const double accepted_dca,
        const double accepted_weight,
        const double voxel_x,
        const double voxel_y,
        const double voxel_z)
    {
      ++accepted_pairs;
      CPMCorrectionAccumulator::add(
          accepted_delta_r,
          accepted_delta_rphi,
          accepted_delta_phi,
          accepted_delta_z,
          accepted_dca,
          accepted_weight,
          voxel_x,
          voxel_y,
          voxel_z);
    }
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

  void hash_mix_byte(std::uint64_t& hash, const unsigned char value)
  {
    constexpr std::uint64_t fnv_prime = 1099511628211ULL;
    hash ^= value;
    hash *= fnv_prime;
  }

  void hash_mix_uint64(std::uint64_t& hash, std::uint64_t value)
  {
    for (int i = 0; i < 8; ++i)
    {
      hash_mix_byte(hash, static_cast<unsigned char>(value & 0xffU));
      value >>= 8U;
    }
  }

  void hash_mix_string(std::uint64_t& hash, const std::string& value)
  {
    for (const auto character : value)
    {
      hash_mix_byte(hash, static_cast<unsigned char>(character));
    }
    hash_mix_byte(hash, 0xffU);
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

  std::uint64_t record_selection_hash(const Record& record)
  {
    std::uint64_t hash = 14695981039346656037ULL;
    hash_mix_string(hash, record.event_track.cluster_source);
    hash_mix_string(hash, record.event_track.track_source);
    hash_mix_uint64(hash, static_cast<std::uint64_t>(record.event_track.run));
    hash_mix_uint64(hash, static_cast<std::uint64_t>(record.event_track.segment));
    hash_mix_uint64(hash, static_cast<std::uint64_t>(record.event_track.sync_event));
    hash_mix_uint64(hash, static_cast<std::uint64_t>(record.event_track.event_sequence));
    hash_mix_uint64(hash, record.event_track.stream_event_ordinal);
    hash_mix_uint64(hash, record.event_track.track_id);
    hash_mix_uint64(hash, record.cluskey);
    return hash;
  }

  double offset_magnitude2(const Record& record)
  {
    return record.offset.x * record.offset.x +
           record.offset.y * record.offset.y +
           record.offset.z * record.offset.z;
  }

  bool is_closer_to_voxel_center(const Record& candidate, const Record& current_best)
  {
    const double candidate_distance2 = offset_magnitude2(candidate);
    const double current_distance2 = offset_magnitude2(current_best);
    if (std::isfinite(candidate_distance2) && std::isfinite(current_distance2) &&
        candidate_distance2 != current_distance2)
    {
      return candidate_distance2 < current_distance2;
    }
    if (std::isfinite(candidate_distance2) != std::isfinite(current_distance2))
    {
      return std::isfinite(candidate_distance2);
    }
    if (candidate.cluskey != current_best.cluskey)
    {
      return candidate.cluskey < current_best.cluskey;
    }
    return candidate.entry < current_best.entry;
  }

  bool same_metadata_value(const double lhs, const double rhs)
  {
    return std::fabs(lhs - rhs) < 1.0e-9;
  }

  bool same_grid_metadata(const GridMetadata& lhs, const GridMetadata& rhs)
  {
    return lhs.phi_bins == rhs.phi_bins &&
           lhs.r_bins == rhs.r_bins &&
           lhs.z_bins == rhs.z_bins &&
           same_metadata_value(lhs.r_min, rhs.r_min) &&
           same_metadata_value(lhs.r_max, rhs.r_max) &&
           same_metadata_value(lhs.z_min, rhs.z_min) &&
           same_metadata_value(lhs.z_max, rhs.z_max);
  }

  GridMetadata load_grid_metadata_file(const std::string& input_file)
  {
    GridMetadata metadata;
    metadata.source_file = input_file;

    TFile input(input_file.c_str(), "READ");
    if (input.IsZombie())
    {
      metadata.message = "could not open " + input_file;
      return metadata;
    }

    auto* tree = dynamic_cast<TTree*>(input.Get("cpm_metadata"));
    if (!tree || tree->GetEntries() <= 0)
    {
      metadata.message = "missing cpm_metadata tree in " + input_file;
      return metadata;
    }

    for (const auto& branch_name : {"phi_bins", "r_bins", "z_bins", "r_min", "r_max", "z_min", "z_max"})
    {
      if (!tree->GetBranch(branch_name))
      {
        metadata.message = std::string("missing cpm_metadata branch ") + branch_name + " in " + input_file;
        return metadata;
      }
    }

    tree->SetBranchAddress("phi_bins", &metadata.phi_bins);
    tree->SetBranchAddress("r_bins", &metadata.r_bins);
    tree->SetBranchAddress("z_bins", &metadata.z_bins);
    tree->SetBranchAddress("r_min", &metadata.r_min);
    tree->SetBranchAddress("r_max", &metadata.r_max);
    tree->SetBranchAddress("z_min", &metadata.z_min);
    tree->SetBranchAddress("z_max", &metadata.z_max);
    tree->GetEntry(0);
    metadata.valid =
        metadata.phi_bins > 0 &&
        metadata.r_bins > 0 &&
        metadata.z_bins > 0 &&
        std::isfinite(metadata.r_min) &&
        std::isfinite(metadata.r_max) &&
        std::isfinite(metadata.z_min) &&
        std::isfinite(metadata.z_max) &&
        metadata.r_min < metadata.r_max &&
        metadata.z_min < metadata.z_max;
    if (!metadata.valid)
    {
      metadata.message = "invalid cpm_metadata values in " + input_file;
    }
    return metadata;
  }

  GridMetadata load_grid_metadata(const std::vector<std::string>& input_files)
  {
    GridMetadata metadata;
    if (input_files.empty())
    {
      metadata.consistent = false;
      metadata.message = "no input files";
      return metadata;
    }

    for (const auto& input_file : input_files)
    {
      const auto current = load_grid_metadata_file(input_file);
      if (!current.valid)
      {
        metadata = current;
        metadata.consistent = false;
        return metadata;
      }

      if (!metadata.valid)
      {
        metadata = current;
        continue;
      }

      if (!same_grid_metadata(metadata, current))
      {
        metadata.consistent = false;
        metadata.message = "cpm_metadata mismatch between " + metadata.source_file + " and " + input_file;
        return metadata;
      }
    }

    metadata.consistent = true;
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
    std::cout << "CPM_QA_B1_ComputePoCA - voxel (iphi,ir,iz)=("
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
              << " duplicate_dropped_records=" << summary.duplicate_dropped_records
              << " cap_dropped_records=" << summary.cap_dropped_records
              << " batches=" << summary.batch_count
              << " batched_charge_records(+,-)=("
              << summary.batched_positive_records << ","
              << summary.batched_negative_records << ")"
              << " raw_record_pairs=" << summary.raw_pairs
              << " selected_record_pairs=" << summary.selected_record_pairs
              << " batched_record_pairs=" << summary.batched_record_pairs
              << " same_charge_pairs=" << summary.same_charge_pairs
              << " opposite_charge_pairs=" << summary.opposite_charge_pairs
              << " batched_opposite_charge_pairs=" << summary.batched_opposite_charge_pairs
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

void CPM_QA_B1_ComputePoCA(
    const std::vector<std::string>& input_files,
    const std::string& output_file = "CPM_QA_B1_poca.root",
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 0,
    const unsigned int min_records_per_charge = 2,
    const bool print_voxel_summaries = true,
    const double min_pair_pt = 0.5,
    const unsigned int max_pair_records_per_voxel = 10,
    const std::string& crossing_solver = "helix",
    const double magnetic_field_z = 1.4,
    const bool write_pair_tree = true)
{
  const bool use_helix_solver = crossing_solver == "helix";
  if (crossing_solver != "line" && crossing_solver != "helix")
  {
    std::cerr << "CPM_QA_B1_ComputePoCA - invalid crossing_solver: "
              << crossing_solver << " (expected line or helix)" << std::endl;
    return;
  }

  const auto grid_metadata = CPMB1::load_grid_metadata(input_files);
  if (!grid_metadata.valid || !grid_metadata.consistent)
  {
    std::cerr << "CPM_QA_B1_ComputePoCA - refusing to merge inputs with invalid or inconsistent cpm_metadata: "
              << grid_metadata.message << std::endl;
    return;
  }

  std::map<CPMB1::VoxelKey, std::vector<CPMB1::Record>> records_by_voxel;

  const auto entries = CPMRecordReader::read_records(
      input_files,
      [&](const CPMRecordReader::Record& input_record)
      {
    CPMB1::Record record;
    record.entry = input_record.entry;
    record.event_track.cluster_source = input_record.cluster_source;
    record.event_track.track_source = input_record.track_source;
    record.event_track.run = input_record.run;
    record.event_track.segment = input_record.segment;
    record.event_track.sync_event = input_record.sync_event;
    record.event_track.event_sequence = input_record.event_sequence;
    record.event_track.stream_event_ordinal = input_record.stream_event_ordinal;
    record.event_track.track_id = input_record.track_id;
    record.cluskey = input_record.cluskey;
    record.charge = input_record.charge;
    record.pt = input_record.pt;
    record.voxel = {input_record.iphi, input_record.ir, input_record.iz};
    record.voxel_center = {
        input_record.voxel_center.x,
        input_record.voxel_center.y,
        input_record.voxel_center.z};
    record.offset = {
        input_record.offset.x,
        input_record.offset.y,
        input_record.offset.z};
    record.state_position = {
        input_record.state_position.x,
        input_record.state_position.y,
        input_record.state_position.z};
    record.state_momentum = {
        input_record.state_momentum.x,
        input_record.state_momentum.y,
        input_record.state_momentum.z};

    if (record.voxel.iphi < 0 || record.voxel.ir < 0 || record.voxel.iz < 0)
    {
      return;
    }

    records_by_voxel[record.voxel].push_back(record);
      });

  TFile output(output_file.c_str(), "RECREATE");

  TTree pairs("cpm_poca_pairs", "CPM crossing-point PoCA pairs");
  int out_iphi = -1;
  int out_ir = -1;
  int out_iz = -1;
  unsigned long long out_batch_id = 0;
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
  int solver_id = 0;
  bool used_line_fallback = false;
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
  pairs.Branch("batch_id", &out_batch_id);
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
  pairs.Branch("solver_id", &solver_id);
  pairs.Branch("used_line_fallback", &used_line_fallback);
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

  TTree batch_corrections("cpm_b1_batch_corrections", "CPM B1 batch-level crossing correction sums");
  int batch_iphi = -1;
  int batch_ir = -1;
  int batch_iz = -1;
  unsigned long long batch_id = 0;
  unsigned int batch_positive_records = 0;
  unsigned int batch_negative_records = 0;
  unsigned long long batch_record_pairs = 0;
  unsigned long long batch_opposite_charge_pairs = 0;
  unsigned long long batch_same_event_track_pairs = 0;
  unsigned long long batch_pt_rejected_pairs = 0;
  unsigned long long batch_candidate_pairs = 0;
  unsigned long long batch_invalid_weight_pairs = 0;
  unsigned long long batch_invalid_poca_pairs = 0;
  unsigned long long batch_dca_rejected_pairs = 0;
  unsigned long long batch_accepted_pairs = 0;
  double batch_voxel_x = std::numeric_limits<double>::quiet_NaN();
  double batch_voxel_y = std::numeric_limits<double>::quiet_NaN();
  double batch_voxel_z = std::numeric_limits<double>::quiet_NaN();
  double batch_sum_pair_weight = std::numeric_limits<double>::quiet_NaN();
  double batch_sum_pair_weight2 = std::numeric_limits<double>::quiet_NaN();
  double batch_mean_pair_weight = std::numeric_limits<double>::quiet_NaN();
  double batch_effective_pair_entries = std::numeric_limits<double>::quiet_NaN();
  double batch_sum_delta_r = 0.0;
  double batch_sum_delta_r2 = 0.0;
  double batch_sum_delta_rphi = 0.0;
  double batch_sum_delta_rphi2 = 0.0;
  double batch_sum_delta_phi = 0.0;
  double batch_sum_delta_phi2 = 0.0;
  double batch_sum_delta_z = 0.0;
  double batch_sum_delta_z2 = 0.0;
  double batch_sum_weighted_delta_r = 0.0;
  double batch_sum_weighted_delta_r2 = 0.0;
  double batch_sum_weighted_delta_rphi = 0.0;
  double batch_sum_weighted_delta_rphi2 = 0.0;
  double batch_sum_weighted_delta_phi = 0.0;
  double batch_sum_weighted_delta_phi2 = 0.0;
  double batch_sum_weighted_delta_z = 0.0;
  double batch_sum_weighted_delta_z2 = 0.0;
  double batch_sum_dca = 0.0;
  double batch_sum_dca2 = 0.0;
  double batch_mean_delta_r = std::numeric_limits<double>::quiet_NaN();
  double batch_rms_delta_r = std::numeric_limits<double>::quiet_NaN();
  double batch_mean_delta_rphi = std::numeric_limits<double>::quiet_NaN();
  double batch_rms_delta_rphi = std::numeric_limits<double>::quiet_NaN();
  double batch_mean_delta_phi = std::numeric_limits<double>::quiet_NaN();
  double batch_rms_delta_phi = std::numeric_limits<double>::quiet_NaN();
  double batch_mean_delta_z = std::numeric_limits<double>::quiet_NaN();
  double batch_rms_delta_z = std::numeric_limits<double>::quiet_NaN();
  double batch_weighted_mean_delta_r = std::numeric_limits<double>::quiet_NaN();
  double batch_weighted_rms_delta_r = std::numeric_limits<double>::quiet_NaN();
  double batch_weighted_mean_delta_rphi = std::numeric_limits<double>::quiet_NaN();
  double batch_weighted_rms_delta_rphi = std::numeric_limits<double>::quiet_NaN();
  double batch_weighted_mean_delta_phi = std::numeric_limits<double>::quiet_NaN();
  double batch_weighted_rms_delta_phi = std::numeric_limits<double>::quiet_NaN();
  double batch_weighted_mean_delta_z = std::numeric_limits<double>::quiet_NaN();
  double batch_weighted_rms_delta_z = std::numeric_limits<double>::quiet_NaN();
  double batch_mean_dca = std::numeric_limits<double>::quiet_NaN();
  double batch_rms_dca = std::numeric_limits<double>::quiet_NaN();

  batch_corrections.Branch("iphi", &batch_iphi);
  batch_corrections.Branch("ir", &batch_ir);
  batch_corrections.Branch("iz", &batch_iz);
  batch_corrections.Branch("batch_id", &batch_id);
  batch_corrections.Branch("positive_records", &batch_positive_records);
  batch_corrections.Branch("negative_records", &batch_negative_records);
  batch_corrections.Branch("record_pairs", &batch_record_pairs);
  batch_corrections.Branch("opposite_charge_pairs", &batch_opposite_charge_pairs);
  batch_corrections.Branch("same_event_track_pairs", &batch_same_event_track_pairs);
  batch_corrections.Branch("pt_rejected_pairs", &batch_pt_rejected_pairs);
  batch_corrections.Branch("candidate_pairs", &batch_candidate_pairs);
  batch_corrections.Branch("invalid_weight_pairs", &batch_invalid_weight_pairs);
  batch_corrections.Branch("invalid_poca_pairs", &batch_invalid_poca_pairs);
  batch_corrections.Branch("dca_rejected_pairs", &batch_dca_rejected_pairs);
  batch_corrections.Branch("accepted_pairs", &batch_accepted_pairs);
  batch_corrections.Branch("voxel_x", &batch_voxel_x);
  batch_corrections.Branch("voxel_y", &batch_voxel_y);
  batch_corrections.Branch("voxel_z", &batch_voxel_z);
  batch_corrections.Branch("sum_pair_weight", &batch_sum_pair_weight);
  batch_corrections.Branch("sum_pair_weight2", &batch_sum_pair_weight2);
  batch_corrections.Branch("mean_pair_weight", &batch_mean_pair_weight);
  batch_corrections.Branch("effective_pair_entries", &batch_effective_pair_entries);
  batch_corrections.Branch("sum_delta_r", &batch_sum_delta_r);
  batch_corrections.Branch("sum_delta_r2", &batch_sum_delta_r2);
  batch_corrections.Branch("sum_delta_rphi", &batch_sum_delta_rphi);
  batch_corrections.Branch("sum_delta_rphi2", &batch_sum_delta_rphi2);
  batch_corrections.Branch("sum_delta_phi", &batch_sum_delta_phi);
  batch_corrections.Branch("sum_delta_phi2", &batch_sum_delta_phi2);
  batch_corrections.Branch("sum_delta_z", &batch_sum_delta_z);
  batch_corrections.Branch("sum_delta_z2", &batch_sum_delta_z2);
  batch_corrections.Branch("sum_weighted_delta_r", &batch_sum_weighted_delta_r);
  batch_corrections.Branch("sum_weighted_delta_r2", &batch_sum_weighted_delta_r2);
  batch_corrections.Branch("sum_weighted_delta_rphi", &batch_sum_weighted_delta_rphi);
  batch_corrections.Branch("sum_weighted_delta_rphi2", &batch_sum_weighted_delta_rphi2);
  batch_corrections.Branch("sum_weighted_delta_phi", &batch_sum_weighted_delta_phi);
  batch_corrections.Branch("sum_weighted_delta_phi2", &batch_sum_weighted_delta_phi2);
  batch_corrections.Branch("sum_weighted_delta_z", &batch_sum_weighted_delta_z);
  batch_corrections.Branch("sum_weighted_delta_z2", &batch_sum_weighted_delta_z2);
  batch_corrections.Branch("sum_dca", &batch_sum_dca);
  batch_corrections.Branch("sum_dca2", &batch_sum_dca2);
  batch_corrections.Branch("mean_delta_r", &batch_mean_delta_r);
  batch_corrections.Branch("rms_delta_r", &batch_rms_delta_r);
  batch_corrections.Branch("mean_delta_rphi", &batch_mean_delta_rphi);
  batch_corrections.Branch("rms_delta_rphi", &batch_rms_delta_rphi);
  batch_corrections.Branch("mean_delta_phi", &batch_mean_delta_phi);
  batch_corrections.Branch("rms_delta_phi", &batch_rms_delta_phi);
  batch_corrections.Branch("mean_delta_z", &batch_mean_delta_z);
  batch_corrections.Branch("rms_delta_z", &batch_rms_delta_z);
  batch_corrections.Branch("weighted_mean_delta_r", &batch_weighted_mean_delta_r);
  batch_corrections.Branch("weighted_rms_delta_r", &batch_weighted_rms_delta_r);
  batch_corrections.Branch("weighted_mean_delta_rphi", &batch_weighted_mean_delta_rphi);
  batch_corrections.Branch("weighted_rms_delta_rphi", &batch_weighted_rms_delta_rphi);
  batch_corrections.Branch("weighted_mean_delta_phi", &batch_weighted_mean_delta_phi);
  batch_corrections.Branch("weighted_rms_delta_phi", &batch_weighted_rms_delta_phi);
  batch_corrections.Branch("weighted_mean_delta_z", &batch_weighted_mean_delta_z);
  batch_corrections.Branch("weighted_rms_delta_z", &batch_weighted_rms_delta_z);
  batch_corrections.Branch("mean_dca", &batch_mean_dca);
  batch_corrections.Branch("rms_dca", &batch_rms_dca);

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
  unsigned long long summary_duplicate_dropped_records = 0;
  unsigned long long summary_cap_dropped_records = 0;
  unsigned long long summary_batch_count = 0;
  unsigned long long summary_batched_positive_records = 0;
  unsigned long long summary_batched_negative_records = 0;
  unsigned long long summary_raw_record_pairs = 0;
  unsigned long long summary_selected_record_pairs = 0;
  unsigned long long summary_batched_record_pairs = 0;
  unsigned long long summary_same_charge_pairs = 0;
  unsigned long long summary_opposite_charge_pairs = 0;
  unsigned long long summary_batched_opposite_charge_pairs = 0;
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
  voxel_summaries.Branch("duplicate_dropped_records", &summary_duplicate_dropped_records);
  voxel_summaries.Branch("cap_dropped_records", &summary_cap_dropped_records);
  voxel_summaries.Branch("batches", &summary_batch_count);
  voxel_summaries.Branch("batched_positive_records", &summary_batched_positive_records);
  voxel_summaries.Branch("batched_negative_records", &summary_batched_negative_records);
  voxel_summaries.Branch("raw_record_pairs", &summary_raw_record_pairs);
  voxel_summaries.Branch("selected_record_pairs", &summary_selected_record_pairs);
  voxel_summaries.Branch("batched_record_pairs", &summary_batched_record_pairs);
  voxel_summaries.Branch("same_charge_pairs", &summary_same_charge_pairs);
  voxel_summaries.Branch("opposite_charge_pairs", &summary_opposite_charge_pairs);
  voxel_summaries.Branch("batched_opposite_charge_pairs", &summary_batched_opposite_charge_pairs);
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
    summary_duplicate_dropped_records = voxel_summary.duplicate_dropped_records;
    summary_cap_dropped_records = voxel_summary.cap_dropped_records;
    summary_batch_count = voxel_summary.batch_count;
    summary_batched_positive_records = voxel_summary.batched_positive_records;
    summary_batched_negative_records = voxel_summary.batched_negative_records;
    summary_raw_record_pairs = voxel_summary.raw_pairs;
    summary_selected_record_pairs = voxel_summary.selected_record_pairs;
    summary_batched_record_pairs = voxel_summary.batched_record_pairs;
    summary_same_charge_pairs = voxel_summary.same_charge_pairs;
    summary_opposite_charge_pairs = voxel_summary.opposite_charge_pairs;
    summary_batched_opposite_charge_pairs = voxel_summary.batched_opposite_charge_pairs;
    summary_candidate_pairs = voxel_summary.candidate_pairs;
    summary_accepted_pairs = voxel_summary.accepted_pairs;
    summary_positive_records = positive_records;
    summary_negative_records = negative_records;
    summary_status = status;
    voxel_summaries.Fill();
  };

  CPMPairOptions pair_options;
  pair_options.solver = use_helix_solver ? CPMPairSolver::Helix : CPMPairSolver::Line;
  pair_options.min_pt = min_pair_pt;
  pair_options.max_pair_dca = max_pair_dca;
  pair_options.magnetic_field_z = magnetic_field_z;
  pair_options.min_sin_angle = min_sin_angle;

  unsigned long long candidate_pairs = 0;
  unsigned long long accepted_pairs = 0;
  unsigned int processed_voxels = 0;
  unsigned int skipped_low_charge_voxels = 0;
  unsigned int skipped_large_voxels = 0;

  for (const auto& [voxel, records] : records_by_voxel)
  {
    CPMB1::VoxelSummary voxel_summary;
    std::map<CPMB1::UniqueTrackId, const CPMB1::Record*> closest_record_by_track;
    unsigned long long good_records = 0;

    for (const auto& record : records)
    {
      if (!cpmPairHasGoodPt(record.charge, record.pt, min_pair_pt))
      {
        ++voxel_summary.pt_rejected_records;
        continue;
      }
      ++good_records;
      const auto track_id = CPMB1::make_unique_track_id(record);
      auto [iter, inserted] = closest_record_by_track.emplace(track_id, &record);
      if (!inserted && CPMB1::is_closer_to_voxel_center(record, *iter->second))
      {
        iter->second = &record;
      }
    }

    voxel_summary.duplicate_dropped_records =
        good_records - static_cast<unsigned long long>(closest_record_by_track.size());

    std::vector<const CPMB1::Record*> selected_records;
    selected_records.reserve(closest_record_by_track.size());
    for (const auto& entry : closest_record_by_track)
    {
      selected_records.push_back(entry.second);
    }

    std::stable_sort(
        selected_records.begin(),
        selected_records.end(),
        [](const CPMB1::Record* lhs, const CPMB1::Record* rhs)
        {
          const auto lhs_hash = CPMB1::record_selection_hash(*lhs);
          const auto rhs_hash = CPMB1::record_selection_hash(*rhs);
          if (lhs_hash != rhs_hash)
          {
            return lhs_hash < rhs_hash;
          }
          if (lhs->event_track.track_id != rhs->event_track.track_id)
          {
            return lhs->event_track.track_id < rhs->event_track.track_id;
          }
          return lhs->cluskey < rhs->cluskey;
        });

    std::vector<const CPMB1::Record*> positive_selected_records;
    std::vector<const CPMB1::Record*> negative_selected_records;
    positive_selected_records.reserve(selected_records.size());
    negative_selected_records.reserve(selected_records.size());
    for (const auto* record : selected_records)
    {
      if (record->charge > 0)
      {
        positive_selected_records.push_back(record);
      }
      else if (record->charge < 0)
      {
        negative_selected_records.push_back(record);
      }
    }
    const auto positive_records =
        static_cast<unsigned int>(positive_selected_records.size());
    const auto negative_records =
        static_cast<unsigned int>(negative_selected_records.size());

    voxel_summary.unique_tracks = CPMB1::count_unique_tracks(records);
    voxel_summary.unique_track_pairs =
        CPMB1::pair_count(static_cast<std::size_t>(voxel_summary.unique_tracks));
    voxel_summary.selected_records = selected_records.size();
    voxel_summary.selected_unique_tracks = CPMB1::count_unique_tracks(selected_records);
    voxel_summary.raw_pairs = CPMB1::pair_count(records.size());
    voxel_summary.selected_record_pairs = CPMB1::pair_count(selected_records.size());
    voxel_summary.same_charge_pairs =
        CPMB1::pair_count(positive_records) + CPMB1::pair_count(negative_records);
    voxel_summary.opposite_charge_pairs =
        static_cast<unsigned long long>(positive_records) * negative_records;

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

    if (positive_records < min_records_per_charge ||
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
            "skipped_low_opposite_charge_records");
      }
      fill_voxel_summary(
          voxel,
          records.size(),
          positive_records,
          negative_records,
          voxel_summary,
          "skipped_low_opposite_charge_records");
      continue;
    }

    const std::size_t batch_charge_limit =
        max_pair_records_per_voxel > 0 ?
        static_cast<std::size_t>(max_pair_records_per_voxel) :
        std::numeric_limits<std::size_t>::max();

    std::size_t positive_index = 0;
    std::size_t negative_index = 0;
    while (positive_index < positive_selected_records.size() &&
           negative_index < negative_selected_records.size())
    {
      const std::size_t positive_remaining =
          positive_selected_records.size() - positive_index;
      const std::size_t negative_remaining =
          negative_selected_records.size() - negative_index;
      const std::size_t positive_take =
          std::min(batch_charge_limit, positive_remaining);
      const std::size_t negative_take =
          std::min(batch_charge_limit, negative_remaining);

      if (positive_take < min_records_per_charge ||
          negative_take < min_records_per_charge)
      {
        break;
      }

      const std::size_t positive_begin = positive_index;
      const std::size_t negative_begin = negative_index;
      positive_index += positive_take;
      negative_index += negative_take;

      const unsigned long long current_batch_id = voxel_summary.batch_count;
      ++voxel_summary.batch_count;
      voxel_summary.batched_positive_records += positive_take;
      voxel_summary.batched_negative_records += negative_take;
      voxel_summary.batched_record_pairs +=
          CPMB1::pair_count(positive_take + negative_take);
      voxel_summary.batched_opposite_charge_pairs +=
          static_cast<unsigned long long>(positive_take) * negative_take;

      CPMB1::BatchAccumulator batch_accumulator;
      for (std::size_t ipos = 0; ipos < positive_take; ++ipos)
      {
        for (std::size_t ineg = 0; ineg < negative_take; ++ineg)
        {
          const auto& a = *positive_selected_records[positive_begin + ipos];
          const auto& b = *negative_selected_records[negative_begin + ineg];
          if (a.event_track == b.event_track)
          {
            ++batch_accumulator.same_event_track_pairs;
            continue;
          }
          const CPMPairInput pair_input_a{
              a.charge,
              a.pt,
              a.state_position,
              a.state_momentum,
              a.offset};
          const CPMPairInput pair_input_b{
              b.charge,
              b.pt,
              b.state_position,
              b.state_momentum,
              b.offset};
          const auto pair_result = computeCPMPair(
              pair_input_a,
              pair_input_b,
              a.voxel_center,
              pair_options);

          if (pair_result.status == CPMPairStatus::PtRejected)
          {
            ++batch_accumulator.pt_rejected_pairs;
            continue;
          }

          ++candidate_pairs;
          ++voxel_summary.candidate_pairs;
          ++batch_accumulator.candidate_pairs;
          if (pair_result.status == CPMPairStatus::InvalidWeight)
          {
            ++batch_accumulator.invalid_weight_pairs;
            continue;
          }

          if (pair_result.status == CPMPairStatus::InvalidPoCA)
          {
            ++batch_accumulator.invalid_poca_pairs;
            continue;
          }
          if (pair_result.status == CPMPairStatus::DcaRejected)
          {
            ++batch_accumulator.dca_rejected_pairs;
            continue;
          }
          if (!pair_result.accepted())
          {
            continue;
          }

          out_iphi = voxel.iphi;
          out_ir = voxel.ir;
          out_iz = voxel.iz;
          out_batch_id = current_batch_id;
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
          pair_weight = pair_result.pair_weight;
          solver_id = use_helix_solver ? 1 : 0;
          used_line_fallback = pair_result.used_line_fallback;
          dca = pair_result.dca;
          s = pair_result.s;
          t = pair_result.t;
          point_ax = pair_result.point_a.x;
          point_ay = pair_result.point_a.y;
          point_az = pair_result.point_a.z;
          point_bx = pair_result.point_b.x;
          point_by = pair_result.point_b.y;
          point_bz = pair_result.point_b.z;
          midpoint_x = pair_result.midpoint.x;
          midpoint_y = pair_result.midpoint.y;
          midpoint_z = pair_result.midpoint.z;
          voxel_center_x = a.voxel_center.x;
          voxel_center_y = a.voxel_center.y;
          voxel_center_z = a.voxel_center.z;
          delta_x = pair_result.delta_x;
          delta_y = pair_result.delta_y;
          delta_z = pair_result.delta_z;
          delta_r = pair_result.delta_r;
          delta_phi = pair_result.delta_phi;
          delta_rphi = pair_result.delta_rphi;

          batch_accumulator.add(
              delta_r,
              delta_rphi,
              delta_phi,
              delta_z,
              dca,
              pair_weight,
              voxel_center_x,
              voxel_center_y,
              voxel_center_z);
          if (write_pair_tree)
          {
            pairs.Fill();
          }
          ++accepted_pairs;
          ++voxel_summary.accepted_pairs;
        }
      }

      if (batch_accumulator.accepted_pairs > 0)
      {
        batch_iphi = voxel.iphi;
        batch_ir = voxel.ir;
        batch_iz = voxel.iz;
        batch_id = current_batch_id;
        batch_positive_records = positive_take;
        batch_negative_records = negative_take;
        batch_record_pairs = CPMB1::pair_count(positive_take + negative_take);
        batch_opposite_charge_pairs =
            static_cast<unsigned long long>(positive_take) * negative_take;
        batch_same_event_track_pairs = batch_accumulator.same_event_track_pairs;
        batch_pt_rejected_pairs = batch_accumulator.pt_rejected_pairs;
        batch_candidate_pairs = batch_accumulator.candidate_pairs;
        batch_invalid_weight_pairs = batch_accumulator.invalid_weight_pairs;
        batch_invalid_poca_pairs = batch_accumulator.invalid_poca_pairs;
        batch_dca_rejected_pairs = batch_accumulator.dca_rejected_pairs;
        batch_accepted_pairs = batch_accumulator.accepted_pairs;
        batch_voxel_x = cpmCorrectionMean(batch_accumulator.sum_voxel_x, batch_accumulator.accepted_pairs);
        batch_voxel_y = cpmCorrectionMean(batch_accumulator.sum_voxel_y, batch_accumulator.accepted_pairs);
        batch_voxel_z = cpmCorrectionMean(batch_accumulator.sum_voxel_z, batch_accumulator.accepted_pairs);
        batch_sum_pair_weight = batch_accumulator.sum_pair_weight;
        batch_sum_pair_weight2 = batch_accumulator.sum_pair_weight2;
        batch_mean_pair_weight = cpmCorrectionMean(batch_accumulator.sum_pair_weight, batch_accumulator.accepted_pairs);
        batch_effective_pair_entries = cpmCorrectionEffectiveEntries(
            batch_accumulator.sum_pair_weight,
            batch_accumulator.sum_pair_weight2);
        batch_sum_delta_r = batch_accumulator.sum_delta_r;
        batch_sum_delta_r2 = batch_accumulator.sum_delta_r2;
        batch_sum_delta_rphi = batch_accumulator.sum_delta_rphi;
        batch_sum_delta_rphi2 = batch_accumulator.sum_delta_rphi2;
        batch_sum_delta_phi = batch_accumulator.sum_delta_phi;
        batch_sum_delta_phi2 = batch_accumulator.sum_delta_phi2;
        batch_sum_delta_z = batch_accumulator.sum_delta_z;
        batch_sum_delta_z2 = batch_accumulator.sum_delta_z2;
        batch_sum_weighted_delta_r = batch_accumulator.sum_weighted_delta_r;
        batch_sum_weighted_delta_r2 = batch_accumulator.sum_weighted_delta_r2;
        batch_sum_weighted_delta_rphi = batch_accumulator.sum_weighted_delta_rphi;
        batch_sum_weighted_delta_rphi2 = batch_accumulator.sum_weighted_delta_rphi2;
        batch_sum_weighted_delta_phi = batch_accumulator.sum_weighted_delta_phi;
        batch_sum_weighted_delta_phi2 = batch_accumulator.sum_weighted_delta_phi2;
        batch_sum_weighted_delta_z = batch_accumulator.sum_weighted_delta_z;
        batch_sum_weighted_delta_z2 = batch_accumulator.sum_weighted_delta_z2;
        batch_sum_dca = batch_accumulator.sum_dca;
        batch_sum_dca2 = batch_accumulator.sum_dca2;
        batch_mean_delta_r = cpmCorrectionMean(batch_accumulator.sum_delta_r, batch_accumulator.accepted_pairs);
        batch_rms_delta_r = cpmCorrectionRms(batch_accumulator.sum_delta_r, batch_accumulator.sum_delta_r2, batch_accumulator.accepted_pairs);
        batch_mean_delta_rphi = cpmCorrectionMean(batch_accumulator.sum_delta_rphi, batch_accumulator.accepted_pairs);
        batch_rms_delta_rphi = cpmCorrectionRms(batch_accumulator.sum_delta_rphi, batch_accumulator.sum_delta_rphi2, batch_accumulator.accepted_pairs);
        batch_mean_delta_phi = cpmCorrectionMean(batch_accumulator.sum_delta_phi, batch_accumulator.accepted_pairs);
        batch_rms_delta_phi = cpmCorrectionRms(batch_accumulator.sum_delta_phi, batch_accumulator.sum_delta_phi2, batch_accumulator.accepted_pairs);
        batch_mean_delta_z = cpmCorrectionMean(batch_accumulator.sum_delta_z, batch_accumulator.accepted_pairs);
        batch_rms_delta_z = cpmCorrectionRms(batch_accumulator.sum_delta_z, batch_accumulator.sum_delta_z2, batch_accumulator.accepted_pairs);
        batch_weighted_mean_delta_r = cpmCorrectionWeightedMean(batch_accumulator.sum_weighted_delta_r, batch_accumulator.sum_pair_weight);
        batch_weighted_rms_delta_r = cpmCorrectionWeightedRms(batch_accumulator.sum_weighted_delta_r, batch_accumulator.sum_weighted_delta_r2, batch_accumulator.sum_pair_weight);
        batch_weighted_mean_delta_rphi = cpmCorrectionWeightedMean(batch_accumulator.sum_weighted_delta_rphi, batch_accumulator.sum_pair_weight);
        batch_weighted_rms_delta_rphi = cpmCorrectionWeightedRms(batch_accumulator.sum_weighted_delta_rphi, batch_accumulator.sum_weighted_delta_rphi2, batch_accumulator.sum_pair_weight);
        batch_weighted_mean_delta_phi = cpmCorrectionWeightedMean(batch_accumulator.sum_weighted_delta_phi, batch_accumulator.sum_pair_weight);
        batch_weighted_rms_delta_phi = cpmCorrectionWeightedRms(batch_accumulator.sum_weighted_delta_phi, batch_accumulator.sum_weighted_delta_phi2, batch_accumulator.sum_pair_weight);
        batch_weighted_mean_delta_z = cpmCorrectionWeightedMean(batch_accumulator.sum_weighted_delta_z, batch_accumulator.sum_pair_weight);
        batch_weighted_rms_delta_z = cpmCorrectionWeightedRms(batch_accumulator.sum_weighted_delta_z, batch_accumulator.sum_weighted_delta_z2, batch_accumulator.sum_pair_weight);
        batch_mean_dca = cpmCorrectionMean(batch_accumulator.sum_dca, batch_accumulator.accepted_pairs);
        batch_rms_dca = cpmCorrectionRms(batch_accumulator.sum_dca, batch_accumulator.sum_dca2, batch_accumulator.accepted_pairs);
        batch_corrections.Fill();
      }
    }

    if (voxel_summary.batch_count == 0)
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
            "skipped_low_batch_charge_records");
      }
      fill_voxel_summary(
          voxel,
          records.size(),
          positive_records,
          negative_records,
          voxel_summary,
          "skipped_low_batch_charge_records");
      continue;
    }

    ++processed_voxels;

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

  TTree summary("cpm_b1_summary", "CPM B1 crossing-point PoCA summary");
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
  unsigned int summary_max_pair_records_per_charge_batch = max_pair_records_per_voxel;
  std::string summary_crossing_solver = crossing_solver;
  bool summary_use_helix_solver = use_helix_solver;
  double summary_magnetic_field_z = magnetic_field_z;
  bool summary_write_pair_tree = write_pair_tree;
  Long64_t summary_output_pair_rows = pairs.GetEntries();
  Long64_t summary_output_batch_rows = batch_corrections.GetEntries();
  int summary_phi_bins = grid_metadata.phi_bins;
  int summary_r_bins = grid_metadata.r_bins;
  int summary_z_bins = grid_metadata.z_bins;
  double summary_r_min = grid_metadata.r_min;
  double summary_r_max = grid_metadata.r_max;
  double summary_z_min = grid_metadata.z_min;
  double summary_z_max = grid_metadata.z_max;
  bool summary_grid_metadata_valid = grid_metadata.valid;
  bool summary_grid_metadata_consistent = grid_metadata.consistent;

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
  summary.Branch("max_pair_records_per_charge_batch", &summary_max_pair_records_per_charge_batch);
  summary.Branch("crossing_solver", &summary_crossing_solver);
  summary.Branch("use_helix_solver", &summary_use_helix_solver);
  summary.Branch("magnetic_field_z", &summary_magnetic_field_z);
  summary.Branch("write_pair_tree", &summary_write_pair_tree);
  summary.Branch("output_pair_rows", &summary_output_pair_rows);
  summary.Branch("output_batch_rows", &summary_output_batch_rows);
  summary.Branch("phi_bins", &summary_phi_bins);
  summary.Branch("r_bins", &summary_r_bins);
  summary.Branch("z_bins", &summary_z_bins);
  summary.Branch("r_min", &summary_r_min);
  summary.Branch("r_max", &summary_r_max);
  summary.Branch("z_min", &summary_z_min);
  summary.Branch("z_max", &summary_z_max);
  summary.Branch("grid_metadata_valid", &summary_grid_metadata_valid);
  summary.Branch("grid_metadata_consistent", &summary_grid_metadata_consistent);
  summary.Fill();

  pairs.Write();
  batch_corrections.Write();
  voxel_summaries.Write();
  summary.Write();
  output.Close();

  std::cout << "CPM_QA_B1_ComputePoCA - input records: " << input_records << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - voxels: " << voxel_count << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - processed voxels: " << processed_voxels << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - skipped large voxels: " << skipped_large_voxels << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - skipped low-charge voxels: " << skipped_low_charge_voxels << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - candidate pairs: " << candidate_pairs << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - accepted pairs: " << accepted_pairs << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - output pair rows: " << summary_output_pair_rows << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - output batch rows: " << summary_output_batch_rows << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - min pair pt: " << min_pair_pt << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - max pair records per charge batch: " << max_pair_records_per_voxel << std::endl;
  std::cout << "CPM_QA_B1_ComputePoCA - crossing solver: " << crossing_solver << std::endl;
  if (use_helix_solver)
  {
    std::cout << "CPM_QA_B1_ComputePoCA - magnetic field z: " << magnetic_field_z << std::endl;
  }
  if (grid_metadata.valid)
  {
    std::cout << "CPM_QA_B1_ComputePoCA - grid bins: ("
              << grid_metadata.phi_bins << ", "
              << grid_metadata.r_bins << ", "
              << grid_metadata.z_bins << ")"
              << " r=[" << grid_metadata.r_min << ", " << grid_metadata.r_max << "]"
              << " z=[" << grid_metadata.z_min << ", " << grid_metadata.z_max << "]"
              << std::endl;
  }
  std::cout << "CPM_QA_B1_ComputePoCA - output: " << output_file << std::endl;
}

void CPM_QA_B1_ComputePoCA(
    const std::string& input_file,
    const std::string& output_file = "CPM_QA_B1_poca.root",
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 0,
    const unsigned int min_records_per_charge = 2,
    const bool print_voxel_summaries = true,
    const double min_pair_pt = 0.5,
    const unsigned int max_pair_records_per_voxel = 10,
    const std::string& crossing_solver = "helix",
    const double magnetic_field_z = 1.4,
    const bool write_pair_tree = true)
{
  CPM_QA_B1_ComputePoCA(
      std::vector<std::string>{input_file},
      output_file,
      max_pair_dca,
      min_sin_angle,
      max_records_per_voxel,
      min_records_per_charge,
      print_voxel_summaries,
      min_pair_pt,
      max_pair_records_per_voxel,
      crossing_solver,
      magnetic_field_z,
      write_pair_tree);
}

void CPM_QA_B1_ComputePoCA(
    const std::string& input_file_or_list,
    const std::string& output_file,
    const bool input_is_list,
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 0,
    const unsigned int min_records_per_charge = 2,
    const bool print_voxel_summaries = true,
    const double min_pair_pt = 0.5,
    const unsigned int max_pair_records_per_voxel = 10,
    const std::string& crossing_solver = "helix",
    const double magnetic_field_z = 1.4,
    const bool write_pair_tree = true)
{
  const auto input_files = input_is_list ?
      CPMB1::read_file_list(input_file_or_list) :
      std::vector<std::string>{input_file_or_list};

  CPM_QA_B1_ComputePoCA(
      input_files,
      output_file,
      max_pair_dca,
      min_sin_angle,
      max_records_per_voxel,
      min_records_per_charge,
      print_voxel_summaries,
      min_pair_pt,
      max_pair_records_per_voxel,
      crossing_solver,
      magnetic_field_z,
      write_pair_tree);
}
