/*
 * CPM production Job B average-correction calculation.
 *
 * This macro is the compact production path:
 *   Job A cpm_records -> opposite-charge crossing pairs -> voxel averages
 *   -> guarded average-correction histograms.
 *
 * It intentionally does not write pair trees, batch trees, or detailed
 * per-voxel QA. Use the CPM_QA_* macros when those diagnostic products are
 * needed.
 */

#include <CPMCorrectionAccumulator.h>
#include <CPMPairUtils.h>

#include <tpccalib/TpcSpaceChargeReconstructionHelper.h>

#include <TChain.h>
#include <TFile.h>
#include <TH3.h>
#include <TH3F.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

R__LOAD_LIBRARY(libcpm.so)
R__LOAD_LIBRARY(libtpccalib.so)

namespace CPMProduction
{
  constexpr double kPi = 3.14159265358979323846;

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
    double phi_min = 0.0;
    double phi_max = 2.0 * kPi;
    double r_min = std::numeric_limits<double>::quiet_NaN();
    double r_max = std::numeric_limits<double>::quiet_NaN();
    double z_min = std::numeric_limits<double>::quiet_NaN();
    double z_max = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
    bool consistent = true;
    std::string source_file;
    std::string message;
  };

  struct CalculationSummary
  {
    unsigned int input_files = 0;
    unsigned long long input_records = 0;
    unsigned int input_voxels = 0;
    unsigned int processed_voxels = 0;
    unsigned int skipped_large_voxels = 0;
    unsigned int skipped_low_selected_voxels = 0;
    unsigned int skipped_low_charge_voxels = 0;
    unsigned int skipped_low_batch_charge_voxels = 0;
    unsigned long long pt_rejected_records = 0;
    unsigned long long duplicate_dropped_records = 0;
    unsigned long long candidate_pairs = 0;
    unsigned long long accepted_pairs = 0;
    unsigned long long same_event_track_pairs = 0;
    unsigned long long invalid_weight_pairs = 0;
    unsigned long long invalid_poca_pairs = 0;
    unsigned long long dca_rejected_pairs = 0;
    unsigned int accumulator_voxels = 0;
    unsigned int filled_voxels = 0;
    unsigned int skipped_low_entry_voxels = 0;
    unsigned int skipped_invalid_voxels = 0;
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
    if (std::isfinite(candidate_distance2) &&
        std::isfinite(current_distance2) &&
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

  std::tuple<std::unique_ptr<TH3>, std::unique_ptr<TH3>> finish_histogram(
      TH3* source,
      const TString& name)
  {
    const auto result = TpcSpaceChargeReconstructionHelper::split(source);
    std::unique_ptr<TH3> hneg(std::get<0>(result));
    std::unique_ptr<TH3> hpos(std::get<1>(result));

    return std::make_tuple(
        std::unique_ptr<TH3>(TpcSpaceChargeReconstructionHelper::add_guarding_bins(hneg.get(), name + "_negz")),
        std::unique_ptr<TH3>(TpcSpaceChargeReconstructionHelper::add_guarding_bins(hpos.get(), name + "_posz")));
  }

  void write_finished_histogram(TFile& output, TH3* source, const TString& name)
  {
    auto finished = finish_histogram(source, name);
    output.cd();
    std::get<0>(finished)->Write();
    std::get<1>(finished)->Write();
  }
}

bool CPM_ComputeAverageCorrection(
    const std::vector<std::string>& input_files,
    const std::string& output_file = "CPM_B3_average_correction_histograms.root",
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1,
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 0,
    const unsigned int min_records_per_charge = 2,
    const double min_pair_pt = 0.5,
    const unsigned int max_pair_records_per_charge_batch = 10,
    const std::string& crossing_solver = "helix",
    const double magnetic_field_z = 1.4,
    const std::string& metadata_file = "")
{
  const bool use_helix_solver = crossing_solver == "helix";
  if (crossing_solver != "line" && crossing_solver != "helix")
  {
    std::cerr << "CPM_ComputeAverageCorrection - invalid crossing_solver: "
              << crossing_solver << " (expected line or helix)" << std::endl;
    return false;
  }

  auto grid_metadata = CPMProduction::load_grid_metadata(input_files);
  if (!grid_metadata.valid || !grid_metadata.consistent)
  {
    std::cerr << "CPM_ComputeAverageCorrection - refusing to merge inputs with invalid or inconsistent cpm_metadata: "
              << grid_metadata.message << std::endl;
    return false;
  }
  if (!metadata_file.empty())
  {
    const auto requested_metadata = CPMProduction::load_grid_metadata_file(metadata_file);
    if (!requested_metadata.valid ||
        !CPMProduction::same_grid_metadata(grid_metadata, requested_metadata))
    {
      std::cerr << "CPM_ComputeAverageCorrection - provided metadata file is invalid or inconsistent with input records: "
                << metadata_file << std::endl;
      return false;
    }
  }

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

  std::map<CPMProduction::VoxelKey, std::vector<CPMProduction::Record>> records_by_voxel;

  CPMProduction::CalculationSummary summary;
  summary.input_files = input_files.size();
  const auto entries = chain.GetEntries();
  summary.input_records = entries;
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    chain.GetEntry(entry);

    CPMProduction::Record record;
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
  summary.input_voxels = records_by_voxel.size();

  CPMPairOptions pair_options;
  pair_options.solver = use_helix_solver ? CPMPairSolver::Helix : CPMPairSolver::Line;
  pair_options.min_pt = min_pair_pt;
  pair_options.max_pair_dca = max_pair_dca;
  pair_options.magnetic_field_z = magnetic_field_z;
  pair_options.min_sin_angle = min_sin_angle;

  std::map<CPMProduction::VoxelKey, CPMCorrectionAccumulator> accumulators;

  for (const auto& [voxel, records] : records_by_voxel)
  {
    if (max_records_per_voxel > 0 && records.size() > max_records_per_voxel)
    {
      ++summary.skipped_large_voxels;
      continue;
    }

    std::map<CPMProduction::UniqueTrackId, const CPMProduction::Record*> closest_record_by_track;
    unsigned long long good_records = 0;
    for (const auto& record : records)
    {
      if (!cpmPairHasGoodPt(record.charge, record.pt, min_pair_pt))
      {
        ++summary.pt_rejected_records;
        continue;
      }
      ++good_records;
      const auto track_id_key = CPMProduction::make_unique_track_id(record);
      auto [iter, inserted] = closest_record_by_track.emplace(track_id_key, &record);
      if (!inserted && CPMProduction::is_closer_to_voxel_center(record, *iter->second))
      {
        iter->second = &record;
      }
    }
    summary.duplicate_dropped_records +=
        good_records - static_cast<unsigned long long>(closest_record_by_track.size());

    std::vector<const CPMProduction::Record*> selected_records;
    selected_records.reserve(closest_record_by_track.size());
    for (const auto& entry : closest_record_by_track)
    {
      selected_records.push_back(entry.second);
    }

    if (selected_records.size() < 2)
    {
      ++summary.skipped_low_selected_voxels;
      continue;
    }

    std::stable_sort(
        selected_records.begin(),
        selected_records.end(),
        [](const CPMProduction::Record* lhs, const CPMProduction::Record* rhs)
        {
          const auto lhs_hash = CPMProduction::record_selection_hash(*lhs);
          const auto rhs_hash = CPMProduction::record_selection_hash(*rhs);
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

    std::vector<const CPMProduction::Record*> positive_selected_records;
    std::vector<const CPMProduction::Record*> negative_selected_records;
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

    if (positive_selected_records.size() < min_records_per_charge ||
        negative_selected_records.size() < min_records_per_charge)
    {
      ++summary.skipped_low_charge_voxels;
      continue;
    }

    const std::size_t batch_charge_limit =
        max_pair_records_per_charge_batch > 0 ?
        static_cast<std::size_t>(max_pair_records_per_charge_batch) :
        std::numeric_limits<std::size_t>::max();

    bool accepted_any_batch = false;
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

      for (std::size_t ipos = 0; ipos < positive_take; ++ipos)
      {
        for (std::size_t ineg = 0; ineg < negative_take; ++ineg)
        {
          const auto& a = *positive_selected_records[positive_begin + ipos];
          const auto& b = *negative_selected_records[negative_begin + ineg];
          if (a.event_track == b.event_track)
          {
            ++summary.same_event_track_pairs;
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
            continue;
          }

          ++summary.candidate_pairs;
          if (pair_result.status == CPMPairStatus::InvalidWeight)
          {
            ++summary.invalid_weight_pairs;
            continue;
          }
          if (pair_result.status == CPMPairStatus::InvalidPoCA)
          {
            ++summary.invalid_poca_pairs;
            continue;
          }
          if (pair_result.status == CPMPairStatus::DcaRejected)
          {
            ++summary.dca_rejected_pairs;
            continue;
          }
          if (!pair_result.accepted())
          {
            continue;
          }

          const double accumulation_weight =
              use_pair_weights ? pair_result.pair_weight : 1.0;
          accumulators[voxel].add(
              pair_result.delta_r,
              pair_result.delta_rphi,
              pair_result.delta_phi,
              pair_result.delta_z,
              pair_result.dca,
              accumulation_weight,
              a.voxel_center.x,
              a.voxel_center.y,
              a.voxel_center.z);

          ++summary.accepted_pairs;
          accepted_any_batch = true;
        }
      }
    }

    if (accepted_any_batch)
    {
      ++summary.processed_voxels;
    }
    else
    {
      ++summary.skipped_low_batch_charge_voxels;
    }
  }

  summary.accumulator_voxels = accumulators.size();

  TH3F hentries_rec(
      "hentries_rec",
      "CPM voxel entries;phi;r;z",
      grid_metadata.phi_bins, grid_metadata.phi_min, grid_metadata.phi_max,
      grid_metadata.r_bins, grid_metadata.r_min, grid_metadata.r_max,
      grid_metadata.z_bins, grid_metadata.z_min, grid_metadata.z_max);
  TH3F hdistortion_r_rec(
      "hDistortionR_rec",
      "CPM radial distortion;phi;r;z",
      grid_metadata.phi_bins, grid_metadata.phi_min, grid_metadata.phi_max,
      grid_metadata.r_bins, grid_metadata.r_min, grid_metadata.r_max,
      grid_metadata.z_bins, grid_metadata.z_min, grid_metadata.z_max);
  TH3F hdistortion_p_rec(
      "hDistortionP_rec",
      "CPM phi distortion;phi;r;z",
      grid_metadata.phi_bins, grid_metadata.phi_min, grid_metadata.phi_max,
      grid_metadata.r_bins, grid_metadata.r_min, grid_metadata.r_max,
      grid_metadata.z_bins, grid_metadata.z_min, grid_metadata.z_max);
  TH3F hdistortion_z_rec(
      "hDistortionZ_rec",
      "CPM z distortion;phi;r;z",
      grid_metadata.phi_bins, grid_metadata.phi_min, grid_metadata.phi_max,
      grid_metadata.r_bins, grid_metadata.r_min, grid_metadata.r_max,
      grid_metadata.z_bins, grid_metadata.z_min, grid_metadata.z_max);

  for (const auto& [voxel, accumulator] : accumulators)
  {
    if (accumulator.entries < min_entries_per_voxel)
    {
      ++summary.skipped_low_entry_voxels;
      continue;
    }
    if (voxel.iphi < 0 || voxel.iphi >= grid_metadata.phi_bins ||
        voxel.ir < 0 || voxel.ir >= grid_metadata.r_bins ||
        voxel.iz < 0 || voxel.iz >= grid_metadata.z_bins)
    {
      ++summary.skipped_invalid_voxels;
      continue;
    }

    const double mean_delta_r = cpmCorrectionWeightedMean(
        accumulator.sum_weighted_delta_r,
        accumulator.sum_pair_weight);
    const double mean_delta_phi = cpmCorrectionWeightedMean(
        accumulator.sum_weighted_delta_phi,
        accumulator.sum_pair_weight);
    const double mean_delta_z = cpmCorrectionWeightedMean(
        accumulator.sum_weighted_delta_z,
        accumulator.sum_pair_weight);

    if (!std::isfinite(mean_delta_r) ||
        !std::isfinite(mean_delta_phi) ||
        !std::isfinite(mean_delta_z))
    {
      ++summary.skipped_invalid_voxels;
      continue;
    }

    hentries_rec.SetBinContent(
        voxel.iphi + 1,
        voxel.ir + 1,
        voxel.iz + 1,
        static_cast<double>(accumulator.entries));
    hdistortion_r_rec.SetBinContent(voxel.iphi + 1, voxel.ir + 1, voxel.iz + 1, mean_delta_r);
    hdistortion_p_rec.SetBinContent(voxel.iphi + 1, voxel.ir + 1, voxel.iz + 1, mean_delta_phi);
    hdistortion_z_rec.SetBinContent(voxel.iphi + 1, voxel.ir + 1, voxel.iz + 1, mean_delta_z);
    ++summary.filled_voxels;
  }

  TFile output(output_file.c_str(), "RECREATE");
  if (output.IsZombie())
  {
    std::cerr << "CPM_ComputeAverageCorrection - could not create output: "
              << output_file << std::endl;
    return false;
  }

  hentries_rec.Write();
  hdistortion_r_rec.Write();
  hdistortion_p_rec.Write();
  hdistortion_z_rec.Write();

  CPMProduction::write_finished_histogram(output, &hentries_rec, "hentries");
  CPMProduction::write_finished_histogram(output, &hdistortion_r_rec, "hIntDistortionR");
  CPMProduction::write_finished_histogram(output, &hdistortion_p_rec, "hIntDistortionP");
  CPMProduction::write_finished_histogram(output, &hdistortion_z_rec, "hIntDistortionZ");

  TTree summary_tree("cpm_b3_summary", "CPM production average-correction summary");
  unsigned int input_files_count = summary.input_files;
  unsigned long long input_records = summary.input_records;
  unsigned int input_voxels = summary.input_voxels;
  unsigned int processed_voxels = summary.processed_voxels;
  unsigned int skipped_large_voxels = summary.skipped_large_voxels;
  unsigned int skipped_low_selected_voxels = summary.skipped_low_selected_voxels;
  unsigned int skipped_low_charge_voxels = summary.skipped_low_charge_voxels;
  unsigned int skipped_low_batch_charge_voxels = summary.skipped_low_batch_charge_voxels;
  unsigned long long pt_rejected_records = summary.pt_rejected_records;
  unsigned long long duplicate_dropped_records = summary.duplicate_dropped_records;
  unsigned long long candidate_pairs = summary.candidate_pairs;
  unsigned long long accepted_pairs = summary.accepted_pairs;
  unsigned long long same_event_track_pairs = summary.same_event_track_pairs;
  unsigned long long invalid_weight_pairs = summary.invalid_weight_pairs;
  unsigned long long invalid_poca_pairs = summary.invalid_poca_pairs;
  unsigned long long dca_rejected_pairs = summary.dca_rejected_pairs;
  unsigned int accumulator_voxels = summary.accumulator_voxels;
  unsigned int filled_voxels = summary.filled_voxels;
  unsigned int skipped_low_entry_voxels = summary.skipped_low_entry_voxels;
  unsigned int skipped_invalid_voxels = summary.skipped_invalid_voxels;
  unsigned int summary_min_entries_per_voxel = min_entries_per_voxel;
  double summary_max_pair_dca = max_pair_dca;
  double summary_min_sin_angle = min_sin_angle;
  unsigned int summary_max_records_per_voxel = max_records_per_voxel;
  unsigned int summary_min_records_per_charge = min_records_per_charge;
  double summary_min_pair_pt = min_pair_pt;
  unsigned int summary_max_pair_records_per_charge_batch = max_pair_records_per_charge_batch;
  bool summary_use_pair_weights = use_pair_weights;
  std::string summary_averaging_mode = use_pair_weights ? "weighted" : "plain";
  std::string summary_crossing_solver = crossing_solver;
  double summary_magnetic_field_z = magnetic_field_z;
  int phi_bins = grid_metadata.phi_bins;
  int r_bins = grid_metadata.r_bins;
  int z_bins = grid_metadata.z_bins;
  double phi_min = grid_metadata.phi_min;
  double phi_max = grid_metadata.phi_max;
  double r_min = grid_metadata.r_min;
  double r_max = grid_metadata.r_max;
  double z_min = grid_metadata.z_min;
  double z_max = grid_metadata.z_max;

  summary_tree.Branch("input_files", &input_files_count);
  summary_tree.Branch("input_records", &input_records);
  summary_tree.Branch("input_voxels", &input_voxels);
  summary_tree.Branch("processed_voxels", &processed_voxels);
  summary_tree.Branch("skipped_large_voxels", &skipped_large_voxels);
  summary_tree.Branch("skipped_low_selected_voxels", &skipped_low_selected_voxels);
  summary_tree.Branch("skipped_low_charge_voxels", &skipped_low_charge_voxels);
  summary_tree.Branch("skipped_low_batch_charge_voxels", &skipped_low_batch_charge_voxels);
  summary_tree.Branch("pt_rejected_records", &pt_rejected_records);
  summary_tree.Branch("duplicate_dropped_records", &duplicate_dropped_records);
  summary_tree.Branch("candidate_pairs", &candidate_pairs);
  summary_tree.Branch("accepted_pairs", &accepted_pairs);
  summary_tree.Branch("same_event_track_pairs", &same_event_track_pairs);
  summary_tree.Branch("invalid_weight_pairs", &invalid_weight_pairs);
  summary_tree.Branch("invalid_poca_pairs", &invalid_poca_pairs);
  summary_tree.Branch("dca_rejected_pairs", &dca_rejected_pairs);
  summary_tree.Branch("accumulator_voxels", &accumulator_voxels);
  summary_tree.Branch("filled_voxels", &filled_voxels);
  summary_tree.Branch("skipped_low_entry_voxels", &skipped_low_entry_voxels);
  summary_tree.Branch("skipped_invalid_voxels", &skipped_invalid_voxels);
  summary_tree.Branch("min_entries_per_voxel", &summary_min_entries_per_voxel);
  summary_tree.Branch("max_pair_dca", &summary_max_pair_dca);
  summary_tree.Branch("min_sin_angle", &summary_min_sin_angle);
  summary_tree.Branch("max_records_per_voxel", &summary_max_records_per_voxel);
  summary_tree.Branch("min_records_per_charge", &summary_min_records_per_charge);
  summary_tree.Branch("min_pair_pt", &summary_min_pair_pt);
  summary_tree.Branch("max_pair_records_per_charge_batch", &summary_max_pair_records_per_charge_batch);
  summary_tree.Branch("use_pair_weights", &summary_use_pair_weights);
  summary_tree.Branch("averaging_mode", &summary_averaging_mode);
  summary_tree.Branch("crossing_solver", &summary_crossing_solver);
  summary_tree.Branch("magnetic_field_z", &summary_magnetic_field_z);
  summary_tree.Branch("phi_bins", &phi_bins);
  summary_tree.Branch("r_bins", &r_bins);
  summary_tree.Branch("z_bins", &z_bins);
  summary_tree.Branch("phi_min", &phi_min);
  summary_tree.Branch("phi_max", &phi_max);
  summary_tree.Branch("r_min", &r_min);
  summary_tree.Branch("r_max", &r_max);
  summary_tree.Branch("z_min", &z_min);
  summary_tree.Branch("z_max", &z_max);
  summary_tree.Fill();
  summary_tree.Write();

  output.Close();

  std::cout << "CPM_ComputeAverageCorrection - input files: " << input_files_count << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - input records: " << input_records << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - input voxels: " << input_voxels << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - processed voxels: " << processed_voxels << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - candidate pairs: " << candidate_pairs << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - accepted pairs: " << accepted_pairs << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - filled voxels: " << filled_voxels << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - averaging mode: " << summary_averaging_mode << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - crossing solver: " << crossing_solver << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - output: " << output_file << std::endl;

  return true;
}

bool CPM_ComputeAverageCorrection(
    const std::string& input_file_or_list,
    const std::string& output_file = "CPM_B3_average_correction_histograms.root",
    const bool input_is_list = false,
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1,
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 0,
    const unsigned int min_records_per_charge = 2,
    const double min_pair_pt = 0.5,
    const unsigned int max_pair_records_per_charge_batch = 10,
    const std::string& crossing_solver = "helix",
    const double magnetic_field_z = 1.4,
    const std::string& metadata_file = "")
{
  const auto input_files = input_is_list ?
      CPMProduction::read_file_list(input_file_or_list) :
      std::vector<std::string>{input_file_or_list};

  return CPM_ComputeAverageCorrection(
      input_files,
      output_file,
      use_pair_weights,
      min_entries_per_voxel,
      max_pair_dca,
      min_sin_angle,
      max_records_per_voxel,
      min_records_per_charge,
      min_pair_pt,
      max_pair_records_per_charge_batch,
      crossing_solver,
      magnetic_field_z,
      metadata_file);
}
