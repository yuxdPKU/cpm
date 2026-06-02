#include "CPMAverageCorrectionReconstruction.h"

#include "CPMReconstructionHelper.h"

#include <TFile.h>
#include <TH3.h>
#include <TH3F.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string_view>
#include <tuple>
#include <vector>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace
{
  void trim_heap_after_chunk()
  {
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
  }

  std::string_view view_string(const std::string& value)
  {
    return {value.data(), value.size()};
  }

  struct EventTrackKeyView
  {
    std::string_view cluster_source;
    std::string_view track_source;
    int run = -1;
    int segment = -1;
    int sync_event = -1;
    int event_sequence = -1;
    unsigned long long stream_event_ordinal = 0;
    unsigned int track_id = 0;

    friend bool operator<(const EventTrackKeyView& lhs, const EventTrackKeyView& rhs)
    {
      return std::tie(lhs.cluster_source, lhs.track_source, lhs.run, lhs.segment,
                      lhs.sync_event, lhs.event_sequence, lhs.stream_event_ordinal,
                      lhs.track_id) <
             std::tie(rhs.cluster_source, rhs.track_source, rhs.run, rhs.segment,
                      rhs.sync_event, rhs.event_sequence, rhs.stream_event_ordinal,
                      rhs.track_id);
    }
  };

  EventTrackKeyView event_track_key_view(const TrackStateRecord& record)
  {
    return {
        view_string(record.event_ref.cluster_source),
        view_string(record.event_ref.track_source),
        record.event_ref.run,
        record.event_ref.segment,
        record.event_ref.sync_event,
        record.event_ref.event_sequence,
        record.event_ref.stream_event_ordinal,
        record.track_ref.track_id};
  }

  bool same_event_track(const TrackStateRecord& lhs, const TrackStateRecord& rhs)
  {
    return lhs.event_ref.cluster_source == rhs.event_ref.cluster_source &&
           lhs.event_ref.track_source == rhs.event_ref.track_source &&
           lhs.event_ref.run == rhs.event_ref.run &&
           lhs.event_ref.segment == rhs.event_ref.segment &&
           lhs.event_ref.sync_event == rhs.event_ref.sync_event &&
           lhs.event_ref.event_sequence == rhs.event_ref.event_sequence &&
           lhs.event_ref.stream_event_ordinal == rhs.event_ref.stream_event_ordinal &&
           lhs.track_ref.track_id == rhs.track_ref.track_id;
  }
}

CPMAverageCorrectionReconstruction::CPMAverageCorrectionReconstruction()
{
}

CPMAverageCorrectionReconstruction::~CPMAverageCorrectionReconstruction() = default;

bool CPMAverageCorrectionReconstruction::add(const CPMVoxelContainer& source)
{
  if (!source.grid().valid())
  {
    std::cout << "CPMAverageCorrectionReconstruction::add - invalid input grid" << std::endl;
    return false;
  }

  if (!m_records.grid().valid())
  {
    m_records.set_grid(
        source.grid().phi_bins,
        source.grid().r_bins,
        source.grid().z_bins,
        source.grid().r_min,
        source.grid().r_max,
        source.grid().z_min,
        source.grid().z_max,
        source.grid().phi_min,
        source.grid().phi_max);
  }
  else if (!CPMReconstructionHelper::same_grid(m_records.grid(), source.grid()))
  {
    std::cout << "CPMAverageCorrectionReconstruction::add - inconsistent input grid" << std::endl;
    return false;
  }

  const auto source_records = source.record_count();
  if (!m_records.add(source))
  {
    return false;
  }

  m_summary.input_records += source_records;
  for (const auto& [voxel, records] : source.records())
  {
    (void) records;
    m_input_voxels.insert(voxel);
  }
  m_summary.input_voxels = static_cast<unsigned int>(m_input_voxels.size());
  return true;
}

bool CPMAverageCorrectionReconstruction::add_from_file(
    const std::string& filename,
    const std::string& objectname)
{
  std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    std::cout << "CPMAverageCorrectionReconstruction::add_from_file - could not open "
              << filename << std::endl;
    return false;
  }

  auto* source = dynamic_cast<CPMVoxelContainer*>(input->Get(objectname.c_str()));
  if (!source)
  {
    std::cout << "CPMAverageCorrectionReconstruction::add_from_file - could not find object "
              << objectname << " in " << filename << std::endl;
    return false;
  }

  const bool ok = add(*source);
  if (ok)
  {
    ++m_summary.input_files;
  }
  return ok;
}

bool CPMAverageCorrectionReconstruction::set_crossing_solver(const std::string& solver)
{
  bool ok = false;
  const auto parsed = CPMReconstructionHelper::solver_from_string(solver, ok);
  if (!ok)
  {
    std::cout << "CPMAverageCorrectionReconstruction::set_crossing_solver - invalid solver "
              << solver << " (expected helix or line)" << std::endl;
    return false;
  }
  set_crossing_solver(parsed);
  return true;
}

void CPMAverageCorrectionReconstruction::reset_output()
{
  m_hentries_rec.reset();
  m_hdistortion_r_rec.reset();
  m_hdistortion_p_rec.reset();
  m_hdistortion_z_rec.reset();
}

void CPMAverageCorrectionReconstruction::reset_calculation_summary()
{
  const unsigned int input_files = m_summary.input_files;
  const unsigned long long input_records = m_summary.input_records;
  const unsigned int input_voxels = m_summary.input_voxels;
  m_summary = {};
  m_summary.input_files = input_files;
  m_summary.input_records = input_records;
  m_summary.input_voxels = input_voxels;
}

bool CPMAverageCorrectionReconstruction::process_loaded_records()
{
  const auto& grid = m_records.grid();
  if (!grid.valid())
  {
    std::cout << "CPMAverageCorrectionReconstruction::process_loaded_records - no valid input grid" << std::endl;
    return false;
  }

  if (m_records.record_count() == 0)
  {
    return true;
  }

  std::vector<VoxelId> voxels;
  voxels.reserve(m_records.voxel_count());
  for (const auto& [voxel, records] : m_records.records())
  {
    (void) records;
    voxels.push_back(voxel);
  }

  std::size_t processed_voxels_in_chunk = 0;
  for (const auto& voxel : voxels)
  {
    const auto records = m_records.take_records(voxel);
    ++processed_voxels_in_chunk;
    if (m_max_records_per_voxel > 0 && records.size() > m_max_records_per_voxel)
    {
      ++m_summary.skipped_large_voxels;
      if (processed_voxels_in_chunk % 256 == 0)
      {
        trim_heap_after_chunk();
      }
      continue;
    }

    std::map<EventTrackKeyView, const TrackStateRecord*> closest_record_by_track;
    unsigned long long good_records = 0;
    for (const auto& record : records)
    {
      if (!CPMReconstructionHelper::pair_has_good_pt(
              record.track.charge,
              record.track.pt,
              m_min_pair_pt))
      {
        ++m_summary.pt_rejected_records;
        continue;
      }
      ++good_records;
      auto [iter, inserted] = closest_record_by_track.emplace(
          event_track_key_view(record),
          &record);
      if (!inserted &&
          CPMReconstructionHelper::is_closer_to_voxel_center(record, *iter->second))
      {
        iter->second = &record;
      }
    }
    m_summary.duplicate_dropped_records +=
        good_records - static_cast<unsigned long long>(closest_record_by_track.size());

    std::vector<const TrackStateRecord*> selected_records;
    selected_records.reserve(closest_record_by_track.size());
    for (const auto& [event_track, record] : closest_record_by_track)
    {
      (void) event_track;
      selected_records.push_back(record);
    }

    if (selected_records.size() < 2)
    {
      ++m_summary.skipped_low_selected_voxels;
      if (processed_voxels_in_chunk % 256 == 0)
      {
        trim_heap_after_chunk();
      }
      continue;
    }

    std::stable_sort(
        selected_records.begin(),
        selected_records.end(),
        [](const TrackStateRecord* lhs, const TrackStateRecord* rhs)
        {
          const auto lhs_hash = CPMReconstructionHelper::record_selection_hash(*lhs);
          const auto rhs_hash = CPMReconstructionHelper::record_selection_hash(*rhs);
          if (lhs_hash != rhs_hash)
          {
            return lhs_hash < rhs_hash;
          }
          if (lhs->track_ref.track_id != rhs->track_ref.track_id)
          {
            return lhs->track_ref.track_id < rhs->track_ref.track_id;
          }
          return lhs->cluster_ref.cluskey < rhs->cluster_ref.cluskey;
        });

    std::vector<const TrackStateRecord*> positive_records;
    std::vector<const TrackStateRecord*> negative_records;
    positive_records.reserve(selected_records.size());
    negative_records.reserve(selected_records.size());
    for (const auto* record : selected_records)
    {
      if (record->track.charge > 0)
      {
        positive_records.push_back(record);
      }
      else if (record->track.charge < 0)
      {
        negative_records.push_back(record);
      }
    }

    if (positive_records.size() < m_min_records_per_charge ||
        negative_records.size() < m_min_records_per_charge)
    {
      ++m_summary.skipped_low_charge_voxels;
      if (processed_voxels_in_chunk % 256 == 0)
      {
        trim_heap_after_chunk();
      }
      continue;
    }

    const std::size_t batch_charge_limit =
        m_max_pair_records_per_charge_batch > 0 ?
        static_cast<std::size_t>(m_max_pair_records_per_charge_batch) :
        std::numeric_limits<std::size_t>::max();

    const bool had_accumulator = m_accumulators.find(voxel) != m_accumulators.end();
    bool accepted_any_batch = false;
    std::size_t positive_index = 0;
    std::size_t negative_index = 0;
    while (positive_index < positive_records.size() &&
           negative_index < negative_records.size())
    {
      const std::size_t positive_take = std::min(
          batch_charge_limit,
          positive_records.size() - positive_index);
      const std::size_t negative_take = std::min(
          batch_charge_limit,
          negative_records.size() - negative_index);

      if (positive_take < m_min_records_per_charge ||
          negative_take < m_min_records_per_charge)
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
          const auto& a = *positive_records[positive_begin + ipos];
          const auto& b = *negative_records[negative_begin + ineg];
          if (same_event_track(a, b))
          {
            ++m_summary.same_event_track_pairs;
            continue;
          }

          const auto pair_result = CPMReconstructionHelper::compute_pair(
              CPMReconstructionHelper::make_pair_input(a),
              CPMReconstructionHelper::make_pair_input(b),
              a.cluster.voxel_center,
              m_pair_options);

          if (pair_result.status == CPMReconstructionHelper::PairStatus::PtRejected)
          {
            continue;
          }

          ++m_summary.candidate_pairs;
          if (pair_result.status == CPMReconstructionHelper::PairStatus::InvalidWeight)
          {
            ++m_summary.invalid_weight_pairs;
            continue;
          }
          if (pair_result.status == CPMReconstructionHelper::PairStatus::InvalidPoCA)
          {
            ++m_summary.invalid_poca_pairs;
            continue;
          }
          if (pair_result.status == CPMReconstructionHelper::PairStatus::DcaRejected)
          {
            ++m_summary.dca_rejected_pairs;
            continue;
          }
          if (!pair_result.accepted())
          {
            continue;
          }

          const double accumulation_weight =
              m_use_pair_weights ? pair_result.pair_weight : 1.0;
          m_accumulators[voxel].add(
              pair_result.delta_r,
              pair_result.delta_rphi,
              pair_result.delta_phi,
              pair_result.delta_z,
              pair_result.dca,
              accumulation_weight,
              a.cluster.voxel_center.X(),
              a.cluster.voxel_center.Y(),
              a.cluster.voxel_center.Z());

          ++m_summary.accepted_pairs;
          accepted_any_batch = true;
        }
      }
    }

    if (accepted_any_batch && !had_accumulator)
    {
      ++m_summary.processed_voxels;
    }
    else
    {
      ++m_summary.skipped_low_batch_charge_voxels;
    }

    if (processed_voxels_in_chunk % 256 == 0)
    {
      trim_heap_after_chunk();
    }
  }

  m_summary.accumulator_voxels = m_accumulators.size();
  m_records.Reset();
  trim_heap_after_chunk();
  return true;
}

bool CPMAverageCorrectionReconstruction::finalize_average_corrections()
{
  const auto& grid = m_records.grid();
  if (!grid.valid())
  {
    std::cout << "CPMAverageCorrectionReconstruction::finalize_average_corrections - no valid input grid" << std::endl;
    return false;
  }

  reset_output();
  m_summary.accumulator_voxels = m_accumulators.size();
  m_summary.filled_voxels = 0;
  m_summary.skipped_low_entry_voxels = 0;
  m_summary.skipped_invalid_voxels = 0;

  m_hentries_rec.reset(new TH3F(
      "hentries_rec",
      "CPM voxel entries;phi;r;z",
      grid.phi_bins, grid.phi_min, grid.phi_max,
      grid.r_bins, grid.r_min, grid.r_max,
      grid.z_bins, grid.z_min, grid.z_max));
  m_hdistortion_r_rec.reset(new TH3F(
      "hDistortionR_rec",
      "CPM radial distortion;phi;r;z",
      grid.phi_bins, grid.phi_min, grid.phi_max,
      grid.r_bins, grid.r_min, grid.r_max,
      grid.z_bins, grid.z_min, grid.z_max));
  m_hdistortion_p_rec.reset(new TH3F(
      "hDistortionP_rec",
      "CPM phi distortion;phi;r;z",
      grid.phi_bins, grid.phi_min, grid.phi_max,
      grid.r_bins, grid.r_min, grid.r_max,
      grid.z_bins, grid.z_min, grid.z_max));
  m_hdistortion_z_rec.reset(new TH3F(
      "hDistortionZ_rec",
      "CPM z distortion;phi;r;z",
      grid.phi_bins, grid.phi_min, grid.phi_max,
      grid.r_bins, grid.r_min, grid.r_max,
      grid.z_bins, grid.z_min, grid.z_max));

  for (const auto& [voxel, accumulator] : m_accumulators)
  {
    if (accumulator.entries < m_min_entries_per_voxel)
    {
      ++m_summary.skipped_low_entry_voxels;
      continue;
    }
    if (voxel.iphi < 0 || voxel.iphi >= grid.phi_bins ||
        voxel.ir < 0 || voxel.ir >= grid.r_bins ||
        voxel.iz < 0 || voxel.iz >= grid.z_bins)
    {
      ++m_summary.skipped_invalid_voxels;
      continue;
    }

    const double mean_delta_r = CPMReconstructionHelper::correction_weighted_mean(
        accumulator.sum_weighted_delta_r,
        accumulator.sum_pair_weight);
    const double mean_delta_phi = CPMReconstructionHelper::correction_weighted_mean(
        accumulator.sum_weighted_delta_phi,
        accumulator.sum_pair_weight);
    const double mean_delta_z = CPMReconstructionHelper::correction_weighted_mean(
        accumulator.sum_weighted_delta_z,
        accumulator.sum_pair_weight);

    if (!std::isfinite(mean_delta_r) ||
        !std::isfinite(mean_delta_phi) ||
        !std::isfinite(mean_delta_z))
    {
      ++m_summary.skipped_invalid_voxels;
      continue;
    }

    m_hentries_rec->SetBinContent(
        voxel.iphi + 1,
        voxel.ir + 1,
        voxel.iz + 1,
        static_cast<double>(accumulator.entries));
    m_hdistortion_r_rec->SetBinContent(voxel.iphi + 1, voxel.ir + 1, voxel.iz + 1, mean_delta_r);
    m_hdistortion_p_rec->SetBinContent(voxel.iphi + 1, voxel.ir + 1, voxel.iz + 1, mean_delta_phi);
    m_hdistortion_z_rec->SetBinContent(voxel.iphi + 1, voxel.ir + 1, voxel.iz + 1, mean_delta_z);
    ++m_summary.filled_voxels;
  }

  return true;
}

bool CPMAverageCorrectionReconstruction::calculate_average_corrections()
{
  reset_output();
  reset_calculation_summary();
  m_accumulators.clear();
  if (!process_loaded_records())
  {
    return false;
  }
  return finalize_average_corrections();
}

bool CPMAverageCorrectionReconstruction::save_average_corrections(
    const std::string& filename) const
{
  if (!m_hentries_rec ||
      !m_hdistortion_r_rec ||
      !m_hdistortion_p_rec ||
      !m_hdistortion_z_rec)
  {
    std::cout << "CPMAverageCorrectionReconstruction::save_average_corrections - no histograms calculated" << std::endl;
    return false;
  }

  TFile output(filename.c_str(), "RECREATE");
  if (output.IsZombie())
  {
    std::cout << "CPMAverageCorrectionReconstruction::save_average_corrections - could not create "
              << filename << std::endl;
    return false;
  }

  m_hentries_rec->Write();
  m_hdistortion_r_rec->Write();
  m_hdistortion_p_rec->Write();
  m_hdistortion_z_rec->Write();

  CPMReconstructionHelper::write_guarded_histograms(output, m_hentries_rec.get(), "hentries");
  CPMReconstructionHelper::write_guarded_histograms(output, m_hdistortion_r_rec.get(), "hIntDistortionR");
  CPMReconstructionHelper::write_guarded_histograms(output, m_hdistortion_p_rec.get(), "hIntDistortionP");
  CPMReconstructionHelper::write_guarded_histograms(output, m_hdistortion_z_rec.get(), "hIntDistortionZ");

  write_summary_tree(output);
  output.Close();
  return true;
}

void CPMAverageCorrectionReconstruction::write_summary_tree(TFile& output) const
{
  output.cd();
  TTree summary_tree("cpm_b3_summary", "CPM production average-correction summary");
  unsigned int input_files = m_summary.input_files;
  unsigned long long input_records = m_summary.input_records;
  unsigned int input_voxels = m_summary.input_voxels;
  unsigned int processed_voxels = m_summary.processed_voxels;
  unsigned int skipped_large_voxels = m_summary.skipped_large_voxels;
  unsigned int skipped_low_selected_voxels = m_summary.skipped_low_selected_voxels;
  unsigned int skipped_low_charge_voxels = m_summary.skipped_low_charge_voxels;
  unsigned int skipped_low_batch_charge_voxels = m_summary.skipped_low_batch_charge_voxels;
  unsigned long long pt_rejected_records = m_summary.pt_rejected_records;
  unsigned long long duplicate_dropped_records = m_summary.duplicate_dropped_records;
  unsigned long long candidate_pairs = m_summary.candidate_pairs;
  unsigned long long accepted_pairs = m_summary.accepted_pairs;
  unsigned long long same_event_track_pairs = m_summary.same_event_track_pairs;
  unsigned long long invalid_weight_pairs = m_summary.invalid_weight_pairs;
  unsigned long long invalid_poca_pairs = m_summary.invalid_poca_pairs;
  unsigned long long dca_rejected_pairs = m_summary.dca_rejected_pairs;
  unsigned int accumulator_voxels = m_summary.accumulator_voxels;
  unsigned int filled_voxels = m_summary.filled_voxels;
  unsigned int skipped_low_entry_voxels = m_summary.skipped_low_entry_voxels;
  unsigned int skipped_invalid_voxels = m_summary.skipped_invalid_voxels;
  unsigned int min_entries_per_voxel = m_min_entries_per_voxel;
  double max_pair_dca = m_pair_options.max_pair_dca;
  double min_sin_angle = m_pair_options.min_sin_angle;
  unsigned int max_records_per_voxel = m_max_records_per_voxel;
  unsigned int min_records_per_charge = m_min_records_per_charge;
  double min_pair_pt = m_min_pair_pt;
  unsigned int max_pair_records_per_charge_batch = m_max_pair_records_per_charge_batch;
  bool use_pair_weights = m_use_pair_weights;
  std::string averaging_mode = m_use_pair_weights ? "weighted" : "plain";
  std::string crossing_solver = CPMReconstructionHelper::solver_name(m_pair_options.solver);
  double magnetic_field_z = m_pair_options.magnetic_field_z;
  int phi_bins = m_records.grid().phi_bins;
  int r_bins = m_records.grid().r_bins;
  int z_bins = m_records.grid().z_bins;
  double phi_min = m_records.grid().phi_min;
  double phi_max = m_records.grid().phi_max;
  double r_min = m_records.grid().r_min;
  double r_max = m_records.grid().r_max;
  double z_min = m_records.grid().z_min;
  double z_max = m_records.grid().z_max;

  summary_tree.Branch("input_files", &input_files);
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
  summary_tree.Branch("min_entries_per_voxel", &min_entries_per_voxel);
  summary_tree.Branch("max_pair_dca", &max_pair_dca);
  summary_tree.Branch("min_sin_angle", &min_sin_angle);
  summary_tree.Branch("max_records_per_voxel", &max_records_per_voxel);
  summary_tree.Branch("min_records_per_charge", &min_records_per_charge);
  summary_tree.Branch("min_pair_pt", &min_pair_pt);
  summary_tree.Branch("max_pair_records_per_charge_batch", &max_pair_records_per_charge_batch);
  summary_tree.Branch("use_pair_weights", &use_pair_weights);
  summary_tree.Branch("averaging_mode", &averaging_mode);
  summary_tree.Branch("crossing_solver", &crossing_solver);
  summary_tree.Branch("magnetic_field_z", &magnetic_field_z);
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
}

void CPMAverageCorrectionReconstruction::print_summary(std::ostream& out) const
{
  out << "CPMAverageCorrectionReconstruction - input files: "
      << m_summary.input_files << std::endl;
  out << "CPMAverageCorrectionReconstruction - input records: "
      << m_summary.input_records << std::endl;
  out << "CPMAverageCorrectionReconstruction - input voxels: "
      << m_summary.input_voxels << std::endl;
  out << "CPMAverageCorrectionReconstruction - processed voxels: "
      << m_summary.processed_voxels << std::endl;
  out << "CPMAverageCorrectionReconstruction - candidate pairs: "
      << m_summary.candidate_pairs << std::endl;
  out << "CPMAverageCorrectionReconstruction - accepted pairs: "
      << m_summary.accepted_pairs << std::endl;
  out << "CPMAverageCorrectionReconstruction - filled voxels: "
      << m_summary.filled_voxels << std::endl;
  out << "CPMAverageCorrectionReconstruction - averaging mode: "
      << (m_use_pair_weights ? "weighted" : "plain") << std::endl;
  out << "CPMAverageCorrectionReconstruction - crossing solver: "
      << CPMReconstructionHelper::solver_name(m_pair_options.solver) << std::endl;
}
