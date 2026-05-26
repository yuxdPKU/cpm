/*
 * CPM QA Job B2 voxel-level accumulator.
 *
 * This macro reads B1 crossing-point PoCA outputs, groups accepted pairs
 * by 3D voxel, and writes one correction-summary row per voxel. It can read
 * either pair-level rows or B1 batch-level running sums. Pairwise crossing
 * estimates can be averaged either with the B1 curvature proxy weight
 * 1/(pt_i pt_j), proportional to (1/R_i)(1/R_j) in a fixed magnetic field, or
 * with unit weights for a simple arithmetic mean.
 * The output is a QA/intermediate product. The input B1 delta convention is the
 * same as the TpcDistortionCorrection distortion convention:
 * voxel center - crossing point.
 */

#include <CPMRecord.h>
#include <CPMReconstructionHelper.h>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace CPMB2
{
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

  bool has_tree_with_entries(
      const std::vector<std::string>& input_files,
      const std::string& tree_name)
  {
    for (const auto& file : input_files)
    {
      TFile input(file.c_str(), "READ");
      if (input.IsZombie())
      {
        continue;
      }
      auto* tree = dynamic_cast<TTree*>(input.Get(tree_name.c_str()));
      if (tree && tree->GetEntries() > 0)
      {
        return true;
      }
    }
    return false;
  }
}

void CPM_QA_B2_AccumulateVoxelCorrections(
    const std::vector<std::string>& input_files,
    const std::string& output_file = "CPM_QA_B2_voxel_corrections.root",
    const unsigned int min_entries_per_voxel = 1,
    const double max_pair_dca = -1.0,
    const bool use_pair_weights = true,
    const std::string& input_mode = "auto")
{
  const bool has_batch_tree =
      CPMB2::has_tree_with_entries(input_files, "cpm_b1_batch_corrections");
  std::string resolved_input_mode = input_mode;
  if (resolved_input_mode == "auto")
  {
    resolved_input_mode = has_batch_tree ? "batches" : "pairs";
  }
  if (resolved_input_mode != "pairs" && resolved_input_mode != "batches")
  {
    std::cerr << "CPM_QA_B2_AccumulateVoxelCorrections - invalid input_mode: "
              << input_mode << " (expected auto, pairs, or batches)" << std::endl;
    return;
  }
  if (resolved_input_mode == "batches" && !has_batch_tree)
  {
    std::cerr << "CPM_QA_B2_AccumulateVoxelCorrections - requested batch input, "
              << "but no cpm_b1_batch_corrections tree with entries was found"
              << std::endl;
    return;
  }
  if (resolved_input_mode == "batches" && max_pair_dca >= 0.0)
  {
    std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - warning: max_pair_dca "
              << "is already applied in B1 batch sums; B2 will not refilter "
              << "individual pairs in batch input mode" << std::endl;
  }

  std::map<VoxelId, CPMReconstructionHelper::CorrectionAccumulator> accumulators;
  unsigned long long accepted_pairs = 0;
  unsigned long long rejected_rows = 0;
  Long64_t input_rows = 0;
  Long64_t input_pair_rows = 0;
  Long64_t input_batch_rows = 0;
  bool has_pair_weight = false;

  if (resolved_input_mode == "pairs")
  {
    TChain chain("cpm_poca_pairs");
    for (const auto& file : input_files)
    {
      chain.Add(file.c_str());
    }

    int iphi = -1;
    int ir = -1;
    int iz = -1;
    double dca = std::numeric_limits<double>::quiet_NaN();
    double voxel_center_x = std::numeric_limits<double>::quiet_NaN();
    double voxel_center_y = std::numeric_limits<double>::quiet_NaN();
    double voxel_center_z = std::numeric_limits<double>::quiet_NaN();
    double pair_weight = 1.0;
    double delta_r = std::numeric_limits<double>::quiet_NaN();
    double delta_rphi = std::numeric_limits<double>::quiet_NaN();
    double delta_phi = std::numeric_limits<double>::quiet_NaN();
    double delta_z = std::numeric_limits<double>::quiet_NaN();

    chain.SetBranchAddress("iphi", &iphi);
    chain.SetBranchAddress("ir", &ir);
    chain.SetBranchAddress("iz", &iz);
    chain.SetBranchAddress("dca", &dca);
    chain.SetBranchAddress("voxel_center_x", &voxel_center_x);
    chain.SetBranchAddress("voxel_center_y", &voxel_center_y);
    chain.SetBranchAddress("voxel_center_z", &voxel_center_z);
    has_pair_weight = chain.GetBranch("pair_weight") != nullptr;
    if (has_pair_weight)
    {
      chain.SetBranchAddress("pair_weight", &pair_weight);
    }
    chain.SetBranchAddress("delta_r", &delta_r);
    chain.SetBranchAddress("delta_rphi", &delta_rphi);
    chain.SetBranchAddress("delta_phi", &delta_phi);
    chain.SetBranchAddress("delta_z", &delta_z);

    input_rows = chain.GetEntries();
    input_pair_rows = input_rows;
    for (Long64_t entry = 0; entry < input_rows; ++entry)
    {
      chain.GetEntry(entry);
      if (!has_pair_weight)
      {
        pair_weight = 1.0;
      }

      if (iphi < 0 || ir < 0 || iz < 0 ||
          !std::isfinite(delta_r) || !std::isfinite(delta_rphi) ||
          !std::isfinite(delta_phi) || !std::isfinite(delta_z) ||
          !std::isfinite(dca) ||
          (use_pair_weights && (!std::isfinite(pair_weight) || pair_weight <= 0.0)))
      {
        ++rejected_rows;
        continue;
      }
      if (max_pair_dca >= 0.0 && dca > max_pair_dca)
      {
        ++rejected_rows;
        continue;
      }

      const double accumulation_weight = use_pair_weights ? pair_weight : 1.0;
      accumulators[{iphi, ir, iz}].add(
          delta_r,
          delta_rphi,
          delta_phi,
          delta_z,
          dca,
          accumulation_weight,
          voxel_center_x,
          voxel_center_y,
          voxel_center_z);
      ++accepted_pairs;
    }
  }
  else
  {
    TChain chain("cpm_b1_batch_corrections");
    for (const auto& file : input_files)
    {
      chain.Add(file.c_str());
    }

    int iphi = -1;
    int ir = -1;
    int iz = -1;
    unsigned long long batch_entries = 0;
    double voxel_x = std::numeric_limits<double>::quiet_NaN();
    double voxel_y = std::numeric_limits<double>::quiet_NaN();
    double voxel_z = std::numeric_limits<double>::quiet_NaN();
    double sum_pair_weight = std::numeric_limits<double>::quiet_NaN();
    double sum_pair_weight2 = std::numeric_limits<double>::quiet_NaN();
    double sum_delta_r = std::numeric_limits<double>::quiet_NaN();
    double sum_delta_r2 = std::numeric_limits<double>::quiet_NaN();
    double sum_delta_rphi = std::numeric_limits<double>::quiet_NaN();
    double sum_delta_rphi2 = std::numeric_limits<double>::quiet_NaN();
    double sum_delta_phi = std::numeric_limits<double>::quiet_NaN();
    double sum_delta_phi2 = std::numeric_limits<double>::quiet_NaN();
    double sum_delta_z = std::numeric_limits<double>::quiet_NaN();
    double sum_delta_z2 = std::numeric_limits<double>::quiet_NaN();
    double sum_weighted_delta_r = std::numeric_limits<double>::quiet_NaN();
    double sum_weighted_delta_r2 = std::numeric_limits<double>::quiet_NaN();
    double sum_weighted_delta_rphi = std::numeric_limits<double>::quiet_NaN();
    double sum_weighted_delta_rphi2 = std::numeric_limits<double>::quiet_NaN();
    double sum_weighted_delta_phi = std::numeric_limits<double>::quiet_NaN();
    double sum_weighted_delta_phi2 = std::numeric_limits<double>::quiet_NaN();
    double sum_weighted_delta_z = std::numeric_limits<double>::quiet_NaN();
    double sum_weighted_delta_z2 = std::numeric_limits<double>::quiet_NaN();
    double sum_dca = std::numeric_limits<double>::quiet_NaN();
    double sum_dca2 = std::numeric_limits<double>::quiet_NaN();

    chain.SetBranchAddress("iphi", &iphi);
    chain.SetBranchAddress("ir", &ir);
    chain.SetBranchAddress("iz", &iz);
    chain.SetBranchAddress("accepted_pairs", &batch_entries);
    chain.SetBranchAddress("voxel_x", &voxel_x);
    chain.SetBranchAddress("voxel_y", &voxel_y);
    chain.SetBranchAddress("voxel_z", &voxel_z);
    chain.SetBranchAddress("sum_pair_weight", &sum_pair_weight);
    chain.SetBranchAddress("sum_pair_weight2", &sum_pair_weight2);
    chain.SetBranchAddress("sum_delta_r", &sum_delta_r);
    chain.SetBranchAddress("sum_delta_r2", &sum_delta_r2);
    chain.SetBranchAddress("sum_delta_rphi", &sum_delta_rphi);
    chain.SetBranchAddress("sum_delta_rphi2", &sum_delta_rphi2);
    chain.SetBranchAddress("sum_delta_phi", &sum_delta_phi);
    chain.SetBranchAddress("sum_delta_phi2", &sum_delta_phi2);
    chain.SetBranchAddress("sum_delta_z", &sum_delta_z);
    chain.SetBranchAddress("sum_delta_z2", &sum_delta_z2);
    chain.SetBranchAddress("sum_weighted_delta_r", &sum_weighted_delta_r);
    chain.SetBranchAddress("sum_weighted_delta_r2", &sum_weighted_delta_r2);
    chain.SetBranchAddress("sum_weighted_delta_rphi", &sum_weighted_delta_rphi);
    chain.SetBranchAddress("sum_weighted_delta_rphi2", &sum_weighted_delta_rphi2);
    chain.SetBranchAddress("sum_weighted_delta_phi", &sum_weighted_delta_phi);
    chain.SetBranchAddress("sum_weighted_delta_phi2", &sum_weighted_delta_phi2);
    chain.SetBranchAddress("sum_weighted_delta_z", &sum_weighted_delta_z);
    chain.SetBranchAddress("sum_weighted_delta_z2", &sum_weighted_delta_z2);
    chain.SetBranchAddress("sum_dca", &sum_dca);
    chain.SetBranchAddress("sum_dca2", &sum_dca2);

    has_pair_weight = true;
    input_rows = chain.GetEntries();
    input_batch_rows = input_rows;
    for (Long64_t entry = 0; entry < input_rows; ++entry)
    {
      chain.GetEntry(entry);
      if (iphi < 0 || ir < 0 || iz < 0 || batch_entries == 0 ||
          !std::isfinite(voxel_x) || !std::isfinite(voxel_y) ||
          !std::isfinite(voxel_z) ||
          !std::isfinite(sum_dca) || !std::isfinite(sum_dca2))
      {
        ++rejected_rows;
        continue;
      }

      double add_sum_pair_weight = static_cast<double>(batch_entries);
      double add_sum_pair_weight2 = static_cast<double>(batch_entries);
      double add_sum_delta_r = sum_delta_r;
      double add_sum_delta_r2 = sum_delta_r2;
      double add_sum_delta_rphi = sum_delta_rphi;
      double add_sum_delta_rphi2 = sum_delta_rphi2;
      double add_sum_delta_phi = sum_delta_phi;
      double add_sum_delta_phi2 = sum_delta_phi2;
      double add_sum_delta_z = sum_delta_z;
      double add_sum_delta_z2 = sum_delta_z2;

      if (use_pair_weights)
      {
        add_sum_pair_weight = sum_pair_weight;
        add_sum_pair_weight2 = sum_pair_weight2;
        add_sum_delta_r = sum_weighted_delta_r;
        add_sum_delta_r2 = sum_weighted_delta_r2;
        add_sum_delta_rphi = sum_weighted_delta_rphi;
        add_sum_delta_rphi2 = sum_weighted_delta_rphi2;
        add_sum_delta_phi = sum_weighted_delta_phi;
        add_sum_delta_phi2 = sum_weighted_delta_phi2;
        add_sum_delta_z = sum_weighted_delta_z;
        add_sum_delta_z2 = sum_weighted_delta_z2;
      }

      if (!std::isfinite(add_sum_pair_weight) || add_sum_pair_weight <= 0.0 ||
          !std::isfinite(add_sum_pair_weight2) ||
          !std::isfinite(add_sum_delta_r) ||
          !std::isfinite(add_sum_delta_r2) ||
          !std::isfinite(add_sum_delta_rphi) ||
          !std::isfinite(add_sum_delta_rphi2) ||
          !std::isfinite(add_sum_delta_phi) ||
          !std::isfinite(add_sum_delta_phi2) ||
          !std::isfinite(add_sum_delta_z) ||
          !std::isfinite(add_sum_delta_z2))
      {
        ++rejected_rows;
        continue;
      }

      accumulators[{iphi, ir, iz}].add_sums(
          batch_entries,
          add_sum_pair_weight,
          add_sum_pair_weight2,
          add_sum_delta_r,
          add_sum_delta_r2,
          add_sum_delta_rphi,
          add_sum_delta_rphi2,
          add_sum_delta_phi,
          add_sum_delta_phi2,
          add_sum_delta_z,
          add_sum_delta_z2,
          sum_dca,
          sum_dca2,
          voxel_x,
          voxel_y,
          voxel_z);
      accepted_pairs += batch_entries;
    }
  }

  TFile output(output_file.c_str(), "RECREATE");

  TTree voxels("cpm_voxel_corrections", "CPM voxel-level correction QA");
  int out_iphi = -1;
  int out_ir = -1;
  int out_iz = -1;
  unsigned long long entries = 0;
  double voxel_x = std::numeric_limits<double>::quiet_NaN();
  double voxel_y = std::numeric_limits<double>::quiet_NaN();
  double voxel_z = std::numeric_limits<double>::quiet_NaN();
  double sum_pair_weight = std::numeric_limits<double>::quiet_NaN();
  double mean_pair_weight = std::numeric_limits<double>::quiet_NaN();
  double effective_pair_entries = std::numeric_limits<double>::quiet_NaN();
  double mean_delta_r = std::numeric_limits<double>::quiet_NaN();
  double rms_delta_r = std::numeric_limits<double>::quiet_NaN();
  double mean_delta_rphi = std::numeric_limits<double>::quiet_NaN();
  double rms_delta_rphi = std::numeric_limits<double>::quiet_NaN();
  double mean_delta_phi = std::numeric_limits<double>::quiet_NaN();
  double rms_delta_phi = std::numeric_limits<double>::quiet_NaN();
  double mean_delta_z = std::numeric_limits<double>::quiet_NaN();
  double rms_delta_z = std::numeric_limits<double>::quiet_NaN();
  double mean_dca = std::numeric_limits<double>::quiet_NaN();
  double rms_dca = std::numeric_limits<double>::quiet_NaN();

  voxels.Branch("iphi", &out_iphi);
  voxels.Branch("ir", &out_ir);
  voxels.Branch("iz", &out_iz);
  voxels.Branch("entries", &entries);
  voxels.Branch("voxel_x", &voxel_x);
  voxels.Branch("voxel_y", &voxel_y);
  voxels.Branch("voxel_z", &voxel_z);
  voxels.Branch("sum_pair_weight", &sum_pair_weight);
  voxels.Branch("mean_pair_weight", &mean_pair_weight);
  voxels.Branch("effective_pair_entries", &effective_pair_entries);
  voxels.Branch("mean_delta_r", &mean_delta_r);
  voxels.Branch("rms_delta_r", &rms_delta_r);
  voxels.Branch("mean_delta_rphi", &mean_delta_rphi);
  voxels.Branch("rms_delta_rphi", &rms_delta_rphi);
  voxels.Branch("mean_delta_phi", &mean_delta_phi);
  voxels.Branch("rms_delta_phi", &rms_delta_phi);
  voxels.Branch("mean_delta_z", &mean_delta_z);
  voxels.Branch("rms_delta_z", &rms_delta_z);
  voxels.Branch("mean_dca", &mean_dca);
  voxels.Branch("rms_dca", &rms_dca);

  unsigned int filled_voxels = 0;
  unsigned int skipped_low_entry_voxels = 0;
  for (const auto& [voxel, accumulator] : accumulators)
  {
    if (accumulator.entries < min_entries_per_voxel)
    {
      ++skipped_low_entry_voxels;
      continue;
    }

    out_iphi = voxel.iphi;
    out_ir = voxel.ir;
    out_iz = voxel.iz;
    entries = accumulator.entries;
    voxel_x = CPMReconstructionHelper::correction_mean(accumulator.sum_voxel_x, entries);
    voxel_y = CPMReconstructionHelper::correction_mean(accumulator.sum_voxel_y, entries);
    voxel_z = CPMReconstructionHelper::correction_mean(accumulator.sum_voxel_z, entries);
    sum_pair_weight = accumulator.sum_pair_weight;
    mean_pair_weight = CPMReconstructionHelper::correction_mean(accumulator.sum_pair_weight, entries);
    effective_pair_entries = CPMReconstructionHelper::correction_effective_entries(
        accumulator.sum_pair_weight,
        accumulator.sum_pair_weight2);
    mean_delta_r = CPMReconstructionHelper::correction_weighted_mean(
        accumulator.sum_weighted_delta_r,
        accumulator.sum_pair_weight);
    rms_delta_r = CPMReconstructionHelper::correction_weighted_rms(
        accumulator.sum_weighted_delta_r,
        accumulator.sum_weighted_delta_r2,
        accumulator.sum_pair_weight);
    mean_delta_rphi = CPMReconstructionHelper::correction_weighted_mean(
        accumulator.sum_weighted_delta_rphi,
        accumulator.sum_pair_weight);
    rms_delta_rphi = CPMReconstructionHelper::correction_weighted_rms(
        accumulator.sum_weighted_delta_rphi,
        accumulator.sum_weighted_delta_rphi2,
        accumulator.sum_pair_weight);
    mean_delta_phi = CPMReconstructionHelper::correction_weighted_mean(
        accumulator.sum_weighted_delta_phi,
        accumulator.sum_pair_weight);
    rms_delta_phi = CPMReconstructionHelper::correction_weighted_rms(
        accumulator.sum_weighted_delta_phi,
        accumulator.sum_weighted_delta_phi2,
        accumulator.sum_pair_weight);
    mean_delta_z = CPMReconstructionHelper::correction_weighted_mean(
        accumulator.sum_weighted_delta_z,
        accumulator.sum_pair_weight);
    rms_delta_z = CPMReconstructionHelper::correction_weighted_rms(
        accumulator.sum_weighted_delta_z,
        accumulator.sum_weighted_delta_z2,
        accumulator.sum_pair_weight);
    mean_dca = CPMReconstructionHelper::correction_mean(accumulator.sum_dca, entries);
    rms_dca = CPMReconstructionHelper::correction_rms(accumulator.sum_dca, accumulator.sum_dca2, entries);

    voxels.Fill();
    ++filled_voxels;
  }

  TTree summary("cpm_b2_summary", "CPM B2 voxel accumulator summary");
  unsigned int input_files_count = input_files.size();
  unsigned int accumulator_voxels = accumulators.size();
  unsigned int summary_min_entries_per_voxel = min_entries_per_voxel;
  double summary_max_pair_dca = max_pair_dca;
  bool summary_has_pair_weight = has_pair_weight;
  bool summary_use_pair_weights = use_pair_weights;
  std::string summary_averaging_mode = use_pair_weights ? "weighted" : "unweighted";
  std::string summary_requested_input_mode = input_mode;
  std::string summary_input_mode = resolved_input_mode;
  Long64_t summary_input_rows = input_rows;
  Long64_t summary_input_pair_rows = input_pair_rows;
  Long64_t summary_input_batch_rows = input_batch_rows;
  unsigned long long summary_input_pairs =
      resolved_input_mode == "pairs" ?
      static_cast<unsigned long long>(input_pair_rows) :
      accepted_pairs;
  unsigned long long summary_rejected_rows = rejected_rows;

  summary.Branch("input_files", &input_files_count);
  summary.Branch("requested_input_mode", &summary_requested_input_mode);
  summary.Branch("input_mode", &summary_input_mode);
  summary.Branch("input_rows", &summary_input_rows);
  summary.Branch("input_pair_rows", &summary_input_pair_rows);
  summary.Branch("input_batch_rows", &summary_input_batch_rows);
  summary.Branch("input_pairs", &summary_input_pairs);
  summary.Branch("accepted_pairs", &accepted_pairs);
  summary.Branch("rejected_rows", &summary_rejected_rows);
  summary.Branch("accumulator_voxels", &accumulator_voxels);
  summary.Branch("filled_voxels", &filled_voxels);
  summary.Branch("skipped_low_entry_voxels", &skipped_low_entry_voxels);
  summary.Branch("min_entries_per_voxel", &summary_min_entries_per_voxel);
  summary.Branch("max_pair_dca", &summary_max_pair_dca);
  summary.Branch("has_pair_weight", &summary_has_pair_weight);
  summary.Branch("use_pair_weights", &summary_use_pair_weights);
  summary.Branch("averaging_mode", &summary_averaging_mode);
  summary.Fill();

  voxels.Write();
  summary.Write();
  output.Close();

  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - input mode: "
            << resolved_input_mode << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - input rows: " << input_rows << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - input pair rows: " << input_pair_rows << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - input batch rows: " << input_batch_rows << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - input pair estimates: "
            << summary_input_pairs << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - accepted pairs: " << accepted_pairs << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - rejected rows: " << rejected_rows << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - accumulator voxels: " << accumulator_voxels << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - filled voxels: " << filled_voxels << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - skipped low-entry voxels: " << skipped_low_entry_voxels << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - averaging mode: "
            << summary_averaging_mode << std::endl;
  std::cout << "CPM_QA_B2_AccumulateVoxelCorrections - output: " << output_file << std::endl;
}

void CPM_QA_B2_AccumulateVoxelCorrections(
    const std::string& input_file,
    const std::string& output_file = "CPM_QA_B2_voxel_corrections.root",
    const unsigned int min_entries_per_voxel = 1,
    const double max_pair_dca = -1.0,
    const bool use_pair_weights = true,
    const std::string& input_mode = "auto")
{
  CPM_QA_B2_AccumulateVoxelCorrections(
      std::vector<std::string>{input_file},
      output_file,
      min_entries_per_voxel,
      max_pair_dca,
      use_pair_weights,
      input_mode);
}

void CPM_QA_B2_AccumulateVoxelCorrectionsFromList(
    const std::string& input_file_or_list,
    const std::string& output_file,
    const unsigned int min_entries_per_voxel = 1,
    const double max_pair_dca = -1.0,
    const bool use_pair_weights = true,
    const std::string& input_mode = "auto")
{
  const auto input_files = CPMB2::read_file_list(input_file_or_list);

  CPM_QA_B2_AccumulateVoxelCorrections(
      input_files,
      output_file,
      min_entries_per_voxel,
      max_pair_dca,
      use_pair_weights,
      input_mode);
}
