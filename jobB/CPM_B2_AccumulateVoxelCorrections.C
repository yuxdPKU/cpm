/*
 * CPM Job B2 voxel-level accumulator.
 *
 * This macro reads B1 local line-line PoCA pair outputs, groups accepted pairs
 * by 3D voxel, and writes one correction-summary row per voxel. Pairwise
 * crossing estimates are averaged with the B1 curvature proxy weight
 * 1/(pt_i pt_j), proportional to (1/R_i)(1/R_j) in a fixed magnetic field.
 * The output is a QA/intermediate product. The input B1 delta convention is the
 * same as the TpcDistortionCorrection distortion convention:
 * voxel center - crossing point.
 */

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
#include <tuple>
#include <vector>

namespace CPMB2
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

  struct Accumulator
  {
    unsigned long long entries = 0;
    double sum_weight = 0.0;
    double sum_weight2 = 0.0;
    double sum_weighted_delta_r = 0.0;
    double sum_weighted_delta_r2 = 0.0;
    double sum_weighted_delta_rphi = 0.0;
    double sum_weighted_delta_rphi2 = 0.0;
    double sum_weighted_delta_phi = 0.0;
    double sum_weighted_delta_phi2 = 0.0;
    double sum_weighted_delta_z = 0.0;
    double sum_weighted_delta_z2 = 0.0;
    double sum_dca = 0.0;
    double sum_dca2 = 0.0;
    double sum_voxel_x = 0.0;
    double sum_voxel_y = 0.0;
    double sum_voxel_z = 0.0;

    void add(
        const double delta_r,
        const double delta_rphi,
        const double delta_phi,
        const double delta_z,
        const double dca,
        const double weight,
        const double voxel_x,
        const double voxel_y,
        const double voxel_z)
    {
      ++entries;
      sum_weight += weight;
      sum_weight2 += weight * weight;
      sum_weighted_delta_r += weight * delta_r;
      sum_weighted_delta_r2 += weight * delta_r * delta_r;
      sum_weighted_delta_rphi += weight * delta_rphi;
      sum_weighted_delta_rphi2 += weight * delta_rphi * delta_rphi;
      sum_weighted_delta_phi += weight * delta_phi;
      sum_weighted_delta_phi2 += weight * delta_phi * delta_phi;
      sum_weighted_delta_z += weight * delta_z;
      sum_weighted_delta_z2 += weight * delta_z * delta_z;
      sum_dca += dca;
      sum_dca2 += dca * dca;
      sum_voxel_x += voxel_x;
      sum_voxel_y += voxel_y;
      sum_voxel_z += voxel_z;
    }
  };

  double mean(const double sum, const unsigned long long entries)
  {
    return entries > 0 ? sum / static_cast<double>(entries) : std::numeric_limits<double>::quiet_NaN();
  }

  double weighted_mean(const double sum_weighted, const double sum_weight)
  {
    return sum_weight > 0.0 ? sum_weighted / sum_weight : std::numeric_limits<double>::quiet_NaN();
  }

  double weighted_rms(const double sum_weighted, const double sum_weighted2, const double sum_weight)
  {
    if (sum_weight <= 0.0)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double average = weighted_mean(sum_weighted, sum_weight);
    const double variance = sum_weighted2 / sum_weight - average * average;
    return std::sqrt(std::max(0.0, variance));
  }

  double rms(const double sum, const double sum2, const unsigned long long entries)
  {
    if (entries == 0)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double average = mean(sum, entries);
    const double variance = sum2 / static_cast<double>(entries) - average * average;
    return std::sqrt(std::max(0.0, variance));
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

void CPM_B2_AccumulateVoxelCorrections(
    const std::vector<std::string>& input_files,
    const std::string& output_file = "CPM_B2_voxel_corrections.root",
    const unsigned int min_entries_per_voxel = 1,
    const double max_pair_dca = -1.0)
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
  const bool has_pair_weight = chain.GetBranch("pair_weight") != nullptr;
  if (has_pair_weight)
  {
    chain.SetBranchAddress("pair_weight", &pair_weight);
  }
  chain.SetBranchAddress("delta_r", &delta_r);
  chain.SetBranchAddress("delta_rphi", &delta_rphi);
  chain.SetBranchAddress("delta_phi", &delta_phi);
  chain.SetBranchAddress("delta_z", &delta_z);

  std::map<CPMB2::VoxelKey, CPMB2::Accumulator> accumulators;
  unsigned long long accepted_pairs = 0;
  unsigned long long rejected_pairs = 0;

  Long64_t input_pairs = chain.GetEntries();
  for (Long64_t entry = 0; entry < input_pairs; ++entry)
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
        !std::isfinite(pair_weight) || pair_weight <= 0.0)
    {
      ++rejected_pairs;
      continue;
    }
    if (max_pair_dca >= 0.0 && dca > max_pair_dca)
    {
      ++rejected_pairs;
      continue;
    }

    accumulators[{iphi, ir, iz}].add(
        delta_r,
        delta_rphi,
        delta_phi,
        delta_z,
        dca,
        pair_weight,
        voxel_center_x,
        voxel_center_y,
        voxel_center_z);
    ++accepted_pairs;
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
    voxel_x = CPMB2::mean(accumulator.sum_voxel_x, entries);
    voxel_y = CPMB2::mean(accumulator.sum_voxel_y, entries);
    voxel_z = CPMB2::mean(accumulator.sum_voxel_z, entries);
    sum_pair_weight = accumulator.sum_weight;
    mean_pair_weight = CPMB2::mean(accumulator.sum_weight, entries);
    effective_pair_entries = accumulator.sum_weight2 > 0.0 ?
        accumulator.sum_weight * accumulator.sum_weight / accumulator.sum_weight2 :
        std::numeric_limits<double>::quiet_NaN();
    mean_delta_r = CPMB2::weighted_mean(accumulator.sum_weighted_delta_r, accumulator.sum_weight);
    rms_delta_r = CPMB2::weighted_rms(accumulator.sum_weighted_delta_r, accumulator.sum_weighted_delta_r2, accumulator.sum_weight);
    mean_delta_rphi = CPMB2::weighted_mean(accumulator.sum_weighted_delta_rphi, accumulator.sum_weight);
    rms_delta_rphi = CPMB2::weighted_rms(accumulator.sum_weighted_delta_rphi, accumulator.sum_weighted_delta_rphi2, accumulator.sum_weight);
    mean_delta_phi = CPMB2::weighted_mean(accumulator.sum_weighted_delta_phi, accumulator.sum_weight);
    rms_delta_phi = CPMB2::weighted_rms(accumulator.sum_weighted_delta_phi, accumulator.sum_weighted_delta_phi2, accumulator.sum_weight);
    mean_delta_z = CPMB2::weighted_mean(accumulator.sum_weighted_delta_z, accumulator.sum_weight);
    rms_delta_z = CPMB2::weighted_rms(accumulator.sum_weighted_delta_z, accumulator.sum_weighted_delta_z2, accumulator.sum_weight);
    mean_dca = CPMB2::mean(accumulator.sum_dca, entries);
    rms_dca = CPMB2::rms(accumulator.sum_dca, accumulator.sum_dca2, entries);

    voxels.Fill();
    ++filled_voxels;
  }

  TTree summary("cpm_b2_summary", "CPM B2 voxel accumulator summary");
  unsigned int input_files_count = input_files.size();
  unsigned int accumulator_voxels = accumulators.size();
  unsigned int summary_min_entries_per_voxel = min_entries_per_voxel;
  double summary_max_pair_dca = max_pair_dca;
  bool summary_has_pair_weight = has_pair_weight;

  summary.Branch("input_files", &input_files_count);
  summary.Branch("input_pairs", &input_pairs);
  summary.Branch("accepted_pairs", &accepted_pairs);
  summary.Branch("rejected_pairs", &rejected_pairs);
  summary.Branch("accumulator_voxels", &accumulator_voxels);
  summary.Branch("filled_voxels", &filled_voxels);
  summary.Branch("skipped_low_entry_voxels", &skipped_low_entry_voxels);
  summary.Branch("min_entries_per_voxel", &summary_min_entries_per_voxel);
  summary.Branch("max_pair_dca", &summary_max_pair_dca);
  summary.Branch("has_pair_weight", &summary_has_pair_weight);
  summary.Fill();

  voxels.Write();
  summary.Write();
  output.Close();

  std::cout << "CPM_B2_AccumulateVoxelCorrections - input pairs: " << input_pairs << std::endl;
  std::cout << "CPM_B2_AccumulateVoxelCorrections - accepted pairs: " << accepted_pairs << std::endl;
  std::cout << "CPM_B2_AccumulateVoxelCorrections - rejected pairs: " << rejected_pairs << std::endl;
  std::cout << "CPM_B2_AccumulateVoxelCorrections - accumulator voxels: " << accumulator_voxels << std::endl;
  std::cout << "CPM_B2_AccumulateVoxelCorrections - filled voxels: " << filled_voxels << std::endl;
  std::cout << "CPM_B2_AccumulateVoxelCorrections - skipped low-entry voxels: " << skipped_low_entry_voxels << std::endl;
  std::cout << "CPM_B2_AccumulateVoxelCorrections - output: " << output_file << std::endl;
}

void CPM_B2_AccumulateVoxelCorrections(
    const std::string& input_file,
    const std::string& output_file = "CPM_B2_voxel_corrections.root",
    const unsigned int min_entries_per_voxel = 1,
    const double max_pair_dca = -1.0)
{
  CPM_B2_AccumulateVoxelCorrections(
      std::vector<std::string>{input_file},
      output_file,
      min_entries_per_voxel,
      max_pair_dca);
}

void CPM_B2_AccumulateVoxelCorrections(
    const std::string& input_file_or_list,
    const std::string& output_file,
    const bool input_is_list,
    const unsigned int min_entries_per_voxel = 1,
    const double max_pair_dca = -1.0)
{
  const auto input_files = input_is_list ?
      CPMB2::read_file_list(input_file_or_list) :
      std::vector<std::string>{input_file_or_list};

  CPM_B2_AccumulateVoxelCorrections(
      input_files,
      output_file,
      min_entries_per_voxel,
      max_pair_dca);
}
