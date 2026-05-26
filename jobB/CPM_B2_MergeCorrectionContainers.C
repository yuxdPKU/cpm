/*
 * CPM Job B2 container merger.
 *
 * Legacy CPM container path. Older Job A files may contain one
 * CPMCorrectionContainer per segment; this macro merges those containers and
 * writes a cpm_voxel_corrections tree for the existing B3 histogram writer.
 * New CPM production uses the record-based B1/B2/B3 path from cpm_records.
 */

#include <CPMCorrectionContainer.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libcpm.so)

namespace CPMB2Container
{
  double mean(const double sum, const int entries)
  {
    return entries > 0 ? sum / static_cast<double>(entries) : std::numeric_limits<double>::quiet_NaN();
  }

  double rms(const double sum, const double sum2, const int entries)
  {
    if (entries <= 0)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double average = mean(sum, entries);
    const double variance = sum2 / static_cast<double>(entries) - average * average;
    return std::sqrt(std::max(0.0, variance));
  }

  double weighted_rms(const double sumWeighted, const double sumWeighted2, const double sumWeight)
  {
    if (sumWeight <= 0.0)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double average = sumWeighted / sumWeight;
    const double variance = sumWeighted2 / sumWeight - average * average;
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

void CPM_B2_MergeCorrectionContainers(
    const std::vector<std::string>& input_files,
    const std::string& output_file = "CPM_B2_voxel_corrections.root",
    const std::string& object_name = "CPMCorrectionContainer",
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1)
{
  CPMCorrectionContainer merged;
  bool initialized = false;
  unsigned int loaded_files = 0;
  unsigned int rejected_files = 0;

  for (const auto& filename : input_files)
  {
    std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
    if (!input || input->IsZombie())
    {
      std::cout << "CPM_B2_MergeCorrectionContainers - could not open: " << filename << std::endl;
      ++rejected_files;
      continue;
    }

    auto* source = dynamic_cast<CPMCorrectionContainer*>(input->Get(object_name.c_str()));
    if (!source)
    {
      std::cout << "CPM_B2_MergeCorrectionContainers - missing object " << object_name
                << " in " << filename << std::endl;
      ++rejected_files;
      continue;
    }

    if (!initialized)
    {
      int phiBins = 0;
      int rBins = 0;
      int zBins = 0;
      double rMin = 0.0;
      double rMax = 0.0;
      double zMin = 0.0;
      double zMax = 0.0;
      source->get_grid_dimensions(phiBins, rBins, zBins);
      source->get_grid_range(rMin, rMax, zMin, zMax);
      merged.set_grid_dimensions(phiBins, rBins, zBins);
      merged.set_grid_range(rMin, rMax, zMin, zMax);
      initialized = true;
    }

    if (!merged.add(*source))
    {
      ++rejected_files;
      continue;
    }
    ++loaded_files;
  }

  if (!initialized || loaded_files == 0)
  {
    std::cout << "CPM_B2_MergeCorrectionContainers - no valid containers loaded" << std::endl;
    return;
  }

  int phiBins = 0;
  int rBins = 0;
  int zBins = 0;
  double rMin = 0.0;
  double rMax = 0.0;
  double zMin = 0.0;
  double zMax = 0.0;
  merged.get_grid_dimensions(phiBins, rBins, zBins);
  merged.get_grid_range(rMin, rMax, zMin, zMax);

  TFile output(output_file.c_str(), "RECREATE");
  if (output.IsZombie())
  {
    std::cout << "CPM_B2_MergeCorrectionContainers - could not create: " << output_file << std::endl;
    return;
  }

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
  constexpr double pi = 3.14159265358979323846;
  for (int iphi = 0; iphi < phiBins; ++iphi)
  {
    for (int ir = 0; ir < rBins; ++ir)
    {
      for (int iz = 0; iz < zBins; ++iz)
      {
        const int cellIndex = merged.get_cell_index(iphi, ir, iz);
        const int cellEntries = merged.get_entries(cellIndex);
        if (cellEntries < static_cast<int>(min_entries_per_voxel))
        {
          ++skipped_low_entry_voxels;
          continue;
        }

        out_iphi = iphi;
        out_ir = ir;
        out_iz = iz;
        entries = cellEntries;

        const double phi = (iphi + 0.5) * 2.0 * pi / phiBins;
        const double radius = rMin + (ir + 0.5) * (rMax - rMin) / rBins;
        voxel_x = radius * std::cos(phi);
        voxel_y = radius * std::sin(phi);
        voxel_z = zMin + (iz + 0.5) * (zMax - zMin) / zBins;

        sum_pair_weight = merged.get_sum_weight(cellIndex);
        mean_pair_weight = CPMB2Container::mean(sum_pair_weight, cellEntries);
        effective_pair_entries = merged.get_sum_weight2(cellIndex) > 0.0 ?
            sum_pair_weight * sum_pair_weight / merged.get_sum_weight2(cellIndex) :
            std::numeric_limits<double>::quiet_NaN();
        mean_delta_r = merged.get_mean_delta_r(cellIndex, use_pair_weights);
        rms_delta_r = use_pair_weights ?
            CPMB2Container::weighted_rms(
                merged.get_sum_weighted_delta_r(cellIndex),
                merged.get_sum_weighted_delta_r2(cellIndex),
                sum_pair_weight) :
            CPMB2Container::rms(
                merged.get_sum_delta_r(cellIndex),
                merged.get_sum_delta_r2(cellIndex),
                cellEntries);
        mean_delta_rphi = merged.get_mean_delta_rphi(cellIndex, use_pair_weights);
        rms_delta_rphi = use_pair_weights ?
            CPMB2Container::weighted_rms(
                merged.get_sum_weighted_delta_rphi(cellIndex),
                merged.get_sum_weighted_delta_rphi2(cellIndex),
                sum_pair_weight) :
            CPMB2Container::rms(
                merged.get_sum_delta_rphi(cellIndex),
                merged.get_sum_delta_rphi2(cellIndex),
                cellEntries);
        mean_delta_phi = merged.get_mean_delta_phi(cellIndex, use_pair_weights);
        rms_delta_phi = use_pair_weights ?
            CPMB2Container::weighted_rms(
                merged.get_sum_weighted_delta_phi(cellIndex),
                merged.get_sum_weighted_delta_phi2(cellIndex),
                sum_pair_weight) :
            CPMB2Container::rms(
                merged.get_sum_delta_phi(cellIndex),
                merged.get_sum_delta_phi2(cellIndex),
                cellEntries);
        mean_delta_z = merged.get_mean_delta_z(cellIndex, use_pair_weights);
        rms_delta_z = use_pair_weights ?
            CPMB2Container::weighted_rms(
                merged.get_sum_weighted_delta_z(cellIndex),
                merged.get_sum_weighted_delta_z2(cellIndex),
                sum_pair_weight) :
            CPMB2Container::rms(
                merged.get_sum_delta_z(cellIndex),
                merged.get_sum_delta_z2(cellIndex),
                cellEntries);
        mean_dca = CPMB2Container::mean(merged.get_sum_dca(cellIndex), cellEntries);
        rms_dca = CPMB2Container::rms(merged.get_sum_dca(cellIndex), merged.get_sum_dca2(cellIndex), cellEntries);

        voxels.Fill();
        ++filled_voxels;
      }
    }
  }

  TTree metadata("cpm_metadata", "CPM merged correction metadata");
  int meta_phi_bins = phiBins;
  int meta_r_bins = rBins;
  int meta_z_bins = zBins;
  double meta_r_min = rMin;
  double meta_r_max = rMax;
  double meta_z_min = zMin;
  double meta_z_max = zMax;
  metadata.Branch("phi_bins", &meta_phi_bins);
  metadata.Branch("r_bins", &meta_r_bins);
  metadata.Branch("z_bins", &meta_z_bins);
  metadata.Branch("r_min", &meta_r_min);
  metadata.Branch("r_max", &meta_r_max);
  metadata.Branch("z_min", &meta_z_min);
  metadata.Branch("z_max", &meta_z_max);
  metadata.Fill();

  TTree summary("cpm_b2_container_summary", "CPM B2 container merge summary");
  unsigned int input_files_count = input_files.size();
  bool summary_use_pair_weights = use_pair_weights;
  unsigned int summary_min_entries_per_voxel = min_entries_per_voxel;
  summary.Branch("input_files", &input_files_count);
  summary.Branch("loaded_files", &loaded_files);
  summary.Branch("rejected_files", &rejected_files);
  summary.Branch("filled_voxels", &filled_voxels);
  summary.Branch("skipped_low_entry_voxels", &skipped_low_entry_voxels);
  summary.Branch("use_pair_weights", &summary_use_pair_weights);
  summary.Branch("min_entries_per_voxel", &summary_min_entries_per_voxel);
  summary.Fill();

  output.cd();
  merged.Write(object_name.c_str());
  voxels.Write();
  metadata.Write();
  summary.Write();
  output.Close();

  std::cout << "CPM_B2_MergeCorrectionContainers - loaded files: " << loaded_files << std::endl;
  std::cout << "CPM_B2_MergeCorrectionContainers - rejected files: " << rejected_files << std::endl;
  std::cout << "CPM_B2_MergeCorrectionContainers - filled voxels: " << filled_voxels << std::endl;
  std::cout << "CPM_B2_MergeCorrectionContainers - averaging mode: "
            << (use_pair_weights ? "weighted" : "plain") << std::endl;
  std::cout << "CPM_B2_MergeCorrectionContainers - output: " << output_file << std::endl;
}

void CPM_B2_MergeCorrectionContainers(
    const std::string& input_file,
    const std::string& output_file = "CPM_B2_voxel_corrections.root",
    const std::string& object_name = "CPMCorrectionContainer",
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1)
{
  CPM_B2_MergeCorrectionContainers(
      std::vector<std::string>{input_file},
      output_file,
      object_name,
      use_pair_weights,
      min_entries_per_voxel);
}

void CPM_B2_MergeCorrectionContainers(
    const std::string& input_file_or_list,
    const std::string& output_file,
    const std::string& object_name,
    const bool use_pair_weights,
    const unsigned int min_entries_per_voxel,
    const bool input_is_list)
{
  CPM_B2_MergeCorrectionContainers(
      input_is_list ? CPMB2Container::read_file_list(input_file_or_list) : std::vector<std::string>{input_file_or_list},
      output_file,
      object_name,
      use_pair_weights,
      min_entries_per_voxel);
}

void CPM_B2_MergeCorrectionContainersFromList(
    const std::string& input_list,
    const std::string& output_file = "CPM_B2_voxel_corrections.root",
    const std::string& object_name = "CPMCorrectionContainer",
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1)
{
  CPM_B2_MergeCorrectionContainers(
      CPMB2Container::read_file_list(input_list),
      output_file,
      object_name,
      use_pair_weights,
      min_entries_per_voxel);
}
