/*
 * CPM Job B3 average-correction histogram writer.
 *
 * This macro converts B2 voxel correction rows into the histogram names and
 * half-TPC/guard-bin layout consumed by TpcDistortionCorrectionContainer.
 * The correction convention is the same as TpcDistortionCorrection: values are
 * subtracted from cluster coordinates, so the B2 input must be
 * voxel center - crossing point.
 */

#include <tpccalib/TpcSpaceChargeReconstructionHelper.h>

#include <TFile.h>
#include <TH3.h>
#include <TH3F.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <tuple>

R__LOAD_LIBRARY(libtpccalib.so)

namespace CPMB3
{
  constexpr double kPi = 3.14159265358979323846;

  struct GridConfig
  {
    int phi_bins = 36;
    int r_bins = 16;
    int z_bins = 80;
    double phi_min = 0.0;
    double phi_max = 2.0 * kPi;
    double r_min = 20.0;
    double r_max = 78.0;
    double z_min = -105.5;
    double z_max = 105.5;
    bool from_metadata = false;
  };

  bool load_grid_metadata(const std::string& metadata_file, GridConfig& grid)
  {
    if (metadata_file.empty())
    {
      return false;
    }

    std::unique_ptr<TFile> input(TFile::Open(metadata_file.c_str(), "READ"));
    if (!input || input->IsZombie())
    {
      std::cout << "CPM_B3_WriteAverageCorrectionHistograms - could not open metadata file: "
                << metadata_file << std::endl;
      return false;
    }

    auto* metadata = dynamic_cast<TTree*>(input->Get("cpm_metadata"));
    if (!metadata || metadata->GetEntries() <= 0)
    {
      std::cout << "CPM_B3_WriteAverageCorrectionHistograms - no cpm_metadata tree in: "
                << metadata_file << std::endl;
      return false;
    }

    metadata->SetBranchAddress("phi_bins", &grid.phi_bins);
    metadata->SetBranchAddress("r_bins", &grid.r_bins);
    metadata->SetBranchAddress("z_bins", &grid.z_bins);
    metadata->SetBranchAddress("r_min", &grid.r_min);
    metadata->SetBranchAddress("r_max", &grid.r_max);
    metadata->SetBranchAddress("z_min", &grid.z_min);
    metadata->SetBranchAddress("z_max", &grid.z_max);
    metadata->GetEntry(0);

    grid.phi_min = 0.0;
    grid.phi_max = 2.0 * kPi;
    grid.from_metadata = true;
    return true;
  }

  std::tuple<std::unique_ptr<TH3>, std::unique_ptr<TH3>> finish_histogram(TH3* source, const TString& name)
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

void CPM_B3_WriteAverageCorrectionHistograms(
    const std::string& input_file = "CPM_B2_voxel_corrections.root",
    const std::string& output_file = "CPM_B3_average_correction_histograms.root",
    const std::string& metadata_file = "",
    const unsigned int min_entries_per_voxel = 1)
{
  CPMB3::GridConfig grid;
  CPMB3::load_grid_metadata(metadata_file, grid);

  if (grid.phi_bins <= 0 || grid.r_bins <= 0 || grid.z_bins <= 0 ||
      grid.phi_min >= grid.phi_max || grid.r_min >= grid.r_max || grid.z_min >= grid.z_max)
  {
    std::cout << "CPM_B3_WriteAverageCorrectionHistograms - invalid grid configuration" << std::endl;
    return;
  }

  std::unique_ptr<TFile> input(TFile::Open(input_file.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    std::cout << "CPM_B3_WriteAverageCorrectionHistograms - could not open input: "
              << input_file << std::endl;
    return;
  }

  auto* tree = dynamic_cast<TTree*>(input->Get("cpm_voxel_corrections"));
  if (!tree)
  {
    std::cout << "CPM_B3_WriteAverageCorrectionHistograms - missing cpm_voxel_corrections in: "
              << input_file << std::endl;
    return;
  }

  int iphi = -1;
  int ir = -1;
  int iz = -1;
  unsigned long long entries = 0;
  double mean_delta_r = std::numeric_limits<double>::quiet_NaN();
  double mean_delta_phi = std::numeric_limits<double>::quiet_NaN();
  double mean_delta_z = std::numeric_limits<double>::quiet_NaN();

  tree->SetBranchAddress("iphi", &iphi);
  tree->SetBranchAddress("ir", &ir);
  tree->SetBranchAddress("iz", &iz);
  tree->SetBranchAddress("entries", &entries);
  tree->SetBranchAddress("mean_delta_r", &mean_delta_r);
  tree->SetBranchAddress("mean_delta_phi", &mean_delta_phi);
  tree->SetBranchAddress("mean_delta_z", &mean_delta_z);

  TH3F hentries_rec(
      "hentries_rec",
      "CPM voxel entries;phi;r;z",
      grid.phi_bins, grid.phi_min, grid.phi_max,
      grid.r_bins, grid.r_min, grid.r_max,
      grid.z_bins, grid.z_min, grid.z_max);
  TH3F hdistortion_r_rec(
      "hDistortionR_rec",
      "CPM radial distortion;phi;r;z",
      grid.phi_bins, grid.phi_min, grid.phi_max,
      grid.r_bins, grid.r_min, grid.r_max,
      grid.z_bins, grid.z_min, grid.z_max);
  TH3F hdistortion_p_rec(
      "hDistortionP_rec",
      "CPM phi distortion;phi;r;z",
      grid.phi_bins, grid.phi_min, grid.phi_max,
      grid.r_bins, grid.r_min, grid.r_max,
      grid.z_bins, grid.z_min, grid.z_max);
  TH3F hdistortion_z_rec(
      "hDistortionZ_rec",
      "CPM z distortion;phi;r;z",
      grid.phi_bins, grid.phi_min, grid.phi_max,
      grid.r_bins, grid.r_min, grid.r_max,
      grid.z_bins, grid.z_min, grid.z_max);

  unsigned long long input_voxels = tree->GetEntries();
  unsigned long long filled_voxels = 0;
  unsigned long long skipped_low_entry_voxels = 0;
  unsigned long long skipped_invalid_voxels = 0;

  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    tree->GetEntry(entry);

    if (entries < min_entries_per_voxel)
    {
      ++skipped_low_entry_voxels;
      continue;
    }
    if (iphi < 0 || iphi >= grid.phi_bins ||
        ir < 0 || ir >= grid.r_bins ||
        iz < 0 || iz >= grid.z_bins ||
        !std::isfinite(mean_delta_r) ||
        !std::isfinite(mean_delta_phi) ||
        !std::isfinite(mean_delta_z))
    {
      ++skipped_invalid_voxels;
      continue;
    }

    hentries_rec.SetBinContent(iphi + 1, ir + 1, iz + 1, static_cast<double>(entries));
    hdistortion_r_rec.SetBinContent(iphi + 1, ir + 1, iz + 1, mean_delta_r);
    hdistortion_p_rec.SetBinContent(iphi + 1, ir + 1, iz + 1, mean_delta_phi);
    hdistortion_z_rec.SetBinContent(iphi + 1, ir + 1, iz + 1, mean_delta_z);
    ++filled_voxels;
  }

  TFile output(output_file.c_str(), "RECREATE");
  if (output.IsZombie())
  {
    std::cout << "CPM_B3_WriteAverageCorrectionHistograms - could not create output: "
              << output_file << std::endl;
    return;
  }

  hentries_rec.Write();
  hdistortion_r_rec.Write();
  hdistortion_p_rec.Write();
  hdistortion_z_rec.Write();

  CPMB3::write_finished_histogram(output, &hentries_rec, "hentries");
  CPMB3::write_finished_histogram(output, &hdistortion_r_rec, "hIntDistortionR");
  CPMB3::write_finished_histogram(output, &hdistortion_p_rec, "hIntDistortionP");
  CPMB3::write_finished_histogram(output, &hdistortion_z_rec, "hIntDistortionZ");

  TTree summary("cpm_b3_summary", "CPM B3 average-correction histogram summary");
  int phi_bins = grid.phi_bins;
  int r_bins = grid.r_bins;
  int z_bins = grid.z_bins;
  double phi_min = grid.phi_min;
  double phi_max = grid.phi_max;
  double r_min = grid.r_min;
  double r_max = grid.r_max;
  double z_min = grid.z_min;
  double z_max = grid.z_max;
  bool grid_from_metadata = grid.from_metadata;
  unsigned int summary_min_entries_per_voxel = min_entries_per_voxel;

  summary.Branch("input_voxels", &input_voxels);
  summary.Branch("filled_voxels", &filled_voxels);
  summary.Branch("skipped_low_entry_voxels", &skipped_low_entry_voxels);
  summary.Branch("skipped_invalid_voxels", &skipped_invalid_voxels);
  summary.Branch("min_entries_per_voxel", &summary_min_entries_per_voxel);
  summary.Branch("phi_bins", &phi_bins);
  summary.Branch("r_bins", &r_bins);
  summary.Branch("z_bins", &z_bins);
  summary.Branch("phi_min", &phi_min);
  summary.Branch("phi_max", &phi_max);
  summary.Branch("r_min", &r_min);
  summary.Branch("r_max", &r_max);
  summary.Branch("z_min", &z_min);
  summary.Branch("z_max", &z_max);
  summary.Branch("grid_from_metadata", &grid_from_metadata);
  summary.Fill();
  summary.Write();

  output.Close();

  std::cout << "CPM_B3_WriteAverageCorrectionHistograms - input voxels: " << input_voxels << std::endl;
  std::cout << "CPM_B3_WriteAverageCorrectionHistograms - filled voxels: " << filled_voxels << std::endl;
  std::cout << "CPM_B3_WriteAverageCorrectionHistograms - skipped low-entry voxels: "
            << skipped_low_entry_voxels << std::endl;
  std::cout << "CPM_B3_WriteAverageCorrectionHistograms - skipped invalid voxels: "
            << skipped_invalid_voxels << std::endl;
  std::cout << "CPM_B3_WriteAverageCorrectionHistograms - grid: ("
            << grid.phi_bins << ", " << grid.r_bins << ", " << grid.z_bins << ")"
            << " r=[" << grid.r_min << ", " << grid.r_max << "]"
            << " z=[" << grid.z_min << ", " << grid.z_max << "]"
            << (grid.from_metadata ? " from metadata" : " default")
            << std::endl;
  std::cout << "CPM_B3_WriteAverageCorrectionHistograms - output: " << output_file << std::endl;
}
