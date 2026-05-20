#include "CPMAverageCorrectionReconstruction.h"

#include <tpccalib/TpcSpaceChargeReconstructionHelper.h>

#include <TFile.h>
#include <TH3.h>
#include <TH3F.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <tuple>

namespace
{
  constexpr double kPi = 3.14159265358979323846;

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

CPMAverageCorrectionReconstruction::CPMAverageCorrectionReconstruction(const std::string& name)
  : Fun4AllBase(name)
{
}

CPMAverageCorrectionReconstruction::~CPMAverageCorrectionReconstruction() = default;

bool CPMAverageCorrectionReconstruction::add_from_file(
    const std::string& filename,
    const std::string& objectname)
{
  const std::string sourceObjectName = objectname.empty() ? m_containerName : objectname;
  std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    std::cout << "CPMAverageCorrectionReconstruction::add_from_file - could not open "
              << filename << std::endl;
    ++m_rejectedFiles;
    return false;
  }

  auto* source = dynamic_cast<CPMCorrectionContainer*>(input->Get(sourceObjectName.c_str()));
  if (!source)
  {
    std::cout << "CPMAverageCorrectionReconstruction::add_from_file - missing object "
              << sourceObjectName << " in " << filename << std::endl;
    ++m_rejectedFiles;
    return false;
  }

  if (!add(*source))
  {
    ++m_rejectedFiles;
    return false;
  }

  ++m_loadedFiles;
  return true;
}

bool CPMAverageCorrectionReconstruction::add(const CPMCorrectionContainer& source)
{
  if (!m_container)
  {
    int phiBins = 0;
    int rBins = 0;
    int zBins = 0;
    double rMin = 0.0;
    double rMax = 0.0;
    double zMin = 0.0;
    double zMax = 0.0;
    source.get_grid_dimensions(phiBins, rBins, zBins);
    source.get_grid_range(rMin, rMax, zMin, zMax);

    m_container = std::make_unique<CPMCorrectionContainer>();
    m_container->set_grid_dimensions(phiBins, rBins, zBins);
    m_container->set_grid_range(rMin, rMax, zMin, zMax);
  }

  reset_calculation();
  return m_container->add(source);
}

void CPMAverageCorrectionReconstruction::Reset()
{
  m_container.reset();
  reset_calculation();
  m_loadedFiles = 0;
  m_rejectedFiles = 0;
}

void CPMAverageCorrectionReconstruction::reset_calculation()
{
  m_voxelCorrections.clear();
  m_hentriesRec.reset();
  m_hdistortionRRec.reset();
  m_hdistortionPRec.reset();
  m_hdistortionZRec.reset();
  m_filledVoxels = 0;
  m_skippedLowEntryVoxels = 0;
  m_scannedVoxels = 0;
}

void CPMAverageCorrectionReconstruction::calculate_average_corrections()
{
  if (!m_container)
  {
    std::cout << "CPMAverageCorrectionReconstruction::calculate_average_corrections - no CPM containers loaded" << std::endl;
    return;
  }

  reset_calculation();

  int phiBins = 0;
  int rBins = 0;
  int zBins = 0;
  double rMin = 0.0;
  double rMax = 0.0;
  double zMin = 0.0;
  double zMax = 0.0;
  m_container->get_grid_dimensions(phiBins, rBins, zBins);
  m_container->get_grid_range(rMin, rMax, zMin, zMax);

  m_hentriesRec = std::make_unique<TH3F>(
      "hentries_rec",
      "CPM voxel entries;phi;r;z",
      phiBins, 0.0, 2.0 * kPi,
      rBins, rMin, rMax,
      zBins, zMin, zMax);
  m_hdistortionRRec = std::make_unique<TH3F>(
      "hDistortionR_rec",
      "CPM radial distortion;phi;r;z",
      phiBins, 0.0, 2.0 * kPi,
      rBins, rMin, rMax,
      zBins, zMin, zMax);
  m_hdistortionPRec = std::make_unique<TH3F>(
      "hDistortionP_rec",
      "CPM phi distortion;phi;r;z",
      phiBins, 0.0, 2.0 * kPi,
      rBins, rMin, rMax,
      zBins, zMin, zMax);
  m_hdistortionZRec = std::make_unique<TH3F>(
      "hDistortionZ_rec",
      "CPM z distortion;phi;r;z",
      phiBins, 0.0, 2.0 * kPi,
      rBins, rMin, rMax,
      zBins, zMin, zMax);

  for (int iphi = 0; iphi < phiBins; ++iphi)
  {
    for (int ir = 0; ir < rBins; ++ir)
    {
      for (int iz = 0; iz < zBins; ++iz)
      {
        ++m_scannedVoxels;
        const int cellIndex = m_container->get_cell_index(iphi, ir, iz);
        const int cellEntries = m_container->get_entries(cellIndex);
        if (cellEntries < static_cast<int>(m_minEntriesPerVoxel))
        {
          ++m_skippedLowEntryVoxels;
          continue;
        }

        VoxelCorrection correction;
        correction.iphi = iphi;
        correction.ir = ir;
        correction.iz = iz;
        correction.entries = cellEntries;

        const double phi = (iphi + 0.5) * 2.0 * kPi / phiBins;
        const double radius = rMin + (ir + 0.5) * (rMax - rMin) / rBins;
        correction.voxel_x = radius * std::cos(phi);
        correction.voxel_y = radius * std::sin(phi);
        correction.voxel_z = zMin + (iz + 0.5) * (zMax - zMin) / zBins;

        correction.sum_pair_weight = m_container->get_sum_weight(cellIndex);
        correction.mean_pair_weight = mean(correction.sum_pair_weight, cellEntries);
        correction.effective_pair_entries = m_container->get_sum_weight2(cellIndex) > 0.0 ?
            correction.sum_pair_weight * correction.sum_pair_weight / m_container->get_sum_weight2(cellIndex) :
            std::numeric_limits<double>::quiet_NaN();
        correction.mean_delta_r = m_container->get_mean_delta_r(cellIndex, m_usePairWeights);
        correction.rms_delta_r = m_usePairWeights ?
            weighted_rms(
                m_container->get_sum_weighted_delta_r(cellIndex),
                m_container->get_sum_weighted_delta_r2(cellIndex),
                correction.sum_pair_weight) :
            rms(
                m_container->get_sum_delta_r(cellIndex),
                m_container->get_sum_delta_r2(cellIndex),
                cellEntries);
        correction.mean_delta_rphi = m_container->get_mean_delta_rphi(cellIndex, m_usePairWeights);
        correction.rms_delta_rphi = m_usePairWeights ?
            weighted_rms(
                m_container->get_sum_weighted_delta_rphi(cellIndex),
                m_container->get_sum_weighted_delta_rphi2(cellIndex),
                correction.sum_pair_weight) :
            rms(
                m_container->get_sum_delta_rphi(cellIndex),
                m_container->get_sum_delta_rphi2(cellIndex),
                cellEntries);
        correction.mean_delta_phi = m_container->get_mean_delta_phi(cellIndex, m_usePairWeights);
        correction.rms_delta_phi = m_usePairWeights ?
            weighted_rms(
                m_container->get_sum_weighted_delta_phi(cellIndex),
                m_container->get_sum_weighted_delta_phi2(cellIndex),
                correction.sum_pair_weight) :
            rms(
                m_container->get_sum_delta_phi(cellIndex),
                m_container->get_sum_delta_phi2(cellIndex),
                cellEntries);
        correction.mean_delta_z = m_container->get_mean_delta_z(cellIndex, m_usePairWeights);
        correction.rms_delta_z = m_usePairWeights ?
            weighted_rms(
                m_container->get_sum_weighted_delta_z(cellIndex),
                m_container->get_sum_weighted_delta_z2(cellIndex),
                correction.sum_pair_weight) :
            rms(
                m_container->get_sum_delta_z(cellIndex),
                m_container->get_sum_delta_z2(cellIndex),
                cellEntries);
        correction.mean_dca = mean(m_container->get_sum_dca(cellIndex), cellEntries);
        correction.rms_dca = rms(m_container->get_sum_dca(cellIndex), m_container->get_sum_dca2(cellIndex), cellEntries);

        if (!std::isfinite(correction.mean_delta_r) ||
            !std::isfinite(correction.mean_delta_phi) ||
            !std::isfinite(correction.mean_delta_z))
        {
          continue;
        }

        m_hentriesRec->SetBinContent(iphi + 1, ir + 1, iz + 1, static_cast<double>(cellEntries));
        m_hdistortionRRec->SetBinContent(iphi + 1, ir + 1, iz + 1, correction.mean_delta_r);
        m_hdistortionPRec->SetBinContent(iphi + 1, ir + 1, iz + 1, correction.mean_delta_phi);
        m_hdistortionZRec->SetBinContent(iphi + 1, ir + 1, iz + 1, correction.mean_delta_z);

        m_voxelCorrections.push_back(correction);
        ++m_filledVoxels;
      }
    }
  }
}

void CPMAverageCorrectionReconstruction::save_average_corrections(const std::string& filename)
{
  if (!m_container)
  {
    std::cout << "CPMAverageCorrectionReconstruction::save_average_corrections - no CPM containers loaded" << std::endl;
    return;
  }

  if (!m_hentriesRec || !m_hdistortionRRec || !m_hdistortionPRec || !m_hdistortionZRec)
  {
    calculate_average_corrections();
  }

  TFile output(filename.c_str(), "RECREATE");
  if (output.IsZombie())
  {
    std::cout << "CPMAverageCorrectionReconstruction::save_average_corrections - could not create "
              << filename << std::endl;
    return;
  }

  output.cd();
  m_container->Write(m_containerName.c_str());
  write_voxel_corrections(output);
  write_metadata(output);
  write_histograms(output);
  write_summary(output);
  output.Close();
}

void CPMAverageCorrectionReconstruction::write_voxel_corrections(TFile& output) const
{
  TTree voxels("cpm_voxel_corrections", "CPM voxel-level correction QA");
  VoxelCorrection row;

  voxels.Branch("iphi", &row.iphi);
  voxels.Branch("ir", &row.ir);
  voxels.Branch("iz", &row.iz);
  voxels.Branch("entries", &row.entries);
  voxels.Branch("voxel_x", &row.voxel_x);
  voxels.Branch("voxel_y", &row.voxel_y);
  voxels.Branch("voxel_z", &row.voxel_z);
  voxels.Branch("sum_pair_weight", &row.sum_pair_weight);
  voxels.Branch("mean_pair_weight", &row.mean_pair_weight);
  voxels.Branch("effective_pair_entries", &row.effective_pair_entries);
  voxels.Branch("mean_delta_r", &row.mean_delta_r);
  voxels.Branch("rms_delta_r", &row.rms_delta_r);
  voxels.Branch("mean_delta_rphi", &row.mean_delta_rphi);
  voxels.Branch("rms_delta_rphi", &row.rms_delta_rphi);
  voxels.Branch("mean_delta_phi", &row.mean_delta_phi);
  voxels.Branch("rms_delta_phi", &row.rms_delta_phi);
  voxels.Branch("mean_delta_z", &row.mean_delta_z);
  voxels.Branch("rms_delta_z", &row.rms_delta_z);
  voxels.Branch("mean_dca", &row.mean_dca);
  voxels.Branch("rms_dca", &row.rms_dca);

  for (const auto& correction : m_voxelCorrections)
  {
    row = correction;
    voxels.Fill();
  }

  output.cd();
  voxels.Write();
}

void CPMAverageCorrectionReconstruction::write_metadata(TFile& output) const
{
  int phiBins = 0;
  int rBins = 0;
  int zBins = 0;
  double rMin = 0.0;
  double rMax = 0.0;
  double zMin = 0.0;
  double zMax = 0.0;
  m_container->get_grid_dimensions(phiBins, rBins, zBins);
  m_container->get_grid_range(rMin, rMax, zMin, zMax);

  TTree metadata("cpm_metadata", "CPM merged correction metadata");
  metadata.Branch("phi_bins", &phiBins);
  metadata.Branch("r_bins", &rBins);
  metadata.Branch("z_bins", &zBins);
  metadata.Branch("r_min", &rMin);
  metadata.Branch("r_max", &rMax);
  metadata.Branch("z_min", &zMin);
  metadata.Branch("z_max", &zMax);
  metadata.Fill();

  output.cd();
  metadata.Write();
}

void CPMAverageCorrectionReconstruction::write_summary(TFile& output) const
{
  TTree summary("cpm_average_correction_summary", "CPM average correction reconstruction summary");
  unsigned int loadedFiles = m_loadedFiles;
  unsigned int rejectedFiles = m_rejectedFiles;
  unsigned int scannedVoxels = m_scannedVoxels;
  unsigned int filledVoxels = m_filledVoxels;
  unsigned int skippedLowEntryVoxels = m_skippedLowEntryVoxels;
  bool usePairWeights = m_usePairWeights;
  unsigned int minEntriesPerVoxel = m_minEntriesPerVoxel;
  std::string averagingMode = m_usePairWeights ? "weighted" : "plain";

  summary.Branch("loaded_files", &loadedFiles);
  summary.Branch("rejected_files", &rejectedFiles);
  summary.Branch("scanned_voxels", &scannedVoxels);
  summary.Branch("filled_voxels", &filledVoxels);
  summary.Branch("skipped_low_entry_voxels", &skippedLowEntryVoxels);
  summary.Branch("use_pair_weights", &usePairWeights);
  summary.Branch("min_entries_per_voxel", &minEntriesPerVoxel);
  summary.Branch("averaging_mode", &averagingMode);
  summary.Fill();

  output.cd();
  summary.Write();
}

void CPMAverageCorrectionReconstruction::write_histograms(TFile& output) const
{
  if (!m_hentriesRec || !m_hdistortionRRec || !m_hdistortionPRec || !m_hdistortionZRec)
  {
    return;
  }

  output.cd();
  m_hentriesRec->Write();
  m_hdistortionRRec->Write();
  m_hdistortionPRec->Write();
  m_hdistortionZRec->Write();

  write_finished_histogram(output, m_hentriesRec.get(), "hentries");
  write_finished_histogram(output, m_hdistortionRRec.get(), "hIntDistortionR");
  write_finished_histogram(output, m_hdistortionPRec.get(), "hIntDistortionP");
  write_finished_histogram(output, m_hdistortionZRec.get(), "hIntDistortionZ");
}

double CPMAverageCorrectionReconstruction::mean(const double sum, const int entries)
{
  return entries > 0 ? sum / static_cast<double>(entries) : std::numeric_limits<double>::quiet_NaN();
}

double CPMAverageCorrectionReconstruction::rms(const double sum, const double sum2, const int entries)
{
  if (entries <= 0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double average = mean(sum, entries);
  const double variance = sum2 / static_cast<double>(entries) - average * average;
  return std::sqrt(std::max(0.0, variance));
}

double CPMAverageCorrectionReconstruction::weighted_rms(
    const double sumWeighted,
    const double sumWeighted2,
    const double sumWeight)
{
  if (sumWeight <= 0.0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double average = sumWeighted / sumWeight;
  const double variance = sumWeighted2 / sumWeight - average * average;
  return std::sqrt(std::max(0.0, variance));
}
