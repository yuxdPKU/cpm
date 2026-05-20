#ifndef CPM_CPMAVERAGECORRECTIONRECONSTRUCTION_H
#define CPM_CPMAVERAGECORRECTIONRECONSTRUCTION_H

#include "CPMCorrectionContainer.h"

#include <fun4all/Fun4AllBase.h>

#include <memory>
#include <string>
#include <vector>

class TFile;
class TH3;

class CPMAverageCorrectionReconstruction : public Fun4AllBase
{
 public:
  explicit CPMAverageCorrectionReconstruction(const std::string& name = "CPMAVERAGECORRECTIONRECONSTRUCTION");
  ~CPMAverageCorrectionReconstruction() override;

  void set_container_name(const std::string& value) { m_containerName = value; }
  void set_use_pair_weights(const bool value) { m_usePairWeights = value; }
  void set_min_entries_per_voxel(const unsigned int value) { m_minEntriesPerVoxel = value; }

  bool add(const CPMCorrectionContainer& source);
  bool add_from_file(const std::string& filename, const std::string& objectname = "");

  void calculate_average_corrections();
  void save_average_corrections(const std::string& filename = "CPM_average_correction.root");
  void Reset();

 private:
  struct VoxelCorrection
  {
    int iphi = -1;
    int ir = -1;
    int iz = -1;
    unsigned long long entries = 0;
    double voxel_x = 0.0;
    double voxel_y = 0.0;
    double voxel_z = 0.0;
    double sum_pair_weight = 0.0;
    double mean_pair_weight = 0.0;
    double effective_pair_entries = 0.0;
    double mean_delta_r = 0.0;
    double rms_delta_r = 0.0;
    double mean_delta_rphi = 0.0;
    double rms_delta_rphi = 0.0;
    double mean_delta_phi = 0.0;
    double rms_delta_phi = 0.0;
    double mean_delta_z = 0.0;
    double rms_delta_z = 0.0;
    double mean_dca = 0.0;
    double rms_dca = 0.0;
  };

  void reset_calculation();
  void write_voxel_corrections(TFile& output) const;
  void write_metadata(TFile& output) const;
  void write_summary(TFile& output) const;
  void write_histograms(TFile& output) const;

  static double mean(double sum, int entries);
  static double rms(double sum, double sum2, int entries);
  static double weighted_rms(double sumWeighted, double sumWeighted2, double sumWeight);

  std::unique_ptr<CPMCorrectionContainer> m_container;
  std::vector<VoxelCorrection> m_voxelCorrections;
  std::unique_ptr<TH3> m_hentriesRec;
  std::unique_ptr<TH3> m_hdistortionRRec;
  std::unique_ptr<TH3> m_hdistortionPRec;
  std::unique_ptr<TH3> m_hdistortionZRec;

  std::string m_containerName = "CPMCorrectionContainer";
  bool m_usePairWeights = true;
  unsigned int m_minEntriesPerVoxel = 1;
  unsigned int m_loadedFiles = 0;
  unsigned int m_rejectedFiles = 0;
  unsigned int m_filledVoxels = 0;
  unsigned int m_skippedLowEntryVoxels = 0;
  unsigned int m_scannedVoxels = 0;
};

#endif
