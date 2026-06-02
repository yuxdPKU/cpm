#ifndef CPM_CPMAVERAGECORRECTIONRECONSTRUCTION_H
#define CPM_CPMAVERAGECORRECTIONRECONSTRUCTION_H

#include "CPMReconstructionHelper.h"
#include "CPMVoxelContainerv1.h"

#include <iostream>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <string>

class TH3;
class TFile;

class CPMAverageCorrectionReconstruction
{
 public:
  struct Summary
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

  CPMAverageCorrectionReconstruction();
  ~CPMAverageCorrectionReconstruction();

  bool add(const CPMVoxelContainer& source);
  bool add_from_file(
      const std::string& filename,
      const std::string& objectname = "CPMVoxelContainer");

  bool process_loaded_records();
  bool finalize_average_corrections();
  bool calculate_average_corrections();
  bool save_average_corrections(
      const std::string& filename = "CPM_B3_average_correction_histograms.root") const;

  void print_summary(std::ostream& out = std::cout) const;

  void set_use_pair_weights(bool value) { m_use_pair_weights = value; }
  void set_min_entries_per_voxel(unsigned int value) { m_min_entries_per_voxel = value; }
  void set_max_pair_dca(double value) { m_pair_options.max_pair_dca = value; }
  void set_min_sin_angle(double value) { m_pair_options.min_sin_angle = value; }
  void set_max_records_per_voxel(unsigned int value) { m_max_records_per_voxel = value; }
  void set_min_records_per_charge(unsigned int value) { m_min_records_per_charge = value; }
  void set_min_pair_pt(double value)
  {
    m_min_pair_pt = value;
    m_pair_options.min_pt = value;
  }
  void set_max_pair_records_per_charge_batch(unsigned int value)
  {
    m_max_pair_records_per_charge_batch = value;
  }
  void set_magnetic_field_z(double value) { m_pair_options.magnetic_field_z = value; }
  void set_crossing_solver(CPMReconstructionHelper::PairSolver solver) { m_pair_options.solver = solver; }
  bool set_crossing_solver(const std::string& solver);

  [[nodiscard]] const Summary& summary() const { return m_summary; }
  [[nodiscard]] const CPMVoxelContainer& records() const { return m_records; }
  [[nodiscard]] const CPMVoxelContainer::Grid& grid() const { return m_records.grid(); }

 private:
  void reset_output();
  void reset_calculation_summary();
  void write_summary_tree(TFile& output) const;

  CPMVoxelContainerv1 m_records;
  std::map<VoxelId, CPMReconstructionHelper::CorrectionAccumulator> m_accumulators;
  std::set<VoxelId> m_input_voxels;
  Summary m_summary;

  bool m_use_pair_weights = true;
  unsigned int m_min_entries_per_voxel = 1;
  unsigned int m_max_records_per_voxel = 0;
  unsigned int m_min_records_per_charge = 2;
  double m_min_pair_pt = 0.5;
  unsigned int m_max_pair_records_per_charge_batch = 10;
  CPMReconstructionHelper::PairOptions m_pair_options;

  std::unique_ptr<TH3> m_hentries_rec;
  std::unique_ptr<TH3> m_hdistortion_r_rec;
  std::unique_ptr<TH3> m_hdistortion_p_rec;
  std::unique_ptr<TH3> m_hdistortion_z_rec;
};

#endif
