/*
 * CPM QA offline-diagnostic driver.
 *
 * This macro keeps the record-based diagnostic path in one user-facing entry
 * point while still writing separate stage ROOT files for debugging:
 *   QA B0 event index, QA B1 PoCA, QA B2 voxel corrections, and B3 histograms.
 */

#include "CPM_QA_B0_BuildEventIndex.C"
#include "CPM_QA_B0_CheckEventIndex.C"
#include "CPM_QA_B1_ComputePoCA.C"
#include "CPM_QA_B2_AccumulateVoxelCorrections.C"
#include "CPM_QA_B3_WriteAverageCorrectionHistograms.C"
#include "CPM_QA_B3_CheckAverageCorrectionHistograms.C"

#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <string>

namespace CPMQADiagnostics
{
  std::string join_path(const std::string& directory, const std::string& filename)
  {
    if (directory.empty() || directory == ".")
    {
      return filename;
    }
    if (directory.back() == '/')
    {
      return directory + filename;
    }
    return directory + "/" + filename;
  }

  std::string first_list_entry(const std::string& input_list)
  {
    std::ifstream input(input_list);
    std::string line;
    while (std::getline(input, line))
    {
      if (line.empty() || line[0] == '#')
      {
        continue;
      }
      return line;
    }
    return "";
  }
}

bool CPM_QA_RunOfflineDiagnostics(
    const std::string& input_file_or_list,
    const std::string& output_dir = ".",
    const std::string& prefix = "CPM",
    const bool input_is_list = false,
    const bool run_b0_qa = false,
    const double b1_max_pair_dca = 2.0,
    const double b1_min_sin_angle = 1.0e-4,
    const unsigned int b1_max_records_per_voxel = 0,
    const unsigned int b1_min_records_per_charge = 2,
    const bool b1_print_voxel_summaries = true,
    const double b1_min_pair_pt = 0.5,
    const unsigned int b1_max_pair_records_per_voxel = 10,
    const std::string& b1_crossing_solver = "helix",
    const double b1_magnetic_field_z = 1.4,
    const bool b1_write_pair_tree = true,
    const unsigned int b2_min_entries_per_voxel = 1,
    const double b2_max_pair_dca = -1.0,
    const bool use_pair_weights = true,
    const std::string& b2_input_mode = "auto",
    const unsigned int b3_min_entries_per_voxel = 1,
    const std::string& metadata_file = "")
{
  if (gSystem && !output_dir.empty())
  {
    gSystem->mkdir(output_dir.c_str(), true);
  }

  const std::string b0_output =
      CPMQADiagnostics::join_path(output_dir, prefix + "_QA_B0_event_index.root");
  const std::string b1_output =
      CPMQADiagnostics::join_path(output_dir, prefix + "_QA_B1_poca.root");
  const std::string b2_output =
      CPMQADiagnostics::join_path(output_dir, prefix + "_QA_B2_voxel_corrections.root");
  const std::string b3_output =
      CPMQADiagnostics::join_path(output_dir, prefix + "_B3_average_correction_histograms.root");

  const std::string resolved_metadata = !metadata_file.empty() ?
      metadata_file :
      (input_is_list ? CPMQADiagnostics::first_list_entry(input_file_or_list) : input_file_or_list);

  bool ok = true;

  if (run_b0_qa)
  {
    CPM_QA_B0_BuildEventIndex(input_file_or_list, b0_output, input_is_list);
    ok = CPM_QA_B0_CheckEventIndex(b0_output) && ok;
  }

  CPM_QA_B1_ComputePoCA(
      input_file_or_list,
      b1_output,
      input_is_list,
      b1_max_pair_dca,
      b1_min_sin_angle,
      b1_max_records_per_voxel,
      b1_min_records_per_charge,
      b1_print_voxel_summaries,
      b1_min_pair_pt,
      b1_max_pair_records_per_voxel,
      b1_crossing_solver,
      b1_magnetic_field_z,
      b1_write_pair_tree);

  CPM_QA_B2_AccumulateVoxelCorrections(
      b1_output,
      b2_output,
      b2_min_entries_per_voxel,
      b2_max_pair_dca,
      use_pair_weights,
      b2_input_mode);

  CPM_QA_B3_WriteAverageCorrectionHistograms(
      b2_output,
      b3_output,
      resolved_metadata,
      b3_min_entries_per_voxel);
  ok = CPM_QA_B3_CheckAverageCorrectionHistograms(b3_output) && ok;

  std::cout << "CPM_QA_RunOfflineDiagnostics - B1 output: " << b1_output << std::endl;
  std::cout << "CPM_QA_RunOfflineDiagnostics - B2 output: " << b2_output << std::endl;
  std::cout << "CPM_QA_RunOfflineDiagnostics - B3 output: " << b3_output << std::endl;
  if (run_b0_qa)
  {
    std::cout << "CPM_QA_RunOfflineDiagnostics - B0 output: " << b0_output << std::endl;
  }
  return ok;
}
