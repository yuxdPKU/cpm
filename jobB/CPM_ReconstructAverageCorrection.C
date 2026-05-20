/*
 * CPM production average-correction reconstruction.
 *
 * This is the matrix-inversion-style entry point for the CPM production path:
 * it loads Job A CPMCorrectionContainer objects, merges them, converts the
 * merged sums into voxel corrections, and writes the final average-correction
 * histograms plus intermediate QA objects to one output ROOT file.
 */

#include <CPMAverageCorrectionReconstruction.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libcpm.so)

namespace CPMAverageCorrectionMacro
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
}

void CPM_ReconstructAverageCorrection(
    const std::vector<std::string>& input_files,
    const std::string& output_file = "CPM_average_correction.root",
    const std::string& object_name = "CPMCorrectionContainer",
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1)
{
  CPMAverageCorrectionReconstruction reconstruction;
  reconstruction.set_container_name(object_name);
  reconstruction.set_use_pair_weights(use_pair_weights);
  reconstruction.set_min_entries_per_voxel(min_entries_per_voxel);

  for (const auto& filename : input_files)
  {
    reconstruction.add_from_file(filename, object_name);
  }

  reconstruction.calculate_average_corrections();
  reconstruction.save_average_corrections(output_file);

  std::cout << "CPM_ReconstructAverageCorrection - input files: " << input_files.size() << std::endl;
  std::cout << "CPM_ReconstructAverageCorrection - averaging mode: "
            << (use_pair_weights ? "weighted" : "plain") << std::endl;
  std::cout << "CPM_ReconstructAverageCorrection - output: " << output_file << std::endl;
}

void CPM_ReconstructAverageCorrection(
    const std::string& input_file,
    const std::string& output_file = "CPM_average_correction.root",
    const std::string& object_name = "CPMCorrectionContainer",
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1)
{
  CPM_ReconstructAverageCorrection(
      std::vector<std::string>{input_file},
      output_file,
      object_name,
      use_pair_weights,
      min_entries_per_voxel);
}

void CPM_ReconstructAverageCorrection(
    const std::string& input_file_or_list,
    const std::string& output_file,
    const std::string& object_name,
    const bool use_pair_weights,
    const unsigned int min_entries_per_voxel,
    const bool input_is_list)
{
  CPM_ReconstructAverageCorrection(
      input_is_list ? CPMAverageCorrectionMacro::read_file_list(input_file_or_list) : std::vector<std::string>{input_file_or_list},
      output_file,
      object_name,
      use_pair_weights,
      min_entries_per_voxel);
}

void CPM_ReconstructAverageCorrectionFromList(
    const std::string& input_list,
    const std::string& output_file = "CPM_average_correction.root",
    const std::string& object_name = "CPMCorrectionContainer",
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1)
{
  CPM_ReconstructAverageCorrection(
      CPMAverageCorrectionMacro::read_file_list(input_list),
      output_file,
      object_name,
      use_pair_weights,
      min_entries_per_voxel);
}
