/*
 * Merge CPM production partial average-correction sums.
 *
 * Each partial file is produced by CPM_ComputeAverageCorrection and contains
 * cpm_voxel_correction_sums plus cpm_b3_summary. This macro merges sums first
 * and then reconstructs the final B3 histograms, so weighted and unweighted
 * runs remain exact.
 */

#include <CPMAverageCorrectionReconstruction.h>
#include <CPMReconstructionHelper.h>

#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libcpm.so)
R__LOAD_LIBRARY(libtpccalib.so)

bool CPM_MergeAverageCorrectionSums(
    const std::vector<std::string>& input_files,
    const std::string& output_file = "CPM_B3_average_correction_histograms.root",
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1)
{
  CPMAverageCorrectionReconstruction reconstruction;
  reconstruction.set_use_pair_weights(use_pair_weights);
  reconstruction.set_min_entries_per_voxel(min_entries_per_voxel);

  std::cout << "CPM_MergeAverageCorrectionSums - input files: "
            << input_files.size() << std::endl;

  for (std::size_t ifile = 0; ifile < input_files.size(); ++ifile)
  {
    const auto& input_file = input_files[ifile];
    const bool print_progress =
        ifile == 0 ||
        (ifile + 1) % 25 == 0 ||
        (ifile + 1) == input_files.size();
    if (print_progress)
    {
      std::cout << "CPM_MergeAverageCorrectionSums - loading partial "
                << (ifile + 1) << "/" << input_files.size()
                << ": " << input_file << std::endl;
    }
    if (!reconstruction.add_accumulators_from_file(input_file))
    {
      std::cout << "CPM_MergeAverageCorrectionSums - failed while loading partial "
                << (ifile + 1) << "/" << input_files.size()
                << ": " << input_file << std::endl;
      return false;
    }
  }

  std::cout << "CPM_MergeAverageCorrectionSums - finalizing average corrections"
            << std::endl;
  if (!reconstruction.finalize_average_corrections())
  {
    return false;
  }

  std::cout << "CPM_MergeAverageCorrectionSums - saving average corrections"
            << std::endl;
  if (!reconstruction.save_average_corrections(output_file))
  {
    return false;
  }

  reconstruction.print_summary();
  std::cout << "CPM_MergeAverageCorrectionSums - output: " << output_file << std::endl;
  return true;
}

bool CPM_MergeAverageCorrectionSumsFromList(
    const std::string& input_list,
    const std::string& output_file = "CPM_B3_average_correction_histograms.root",
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1)
{
  return CPM_MergeAverageCorrectionSums(
      CPMReconstructionHelper::read_file_list(input_list),
      output_file,
      use_pair_weights,
      min_entries_per_voxel);
}
