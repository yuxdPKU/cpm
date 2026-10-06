/*
 * CPM production Job B average-correction calculation.
 *
 * This macro only configures and runs the module-level reconstruction class.
 * The CPM pair selection, PoCA calculation, voxel accumulation, and histogram
 * writing live in CPMAverageCorrectionReconstruction.
 */

#include <CPMAverageCorrectionReconstruction.h>
#include <CPMReconstructionHelper.h>
#include <CPMVoxelContainer.h>

#include <TFile.h>

#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libcpm.so)
R__LOAD_LIBRARY(libtpccalib.so)

bool CPM_ComputeAverageCorrection(
    const std::vector<std::string>& input_files,
    const std::string& output_file = "CPM_B3_average_correction_histograms.root",
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1,
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 0,
    const unsigned int min_records_per_charge = 2,
    const double min_pair_pt = 0.5,
    const unsigned int max_pair_records_per_charge_batch = 10,
    const std::string& crossing_solver = "helix",
    const double magnetic_field_z = 1.4,
    const std::string& metadata_file = "",
    const unsigned long long max_input_records_per_chunk = 500000,
    const bool skip_invalid_inputs = false,
    const double max_abs_path = 5.0,
    const double max_midpoint_distance = 5.0,
    const bool allow_line_fallback = false)
{
  CPMAverageCorrectionReconstruction reconstruction;
  reconstruction.set_use_pair_weights(use_pair_weights);
  reconstruction.set_min_entries_per_voxel(min_entries_per_voxel);
  reconstruction.set_max_pair_dca(max_pair_dca);
  reconstruction.set_min_sin_angle(min_sin_angle);
  reconstruction.set_max_records_per_voxel(max_records_per_voxel);
  reconstruction.set_min_records_per_charge(min_records_per_charge);
  reconstruction.set_min_pair_pt(min_pair_pt);
  reconstruction.set_max_pair_records_per_charge_batch(max_pair_records_per_charge_batch);
  reconstruction.set_magnetic_field_z(magnetic_field_z);
  reconstruction.set_max_abs_path(max_abs_path);
  reconstruction.set_max_midpoint_distance(max_midpoint_distance);
  reconstruction.set_allow_line_fallback(allow_line_fallback);
  if (!reconstruction.set_crossing_solver(crossing_solver))
  {
    return false;
  }

  std::cout << "CPM_ComputeAverageCorrection - input files: "
            << input_files.size() << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - max input records per chunk: "
            << max_input_records_per_chunk << std::endl;
  std::cout << "CPM_ComputeAverageCorrection - skip invalid inputs: "
            << skip_invalid_inputs << std::endl;

  std::vector<std::string> skipped_invalid_input_files;
  std::vector<std::string> skipped_empty_input_files;
  std::size_t loaded_input_files = 0;

  for (std::size_t ifile = 0; ifile < input_files.size(); ++ifile)
  {
    const auto& input_file = input_files[ifile];
    const bool print_progress =
        ifile == 0 ||
        (ifile + 1) % 25 == 0 ||
        (ifile + 1) == input_files.size();

    if (print_progress)
    {
      std::cout << "CPM_ComputeAverageCorrection - loading input "
                << (ifile + 1) << "/" << input_files.size()
                << ": " << input_file << std::endl;
    }

    std::unique_ptr<TFile> input(TFile::Open(input_file.c_str(), "READ"));
    auto* source = input && !input->IsZombie() ?
        dynamic_cast<CPMVoxelContainer*>(input->Get("CPMVoxelContainer")) : nullptr;

    if (!source)
    {
      std::cout << "CPM_ComputeAverageCorrection - failed while loading input "
                << (ifile + 1) << "/" << input_files.size()
                << ": " << input_file << std::endl;
      if (!skip_invalid_inputs)
      {
        return false;
      }
      skipped_invalid_input_files.push_back(input_file);
      std::cout << "CPM_ComputeAverageCorrection - skipping invalid input: "
                << input_file << std::endl;
      continue;
    }

    if (source->empty())
    {
      skipped_empty_input_files.push_back(input_file);
      std::cout << "CPM_ComputeAverageCorrection - skipping empty input: "
                << input_file << std::endl;
      continue;
    }

    if (!reconstruction.add(*source))
    {
      std::cout << "CPM_ComputeAverageCorrection - failed while adding input "
                << (ifile + 1) << "/" << input_files.size()
                << ": " << input_file << std::endl;
      if (!skip_invalid_inputs)
      {
        return false;
      }
      skipped_invalid_input_files.push_back(input_file);
      std::cout << "CPM_ComputeAverageCorrection - skipping invalid input: "
                << input_file << std::endl;
      continue;
    }
    ++loaded_input_files;

    if (print_progress)
    {
      std::cout << "CPM_ComputeAverageCorrection - loaded inputs: "
                << (ifile + 1) << "/" << input_files.size()
                << " records: " << reconstruction.records().record_count()
                << " voxels: " << reconstruction.records().voxel_count()
                << std::endl;
    }

    if (max_input_records_per_chunk > 0 &&
        reconstruction.records().record_count() >= max_input_records_per_chunk)
    {
      std::cout << "CPM_ComputeAverageCorrection - processing loaded records after input "
                << (ifile + 1) << "/" << input_files.size()
                << " records: " << reconstruction.records().record_count()
                << " voxels: " << reconstruction.records().voxel_count()
                << std::endl;
      if (!reconstruction.process_loaded_records())
      {
        return false;
      }
      std::cout << "CPM_ComputeAverageCorrection - accumulated accepted pairs: "
                << reconstruction.summary().accepted_pairs
                << " accumulator voxels: "
                << reconstruction.summary().accumulator_voxels
                << std::endl;
    }
  }

  if (!skipped_empty_input_files.empty())
  {
    std::cout << "CPM_ComputeAverageCorrection - skipped empty input files: "
              << skipped_empty_input_files.size() << std::endl;
    for (const auto& filename : skipped_empty_input_files)
    {
      std::cout << "  " << filename << std::endl;
    }
  }

  if (!skipped_invalid_input_files.empty())
  {
    std::cout << "CPM_ComputeAverageCorrection - skipped invalid input files: "
              << skipped_invalid_input_files.size() << std::endl;
    for (const auto& filename : skipped_invalid_input_files)
    {
      std::cout << "  " << filename << std::endl;
    }
  }

  if (loaded_input_files == 0)
  {
    std::cout << "CPM_ComputeAverageCorrection - no valid input files were loaded"
              << std::endl;
    return false;
  }

  if (!metadata_file.empty())
  {
    CPMAverageCorrectionReconstruction metadata_check;
    if (!metadata_check.add_from_file(metadata_file))
    {
      if (!skip_invalid_inputs)
      {
        return false;
      }
      std::cout << "CPM_ComputeAverageCorrection - ignoring invalid metadata file: "
                << metadata_file << std::endl;
    }
    else if (!CPMReconstructionHelper::same_grid(reconstruction.grid(), metadata_check.grid()))
    {
      std::cout << "CPM_ComputeAverageCorrection - provided metadata file has inconsistent CPMVoxelContainer grid: "
                << metadata_file << std::endl;
      return false;
    }
  }

  std::cout << "CPM_ComputeAverageCorrection - processing final loaded records"
            << " records: " << reconstruction.records().record_count()
            << " voxels: " << reconstruction.records().voxel_count()
            << std::endl;
  if (!reconstruction.process_loaded_records())
  {
    return false;
  }

  std::cout << "CPM_ComputeAverageCorrection - finalizing average corrections"
            << std::endl;
  if (!reconstruction.finalize_average_corrections())
  {
    return false;
  }

  std::cout << "CPM_ComputeAverageCorrection - saving average corrections"
            << std::endl;
  if (!reconstruction.save_average_corrections(output_file))
  {
    return false;
  }

  reconstruction.print_summary();
  std::cout << "CPM_ComputeAverageCorrection - output: " << output_file << std::endl;
  return true;
}

bool CPM_ComputeAverageCorrection(
    const std::string& input_file_or_list,
    const std::string& output_file = "CPM_B3_average_correction_histograms.root",
    const bool input_is_list = false,
    const bool use_pair_weights = true,
    const unsigned int min_entries_per_voxel = 1,
    const double max_pair_dca = 2.0,
    const double min_sin_angle = 1.0e-4,
    const unsigned int max_records_per_voxel = 0,
    const unsigned int min_records_per_charge = 2,
    const double min_pair_pt = 0.5,
    const unsigned int max_pair_records_per_charge_batch = 10,
    const std::string& crossing_solver = "helix",
    const double magnetic_field_z = 1.4,
    const std::string& metadata_file = "",
    const unsigned long long max_input_records_per_chunk = 500000,
    const bool skip_invalid_inputs = false,
    const double max_abs_path = 5.0,
    const double max_midpoint_distance = 5.0,
    const bool allow_line_fallback = false)
{
  const auto input_files = input_is_list ?
      CPMReconstructionHelper::read_file_list(input_file_or_list) :
      std::vector<std::string>{input_file_or_list};

  return CPM_ComputeAverageCorrection(
      input_files,
      output_file,
      use_pair_weights,
      min_entries_per_voxel,
      max_pair_dca,
      min_sin_angle,
      max_records_per_voxel,
      min_records_per_charge,
      min_pair_pt,
      max_pair_records_per_charge_batch,
      crossing_solver,
      magnetic_field_z,
      metadata_file,
      max_input_records_per_chunk,
      skip_invalid_inputs,
      max_abs_path,
      max_midpoint_distance,
      allow_line_fallback);
}
