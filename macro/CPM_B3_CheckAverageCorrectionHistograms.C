/*
 * Light QA for CPM B3 average-correction histogram outputs.
 */

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

bool CPM_B3_CheckAverageCorrectionHistograms(
    const std::string& input_file = "CPM_B3_average_correction_histograms.root")
{
  TFile input(input_file.c_str(), "READ");
  if (input.IsZombie())
  {
    std::cout << "CPM_B3_CheckAverageCorrectionHistograms - could not open: "
              << input_file << std::endl;
    return false;
  }

  const std::vector<std::string> required_histograms = {
      "hentries_negz",
      "hentries_posz",
      "hIntDistortionR_negz",
      "hIntDistortionR_posz",
      "hIntDistortionP_negz",
      "hIntDistortionP_posz",
      "hIntDistortionZ_negz",
      "hIntDistortionZ_posz"};

  unsigned int missing_histograms = 0;
  unsigned int invalid_dimensions = 0;
  for (const auto& name : required_histograms)
  {
    auto* hist = dynamic_cast<TH1*>(input.Get(name.c_str()));
    if (!hist)
    {
      std::cout << "CPM_B3_CheckAverageCorrectionHistograms - missing: "
                << name << std::endl;
      ++missing_histograms;
      continue;
    }
    if (hist->GetDimension() != 3)
    {
      std::cout << "CPM_B3_CheckAverageCorrectionHistograms - invalid dimension for "
                << name << ": " << hist->GetDimension() << std::endl;
      ++invalid_dimensions;
    }
  }

  auto* summary = dynamic_cast<TTree*>(input.Get("cpm_b3_summary"));
  const bool has_summary = summary != nullptr && summary->GetEntries() > 0;

  std::cout << "CPM_B3_CheckAverageCorrectionHistograms - file: " << input_file << std::endl;
  std::cout << "CPM_B3_CheckAverageCorrectionHistograms - required histograms: "
            << required_histograms.size() << std::endl;
  std::cout << "CPM_B3_CheckAverageCorrectionHistograms - missing histograms: "
            << missing_histograms << std::endl;
  std::cout << "CPM_B3_CheckAverageCorrectionHistograms - invalid dimensions: "
            << invalid_dimensions << std::endl;
  std::cout << "CPM_B3_CheckAverageCorrectionHistograms - has summary: "
            << has_summary << std::endl;

  const bool ok = missing_histograms == 0 && invalid_dimensions == 0 && has_summary;
  std::cout << "CPM_B3_CheckAverageCorrectionHistograms - status: "
            << (ok ? "OK" : "FAILED") << std::endl;
  return ok;
}

