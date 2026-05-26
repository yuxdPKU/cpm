#ifndef CPM_PHCPMTPCCALIBRATIONQA_H
#define CPM_PHCPMTPCCALIBRATIONQA_H

#include "CPMVoxelContainerv1.h"

#include <string>

class TFile;

class PHCPMTpcCalibrationQA
{
 public:
  static void write_records(
      TFile& output,
      const VoxelContainer& records,
      const std::string& tree_name = "cpm_qa_records");
};

#endif
