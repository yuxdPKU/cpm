#ifndef CPM_CPMCORRECTIONACCUMULATOR_H
#define CPM_CPMCORRECTIONACCUMULATOR_H

#include <algorithm>
#include <cmath>
#include <limits>

struct CPMCorrectionAccumulator
{
  unsigned long long entries = 0;
  double sum_pair_weight = 0.0;
  double sum_pair_weight2 = 0.0;
  double sum_delta_r = 0.0;
  double sum_delta_r2 = 0.0;
  double sum_delta_rphi = 0.0;
  double sum_delta_rphi2 = 0.0;
  double sum_delta_phi = 0.0;
  double sum_delta_phi2 = 0.0;
  double sum_delta_z = 0.0;
  double sum_delta_z2 = 0.0;
  double sum_weighted_delta_r = 0.0;
  double sum_weighted_delta_r2 = 0.0;
  double sum_weighted_delta_rphi = 0.0;
  double sum_weighted_delta_rphi2 = 0.0;
  double sum_weighted_delta_phi = 0.0;
  double sum_weighted_delta_phi2 = 0.0;
  double sum_weighted_delta_z = 0.0;
  double sum_weighted_delta_z2 = 0.0;
  double sum_dca = 0.0;
  double sum_dca2 = 0.0;
  double sum_voxel_x = 0.0;
  double sum_voxel_y = 0.0;
  double sum_voxel_z = 0.0;

  void add(
      double delta_r,
      double delta_rphi,
      double delta_phi,
      double delta_z,
      double dca,
      double weight,
      double voxel_x,
      double voxel_y,
      double voxel_z)
  {
    ++entries;
    sum_pair_weight += weight;
    sum_pair_weight2 += weight * weight;
    sum_delta_r += delta_r;
    sum_delta_r2 += delta_r * delta_r;
    sum_delta_rphi += delta_rphi;
    sum_delta_rphi2 += delta_rphi * delta_rphi;
    sum_delta_phi += delta_phi;
    sum_delta_phi2 += delta_phi * delta_phi;
    sum_delta_z += delta_z;
    sum_delta_z2 += delta_z * delta_z;
    sum_weighted_delta_r += weight * delta_r;
    sum_weighted_delta_r2 += weight * delta_r * delta_r;
    sum_weighted_delta_rphi += weight * delta_rphi;
    sum_weighted_delta_rphi2 += weight * delta_rphi * delta_rphi;
    sum_weighted_delta_phi += weight * delta_phi;
    sum_weighted_delta_phi2 += weight * delta_phi * delta_phi;
    sum_weighted_delta_z += weight * delta_z;
    sum_weighted_delta_z2 += weight * delta_z * delta_z;
    sum_dca += dca;
    sum_dca2 += dca * dca;
    sum_voxel_x += voxel_x;
    sum_voxel_y += voxel_y;
    sum_voxel_z += voxel_z;
  }

  void add_sums(
      unsigned long long added_entries,
      double added_sum_pair_weight,
      double added_sum_pair_weight2,
      double added_sum_delta_r,
      double added_sum_delta_r2,
      double added_sum_delta_rphi,
      double added_sum_delta_rphi2,
      double added_sum_delta_phi,
      double added_sum_delta_phi2,
      double added_sum_delta_z,
      double added_sum_delta_z2,
      double added_sum_dca,
      double added_sum_dca2,
      double voxel_x,
      double voxel_y,
      double voxel_z)
  {
    entries += added_entries;
    sum_pair_weight += added_sum_pair_weight;
    sum_pair_weight2 += added_sum_pair_weight2;
    sum_weighted_delta_r += added_sum_delta_r;
    sum_weighted_delta_r2 += added_sum_delta_r2;
    sum_weighted_delta_rphi += added_sum_delta_rphi;
    sum_weighted_delta_rphi2 += added_sum_delta_rphi2;
    sum_weighted_delta_phi += added_sum_delta_phi;
    sum_weighted_delta_phi2 += added_sum_delta_phi2;
    sum_weighted_delta_z += added_sum_delta_z;
    sum_weighted_delta_z2 += added_sum_delta_z2;
    sum_dca += added_sum_dca;
    sum_dca2 += added_sum_dca2;
    sum_voxel_x += voxel_x * static_cast<double>(added_entries);
    sum_voxel_y += voxel_y * static_cast<double>(added_entries);
    sum_voxel_z += voxel_z * static_cast<double>(added_entries);
  }
};

inline double cpmCorrectionMean(
    const double sum,
    const unsigned long long entries)
{
  return entries > 0 ? sum / static_cast<double>(entries) : std::numeric_limits<double>::quiet_NaN();
}

inline double cpmCorrectionWeightedMean(
    const double sum_weighted,
    const double sum_weight)
{
  return sum_weight > 0.0 ? sum_weighted / sum_weight : std::numeric_limits<double>::quiet_NaN();
}

inline double cpmCorrectionRms(
    const double sum,
    const double sum2,
    const unsigned long long entries)
{
  if (entries == 0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double average = cpmCorrectionMean(sum, entries);
  const double variance = sum2 / static_cast<double>(entries) - average * average;
  return std::sqrt(std::max(0.0, variance));
}

inline double cpmCorrectionWeightedRms(
    const double sum_weighted,
    const double sum_weighted2,
    const double sum_weight)
{
  if (sum_weight <= 0.0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double average = cpmCorrectionWeightedMean(sum_weighted, sum_weight);
  const double variance = sum_weighted2 / sum_weight - average * average;
  return std::sqrt(std::max(0.0, variance));
}

inline double cpmCorrectionEffectiveEntries(
    const double sum_weight,
    const double sum_weight2)
{
  return sum_weight2 > 0.0 ?
      sum_weight * sum_weight / sum_weight2 :
      std::numeric_limits<double>::quiet_NaN();
}

#endif
