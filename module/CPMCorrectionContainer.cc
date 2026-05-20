#include "CPMCorrectionContainer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

CPMCorrectionContainer::CPMCorrectionContainer()
{
  Reset();
}

void CPMCorrectionContainer::identify(std::ostream& out) const
{
  out << "CPMCorrectionContainer" << std::endl;
  out << "  bins(phi,r,z): " << m_phiBins << ", " << m_rBins << ", " << m_zBins << std::endl;
  out << "  r range: [" << m_rMin << ", " << m_rMax << "]" << std::endl;
  out << "  z range: [" << m_zMin << ", " << m_zMax << "]" << std::endl;
}

void CPMCorrectionContainer::Reset()
{
  const int gridSize = get_grid_size();
  m_entries.assign(gridSize, 0);
  m_sum_weight.assign(gridSize, 0.0);
  m_sum_weight2.assign(gridSize, 0.0);
  m_sum_delta_r.assign(gridSize, 0.0);
  m_sum_delta_r2.assign(gridSize, 0.0);
  m_sum_delta_phi.assign(gridSize, 0.0);
  m_sum_delta_phi2.assign(gridSize, 0.0);
  m_sum_delta_rphi.assign(gridSize, 0.0);
  m_sum_delta_rphi2.assign(gridSize, 0.0);
  m_sum_delta_z.assign(gridSize, 0.0);
  m_sum_delta_z2.assign(gridSize, 0.0);
  m_sum_weighted_delta_r.assign(gridSize, 0.0);
  m_sum_weighted_delta_r2.assign(gridSize, 0.0);
  m_sum_weighted_delta_phi.assign(gridSize, 0.0);
  m_sum_weighted_delta_phi2.assign(gridSize, 0.0);
  m_sum_weighted_delta_rphi.assign(gridSize, 0.0);
  m_sum_weighted_delta_rphi2.assign(gridSize, 0.0);
  m_sum_weighted_delta_z.assign(gridSize, 0.0);
  m_sum_weighted_delta_z2.assign(gridSize, 0.0);
  m_sum_dca.assign(gridSize, 0.0);
  m_sum_dca2.assign(gridSize, 0.0);
}

void CPMCorrectionContainer::set_grid_dimensions(const int phiBins, const int rBins, const int zBins)
{
  if (phiBins <= 0 || rBins <= 0 || zBins <= 0)
  {
    std::cout << "CPMCorrectionContainer::set_grid_dimensions - invalid grid dimensions" << std::endl;
    return;
  }
  m_phiBins = phiBins;
  m_rBins = rBins;
  m_zBins = zBins;
  Reset();
}

void CPMCorrectionContainer::get_grid_dimensions(int& phiBins, int& rBins, int& zBins) const
{
  phiBins = m_phiBins;
  rBins = m_rBins;
  zBins = m_zBins;
}

int CPMCorrectionContainer::get_grid_size() const
{
  return m_phiBins * m_rBins * m_zBins;
}

int CPMCorrectionContainer::get_cell_index(const int iphi, const int ir, const int iz) const
{
  if (iphi < 0 || iphi >= m_phiBins)
  {
    return -1;
  }
  if (ir < 0 || ir >= m_rBins)
  {
    return -1;
  }
  if (iz < 0 || iz >= m_zBins)
  {
    return -1;
  }
  return iz + m_zBins * (ir + m_rBins * iphi);
}

void CPMCorrectionContainer::set_grid_range(
    const double rMin,
    const double rMax,
    const double zMin,
    const double zMax)
{
  m_rMin = rMin;
  m_rMax = rMax;
  m_zMin = zMin;
  m_zMax = zMax;
}

void CPMCorrectionContainer::get_grid_range(double& rMin, double& rMax, double& zMin, double& zMax) const
{
  rMin = m_rMin;
  rMax = m_rMax;
  zMin = m_zMin;
  zMax = m_zMax;
}

int CPMCorrectionContainer::get_entries(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_entries[cellIndex] : 0;
}

double CPMCorrectionContainer::get_sum_weight(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_weight[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_weight2(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_weight2[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_delta_r(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_delta_r[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_delta_r2(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_delta_r2[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_delta_phi(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_delta_phi[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_delta_phi2(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_delta_phi2[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_delta_rphi(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_delta_rphi[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_delta_rphi2(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_delta_rphi2[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_delta_z(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_delta_z[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_delta_z2(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_delta_z2[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_weighted_delta_r(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_weighted_delta_r[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_weighted_delta_r2(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_weighted_delta_r2[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_weighted_delta_phi(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_weighted_delta_phi[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_weighted_delta_phi2(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_weighted_delta_phi2[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_weighted_delta_rphi(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_weighted_delta_rphi[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_weighted_delta_rphi2(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_weighted_delta_rphi2[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_weighted_delta_z(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_weighted_delta_z[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_weighted_delta_z2(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_weighted_delta_z2[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_dca(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_dca[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_sum_dca2(const int cellIndex) const
{
  return bound_check(cellIndex) ? m_sum_dca2[cellIndex] : 0.0;
}

double CPMCorrectionContainer::get_mean_delta_r(const int cellIndex, const bool weighted) const
{
  return weighted ?
      weighted_mean(get_sum_weighted_delta_r(cellIndex), get_sum_weight(cellIndex)) :
      mean(get_sum_delta_r(cellIndex), get_entries(cellIndex));
}

double CPMCorrectionContainer::get_mean_delta_phi(const int cellIndex, const bool weighted) const
{
  return weighted ?
      weighted_mean(get_sum_weighted_delta_phi(cellIndex), get_sum_weight(cellIndex)) :
      mean(get_sum_delta_phi(cellIndex), get_entries(cellIndex));
}

double CPMCorrectionContainer::get_mean_delta_rphi(const int cellIndex, const bool weighted) const
{
  return weighted ?
      weighted_mean(get_sum_weighted_delta_rphi(cellIndex), get_sum_weight(cellIndex)) :
      mean(get_sum_delta_rphi(cellIndex), get_entries(cellIndex));
}

double CPMCorrectionContainer::get_mean_delta_z(const int cellIndex, const bool weighted) const
{
  return weighted ?
      weighted_mean(get_sum_weighted_delta_z(cellIndex), get_sum_weight(cellIndex)) :
      mean(get_sum_delta_z(cellIndex), get_entries(cellIndex));
}

void CPMCorrectionContainer::add_sample(
    const int cellIndex,
    const double deltaR,
    const double deltaPhi,
    const double deltaRPhi,
    const double deltaZ,
    const double dca,
    const double weight)
{
  if (!bound_check(cellIndex) ||
      !std::isfinite(deltaR) ||
      !std::isfinite(deltaPhi) ||
      !std::isfinite(deltaRPhi) ||
      !std::isfinite(deltaZ) ||
      !std::isfinite(dca) ||
      !std::isfinite(weight) ||
      weight <= 0.0)
  {
    return;
  }

  ++m_entries[cellIndex];
  m_sum_weight[cellIndex] += weight;
  m_sum_weight2[cellIndex] += weight * weight;
  m_sum_delta_r[cellIndex] += deltaR;
  m_sum_delta_r2[cellIndex] += deltaR * deltaR;
  m_sum_delta_phi[cellIndex] += deltaPhi;
  m_sum_delta_phi2[cellIndex] += deltaPhi * deltaPhi;
  m_sum_delta_rphi[cellIndex] += deltaRPhi;
  m_sum_delta_rphi2[cellIndex] += deltaRPhi * deltaRPhi;
  m_sum_delta_z[cellIndex] += deltaZ;
  m_sum_delta_z2[cellIndex] += deltaZ * deltaZ;
  m_sum_weighted_delta_r[cellIndex] += weight * deltaR;
  m_sum_weighted_delta_r2[cellIndex] += weight * deltaR * deltaR;
  m_sum_weighted_delta_phi[cellIndex] += weight * deltaPhi;
  m_sum_weighted_delta_phi2[cellIndex] += weight * deltaPhi * deltaPhi;
  m_sum_weighted_delta_rphi[cellIndex] += weight * deltaRPhi;
  m_sum_weighted_delta_rphi2[cellIndex] += weight * deltaRPhi * deltaRPhi;
  m_sum_weighted_delta_z[cellIndex] += weight * deltaZ;
  m_sum_weighted_delta_z2[cellIndex] += weight * deltaZ * deltaZ;
  m_sum_dca[cellIndex] += dca;
  m_sum_dca2[cellIndex] += dca * dca;
}

bool CPMCorrectionContainer::add(const CPMCorrectionContainer& other)
{
  if (!has_same_grid(other))
  {
    std::cout << "CPMCorrectionContainer::add - inconsistent grids" << std::endl;
    return false;
  }

  for (int cellIndex = 0; cellIndex < get_grid_size(); ++cellIndex)
  {
    m_entries[cellIndex] += other.m_entries[cellIndex];
    m_sum_weight[cellIndex] += other.m_sum_weight[cellIndex];
    m_sum_weight2[cellIndex] += other.m_sum_weight2[cellIndex];
    m_sum_delta_r[cellIndex] += other.m_sum_delta_r[cellIndex];
    m_sum_delta_r2[cellIndex] += other.m_sum_delta_r2[cellIndex];
    m_sum_delta_phi[cellIndex] += other.m_sum_delta_phi[cellIndex];
    m_sum_delta_phi2[cellIndex] += other.m_sum_delta_phi2[cellIndex];
    m_sum_delta_rphi[cellIndex] += other.m_sum_delta_rphi[cellIndex];
    m_sum_delta_rphi2[cellIndex] += other.m_sum_delta_rphi2[cellIndex];
    m_sum_delta_z[cellIndex] += other.m_sum_delta_z[cellIndex];
    m_sum_delta_z2[cellIndex] += other.m_sum_delta_z2[cellIndex];
    m_sum_weighted_delta_r[cellIndex] += other.m_sum_weighted_delta_r[cellIndex];
    m_sum_weighted_delta_r2[cellIndex] += other.m_sum_weighted_delta_r2[cellIndex];
    m_sum_weighted_delta_phi[cellIndex] += other.m_sum_weighted_delta_phi[cellIndex];
    m_sum_weighted_delta_phi2[cellIndex] += other.m_sum_weighted_delta_phi2[cellIndex];
    m_sum_weighted_delta_rphi[cellIndex] += other.m_sum_weighted_delta_rphi[cellIndex];
    m_sum_weighted_delta_rphi2[cellIndex] += other.m_sum_weighted_delta_rphi2[cellIndex];
    m_sum_weighted_delta_z[cellIndex] += other.m_sum_weighted_delta_z[cellIndex];
    m_sum_weighted_delta_z2[cellIndex] += other.m_sum_weighted_delta_z2[cellIndex];
    m_sum_dca[cellIndex] += other.m_sum_dca[cellIndex];
    m_sum_dca2[cellIndex] += other.m_sum_dca2[cellIndex];
  }

  return true;
}

bool CPMCorrectionContainer::bound_check(const int cellIndex) const
{
  return cellIndex >= 0 && cellIndex < static_cast<int>(m_entries.size());
}

bool CPMCorrectionContainer::has_same_grid(const CPMCorrectionContainer& other) const
{
  return m_phiBins == other.m_phiBins &&
         m_rBins == other.m_rBins &&
         m_zBins == other.m_zBins &&
         same_value(m_rMin, other.m_rMin) &&
         same_value(m_rMax, other.m_rMax) &&
         same_value(m_zMin, other.m_zMin) &&
         same_value(m_zMax, other.m_zMax);
}

bool CPMCorrectionContainer::same_value(const double lhs, const double rhs)
{
  return std::fabs(lhs - rhs) < 1.0e-9;
}

double CPMCorrectionContainer::mean(const double sum, const int entries)
{
  return entries > 0 ? sum / static_cast<double>(entries) : std::numeric_limits<double>::quiet_NaN();
}

double CPMCorrectionContainer::weighted_mean(const double sumWeighted, const double sumWeight)
{
  return sumWeight > 0.0 ? sumWeighted / sumWeight : std::numeric_limits<double>::quiet_NaN();
}
