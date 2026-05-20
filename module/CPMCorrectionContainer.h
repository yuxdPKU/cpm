#ifndef CPM_CPMCORRECTIONCONTAINER_H
#define CPM_CPMCORRECTIONCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>
#include <vector>

class CPMCorrectionContainer : public PHObject
{
 public:
  CPMCorrectionContainer();
  ~CPMCorrectionContainer() override = default;

  void identify(std::ostream& out = std::cout) const override;
  void Reset() override;

  void set_grid_dimensions(int phiBins, int rBins, int zBins);
  void get_grid_dimensions(int& phiBins, int& rBins, int& zBins) const;
  int get_grid_size() const;
  int get_cell_index(int iphi, int ir, int iz) const;

  void set_grid_range(double rMin, double rMax, double zMin, double zMax);
  void get_grid_range(double& rMin, double& rMax, double& zMin, double& zMax) const;

  int get_entries(int cellIndex) const;
  double get_sum_weight(int cellIndex) const;
  double get_sum_weight2(int cellIndex) const;
  double get_sum_delta_r(int cellIndex) const;
  double get_sum_delta_r2(int cellIndex) const;
  double get_sum_delta_phi(int cellIndex) const;
  double get_sum_delta_phi2(int cellIndex) const;
  double get_sum_delta_rphi(int cellIndex) const;
  double get_sum_delta_rphi2(int cellIndex) const;
  double get_sum_delta_z(int cellIndex) const;
  double get_sum_delta_z2(int cellIndex) const;
  double get_sum_weighted_delta_r(int cellIndex) const;
  double get_sum_weighted_delta_r2(int cellIndex) const;
  double get_sum_weighted_delta_phi(int cellIndex) const;
  double get_sum_weighted_delta_phi2(int cellIndex) const;
  double get_sum_weighted_delta_rphi(int cellIndex) const;
  double get_sum_weighted_delta_rphi2(int cellIndex) const;
  double get_sum_weighted_delta_z(int cellIndex) const;
  double get_sum_weighted_delta_z2(int cellIndex) const;
  double get_sum_dca(int cellIndex) const;
  double get_sum_dca2(int cellIndex) const;

  double get_mean_delta_r(int cellIndex, bool weighted) const;
  double get_mean_delta_phi(int cellIndex, bool weighted) const;
  double get_mean_delta_rphi(int cellIndex, bool weighted) const;
  double get_mean_delta_z(int cellIndex, bool weighted) const;

  void add_sample(
      int cellIndex,
      double deltaR,
      double deltaPhi,
      double deltaRPhi,
      double deltaZ,
      double dca,
      double weight);

  bool add(const CPMCorrectionContainer& other);

 private:
  bool bound_check(int cellIndex) const;
  bool has_same_grid(const CPMCorrectionContainer& other) const;
  static bool same_value(double lhs, double rhs);
  static double mean(double sum, int entries);
  static double weighted_mean(double sumWeighted, double sumWeight);

  int m_phiBins = 36;
  int m_rBins = 16;
  int m_zBins = 80;
  double m_rMin = 20.0;
  double m_rMax = 78.0;
  double m_zMin = -105.5;
  double m_zMax = 105.5;

  std::vector<int> m_entries;
  std::vector<double> m_sum_weight;
  std::vector<double> m_sum_weight2;
  std::vector<double> m_sum_delta_r;
  std::vector<double> m_sum_delta_r2;
  std::vector<double> m_sum_delta_phi;
  std::vector<double> m_sum_delta_phi2;
  std::vector<double> m_sum_delta_rphi;
  std::vector<double> m_sum_delta_rphi2;
  std::vector<double> m_sum_delta_z;
  std::vector<double> m_sum_delta_z2;
  std::vector<double> m_sum_weighted_delta_r;
  std::vector<double> m_sum_weighted_delta_r2;
  std::vector<double> m_sum_weighted_delta_phi;
  std::vector<double> m_sum_weighted_delta_phi2;
  std::vector<double> m_sum_weighted_delta_rphi;
  std::vector<double> m_sum_weighted_delta_rphi2;
  std::vector<double> m_sum_weighted_delta_z;
  std::vector<double> m_sum_weighted_delta_z2;
  std::vector<double> m_sum_dca;
  std::vector<double> m_sum_dca2;

  ClassDefOverride(CPMCorrectionContainer, 1)
};

#endif
