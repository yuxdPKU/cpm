/*
 * Draw CPM closure diagnostics for the same fixed-phi/fixed-z R scan used by
 * jobB/plot/draw1D_r_from3D.C.
 *
 * The macro reads B1 cpm_poca_pairs and compares, per R bin,
 *   voxel center - PoCA point from the positive track,
 *   voxel center - PoCA point from the negative track,
 *   voxel center - PoCA midpoint.
 *
 * It also refilters the already-written pair tree with several DCA thresholds
 * to test whether a tighter pair-DCA requirement improves no-distortion
 * closure. It writes only to a closure-specific plot directory and does not
 * touch the original CPM QA/B3 files.
 */

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace CPMClosureSlice
{
  constexpr double kPi = 3.14159265358979323846;

  struct Grid
  {
    int phi_bins = 36;
    int r_bins = 16;
    int z_bins = 80;
    double phi_min = 0.0;
    double phi_max = 2.0 * kPi;
    double r_min = 20.0;
    double r_max = 78.0;
    double z_min = -102.605;
    double z_max = 102.605;
    bool valid = false;
  };

  struct SourceList
  {
    std::vector<std::string> files;
    std::string description;
  };

  struct Branches
  {
    int iphi = -1;
    int ir = -1;
    int iz = -1;
    int charge_a = 0;
    int charge_b = 0;
    double dca = std::numeric_limits<double>::quiet_NaN();
    double point_ax = std::numeric_limits<double>::quiet_NaN();
    double point_ay = std::numeric_limits<double>::quiet_NaN();
    double point_az = std::numeric_limits<double>::quiet_NaN();
    double point_bx = std::numeric_limits<double>::quiet_NaN();
    double point_by = std::numeric_limits<double>::quiet_NaN();
    double point_bz = std::numeric_limits<double>::quiet_NaN();
    double midpoint_x = std::numeric_limits<double>::quiet_NaN();
    double midpoint_y = std::numeric_limits<double>::quiet_NaN();
    double midpoint_z = std::numeric_limits<double>::quiet_NaN();
    double voxel_center_x = std::numeric_limits<double>::quiet_NaN();
    double voxel_center_y = std::numeric_limits<double>::quiet_NaN();
    double voxel_center_z = std::numeric_limits<double>::quiet_NaN();
    double delta_r = std::numeric_limits<double>::quiet_NaN();
    double delta_phi = std::numeric_limits<double>::quiet_NaN();
    double delta_z = std::numeric_limits<double>::quiet_NaN();
  };

  struct Delta
  {
    double r = std::numeric_limits<double>::quiet_NaN();
    double phi = std::numeric_limits<double>::quiet_NaN();
    double z = std::numeric_limits<double>::quiet_NaN();
  };

  struct ThresholdProfiles
  {
    double dca_threshold = 0.0;
    TProfile* point_a_r = nullptr;
    TProfile* point_a_phi = nullptr;
    TProfile* point_a_z = nullptr;
    TProfile* point_b_r = nullptr;
    TProfile* point_b_phi = nullptr;
    TProfile* point_b_z = nullptr;
    TProfile* midpoint_r = nullptr;
    TProfile* midpoint_phi = nullptr;
    TProfile* midpoint_z = nullptr;
    TProfile* dca = nullptr;
    TH1F* entries = nullptr;
    unsigned long long accepted_pairs = 0;
  };

  std::string sanitize_number(const double value)
  {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(2) << value;
    std::string out = stream.str();
    std::replace(out.begin(), out.end(), '-', 'm');
    std::replace(out.begin(), out.end(), '.', 'p');
    return out;
  }

  std::string basename_without_root(const std::string& value)
  {
    const auto slash = value.find_last_of('/');
    std::string out = slash == std::string::npos ? value : value.substr(slash + 1);
    const std::string suffix = ".root";
    if (out.size() >= suffix.size() &&
        out.substr(out.size() - suffix.size()) == suffix)
    {
      out.resize(out.size() - suffix.size());
    }
    return out;
  }

  std::vector<double> parse_dca_thresholds(const std::string& csv)
  {
    std::vector<double> out;
    std::stringstream stream(csv);
    std::string item;
    while (std::getline(stream, item, ','))
    {
      if (item.empty())
      {
        continue;
      }
      const double value = std::stod(item);
      if (std::isfinite(value) && value > 0.0)
      {
        out.push_back(value);
      }
    }
    if (out.empty())
    {
      out = {2.0, 1.0, 0.5, 0.2};
    }
    std::sort(out.begin(), out.end(), std::greater<double>());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  double wrap_delta_phi(double value)
  {
    while (value > kPi)
    {
      value -= 2.0 * kPi;
    }
    while (value <= -kPi)
    {
      value += 2.0 * kPi;
    }
    return value;
  }

  double normalize_phi(double value)
  {
    while (value < 0.0)
    {
      value += 2.0 * kPi;
    }
    while (value >= 2.0 * kPi)
    {
      value -= 2.0 * kPi;
    }
    return value;
  }

  int find_bin_index(const double value, const double min, const double max, const int bins)
  {
    if (bins <= 0 || min >= max || value < min || value >= max)
    {
      return -1;
    }
    int index = static_cast<int>(bins * (value - min) / (max - min));
    if (index < 0)
    {
      index = 0;
    }
    if (index >= bins)
    {
      index = bins - 1;
    }
    return index;
  }

  double bin_center(const int index, const double min, const double max, const int bins)
  {
    return min + (index + 0.5) * (max - min) / bins;
  }

  bool file_exists(const std::string& path)
  {
    return !path.empty() && !gSystem->AccessPathName(path.c_str());
  }

  SourceList read_source_list(
      const std::string& qa_dir,
      const std::string& prefix,
      const std::string& b1_input,
      const bool b1_input_is_list)
  {
    SourceList out;
    if (!b1_input.empty())
    {
      if (b1_input_is_list)
      {
        std::ifstream input(b1_input);
        std::string line;
        while (std::getline(input, line))
        {
          if (!line.empty() && line[0] != '#')
          {
            out.files.push_back(line);
          }
        }
        out.description = b1_input;
      }
      else
      {
        out.files.push_back(b1_input);
        out.description = b1_input;
      }
      return out;
    }

    const std::string chunk_list = qa_dir + "/" + prefix + "_QA_B1_chunks.txt";
    if (file_exists(chunk_list))
    {
      std::ifstream input(chunk_list);
      std::string line;
      while (std::getline(input, line))
      {
        if (!line.empty() && line[0] != '#')
        {
          out.files.push_back(line);
        }
      }
      out.description = chunk_list;
      return out;
    }

    const std::string single_file = qa_dir + "/" + prefix + "_QA_B1_poca.root";
    if (file_exists(single_file))
    {
      out.files.push_back(single_file);
      out.description = single_file;
    }
    return out;
  }

  bool load_grid_from_b1_summary(const std::string& file_name, Grid& grid)
  {
    std::unique_ptr<TFile> input(TFile::Open(file_name.c_str(), "READ"));
    if (!input || input->IsZombie())
    {
      return false;
    }

    auto* summary = dynamic_cast<TTree*>(input->Get("cpm_b1_summary"));
    if (!summary || summary->GetEntries() <= 0)
    {
      return false;
    }

    summary->SetBranchAddress("phi_bins", &grid.phi_bins);
    summary->SetBranchAddress("r_bins", &grid.r_bins);
    summary->SetBranchAddress("z_bins", &grid.z_bins);
    summary->SetBranchAddress("r_min", &grid.r_min);
    summary->SetBranchAddress("r_max", &grid.r_max);
    summary->SetBranchAddress("z_min", &grid.z_min);
    summary->SetBranchAddress("z_max", &grid.z_max);
    summary->GetEntry(0);
    grid.phi_min = 0.0;
    grid.phi_max = 2.0 * kPi;
    grid.valid =
        grid.phi_bins > 0 &&
        grid.r_bins > 0 &&
        grid.z_bins > 0 &&
        grid.r_min < grid.r_max &&
        grid.z_min < grid.z_max;
    return grid.valid;
  }

  bool setup_branches(TChain& chain, Branches& branches)
  {
    const std::vector<std::string> required = {
        "iphi",
        "ir",
        "iz",
        "charge_a",
        "charge_b",
        "dca",
        "point_ax",
        "point_ay",
        "point_az",
        "point_bx",
        "point_by",
        "point_bz",
        "midpoint_x",
        "midpoint_y",
        "midpoint_z",
        "voxel_center_x",
        "voxel_center_y",
        "voxel_center_z",
        "delta_r",
        "delta_phi",
        "delta_z"};
    for (const auto& branch : required)
    {
      if (!chain.GetBranch(branch.c_str()))
      {
        std::cout << "CPM_QA_DrawClosureSlice - missing branch " << branch << std::endl;
        return false;
      }
    }

    chain.SetBranchAddress("iphi", &branches.iphi);
    chain.SetBranchAddress("ir", &branches.ir);
    chain.SetBranchAddress("iz", &branches.iz);
    chain.SetBranchAddress("charge_a", &branches.charge_a);
    chain.SetBranchAddress("charge_b", &branches.charge_b);
    chain.SetBranchAddress("dca", &branches.dca);
    chain.SetBranchAddress("point_ax", &branches.point_ax);
    chain.SetBranchAddress("point_ay", &branches.point_ay);
    chain.SetBranchAddress("point_az", &branches.point_az);
    chain.SetBranchAddress("point_bx", &branches.point_bx);
    chain.SetBranchAddress("point_by", &branches.point_by);
    chain.SetBranchAddress("point_bz", &branches.point_bz);
    chain.SetBranchAddress("midpoint_x", &branches.midpoint_x);
    chain.SetBranchAddress("midpoint_y", &branches.midpoint_y);
    chain.SetBranchAddress("midpoint_z", &branches.midpoint_z);
    chain.SetBranchAddress("voxel_center_x", &branches.voxel_center_x);
    chain.SetBranchAddress("voxel_center_y", &branches.voxel_center_y);
    chain.SetBranchAddress("voxel_center_z", &branches.voxel_center_z);
    chain.SetBranchAddress("delta_r", &branches.delta_r);
    chain.SetBranchAddress("delta_phi", &branches.delta_phi);
    chain.SetBranchAddress("delta_z", &branches.delta_z);
    return true;
  }

  Delta delta_to_voxel(
      const double voxel_x,
      const double voxel_y,
      const double voxel_z,
      const double point_x,
      const double point_y,
      const double point_z)
  {
    Delta out;
    const double voxel_r = std::hypot(voxel_x, voxel_y);
    const double point_r = std::hypot(point_x, point_y);
    const double voxel_phi = std::atan2(voxel_y, voxel_x);
    const double point_phi = std::atan2(point_y, point_x);
    out.r = voxel_r - point_r;
    out.phi = wrap_delta_phi(voxel_phi - point_phi);
    out.z = voxel_z - point_z;
    return out;
  }

  TProfile* make_profile(
      const std::string& name,
      const std::string& title,
      const Grid& grid,
      const double ymin,
      const double ymax)
  {
    auto* profile = new TProfile(
        name.c_str(),
        title.c_str(),
        grid.r_bins,
        grid.r_min,
        grid.r_max);
    profile->SetDirectory(nullptr);
    profile->SetMinimum(ymin);
    profile->SetMaximum(ymax);
    profile->SetStats(false);
    return profile;
  }

  TH1F* make_entries(
      const std::string& name,
      const std::string& title,
      const Grid& grid)
  {
    auto* hist = new TH1F(
        name.c_str(),
        title.c_str(),
        grid.r_bins,
        grid.r_min,
        grid.r_max);
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
    return hist;
  }

  ThresholdProfiles make_threshold_profiles(
      const std::string& prefix,
      const std::string& side_label,
      const double dca_threshold,
      const Grid& grid)
  {
    const std::string dca_tag = sanitize_number(dca_threshold);
    const std::string base = prefix + "_" + side_label + "_dca" + dca_tag;

    ThresholdProfiles out;
    out.dca_threshold = dca_threshold;
    out.point_a_r = make_profile(base + "_point_pos_delta_r", "positive-track PoCA;R [cm];#Deltar [cm]", grid, -1.0, 1.0);
    out.point_a_phi = make_profile(base + "_point_pos_delta_phi", "positive-track PoCA;R [cm];#Delta#phi [rad]", grid, -0.02, 0.02);
    out.point_a_z = make_profile(base + "_point_pos_delta_z", "positive-track PoCA;R [cm];#Deltaz [cm]", grid, -0.5, 0.5);
    out.point_b_r = make_profile(base + "_point_neg_delta_r", "negative-track PoCA;R [cm];#Deltar [cm]", grid, -1.0, 1.0);
    out.point_b_phi = make_profile(base + "_point_neg_delta_phi", "negative-track PoCA;R [cm];#Delta#phi [rad]", grid, -0.02, 0.02);
    out.point_b_z = make_profile(base + "_point_neg_delta_z", "negative-track PoCA;R [cm];#Deltaz [cm]", grid, -0.5, 0.5);
    out.midpoint_r = make_profile(base + "_midpoint_delta_r", "midpoint;R [cm];#Deltar [cm]", grid, -1.0, 1.0);
    out.midpoint_phi = make_profile(base + "_midpoint_delta_phi", "midpoint;R [cm];#Delta#phi [rad]", grid, -0.02, 0.02);
    out.midpoint_z = make_profile(base + "_midpoint_delta_z", "midpoint;R [cm];#Deltaz [cm]", grid, -0.5, 0.5);
    out.dca = make_profile(base + "_mean_dca", "pair DCA;R [cm];DCA [cm]", grid, 0.0, std::max(2.0, dca_threshold));
    out.entries = make_entries(base + "_entries", "accepted pairs;R [cm];pairs", grid);
    return out;
  }

  void set_line_style(TProfile* profile, const int color, const int style, const int width = 2)
  {
    profile->SetLineColor(color);
    profile->SetMarkerColor(color);
    profile->SetLineStyle(style);
    profile->SetLineWidth(width);
    profile->SetMarkerStyle(20);
    profile->SetMarkerSize(0.6);
  }

  void set_hist_style(TH1* hist, const int color, const int style, const int width = 2)
  {
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetLineStyle(style);
    hist->SetLineWidth(width);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.6);
  }

  TH1* extract_b3_slice(
      const std::string& b3_file,
      const std::string& hist_name,
      const double select_phi,
      const double signed_select_z,
      const std::string& output_name)
  {
    if (!file_exists(b3_file))
    {
      return nullptr;
    }

    std::unique_ptr<TFile> input(TFile::Open(b3_file.c_str(), "READ"));
    if (!input || input->IsZombie())
    {
      return nullptr;
    }

    auto* source = dynamic_cast<TH3*>(input->Get(hist_name.c_str()));
    if (!source)
    {
      return nullptr;
    }

    auto* out = new TH1F(
        output_name.c_str(),
        Form("%s;R [cm];value", hist_name.c_str()),
        source->GetYaxis()->GetNbins(),
        source->GetYaxis()->GetXmin(),
        source->GetYaxis()->GetXmax());
    out->SetDirectory(nullptr);

    const int phi_bin = source->GetXaxis()->FindBin(select_phi);
    const int z_bin = source->GetZaxis()->FindBin(signed_select_z);
    for (int ir = 1; ir <= source->GetNbinsY(); ++ir)
    {
      out->SetBinContent(ir, source->GetBinContent(phi_bin, ir, z_bin));
      out->SetBinError(ir, source->GetBinError(phi_bin, ir, z_bin));
    }
    return out;
  }

  void draw_zero_line()
  {
    auto* line = new TLine(gPad->GetUxmin(), 0.0, gPad->GetUxmax(), 0.0);
    line->SetLineColor(kGray + 2);
    line->SetLineStyle(2);
    line->Draw("same");
  }

  void draw_point_components(
      ThresholdProfiles& profiles,
      const std::string& title,
      const std::string& output_pdf)
  {
    set_line_style(profiles.midpoint_r, kBlack, 1, 3);
    set_line_style(profiles.point_a_r, kRed + 1, 1);
    set_line_style(profiles.point_b_r, kBlue + 1, 1);
    set_line_style(profiles.midpoint_phi, kBlack, 1, 3);
    set_line_style(profiles.point_a_phi, kRed + 1, 1);
    set_line_style(profiles.point_b_phi, kBlue + 1, 1);
    set_line_style(profiles.midpoint_z, kBlack, 1, 3);
    set_line_style(profiles.point_a_z, kRed + 1, 1);
    set_line_style(profiles.point_b_z, kBlue + 1, 1);

    TCanvas canvas("c_closure_points", title.c_str(), 1800, 650);
    canvas.Divide(3, 1);

    canvas.cd(1);
    profiles.midpoint_r->SetTitle((title + ";R [cm];voxel - point #Deltar [cm]").c_str());
    profiles.midpoint_r->Draw("hist");
    profiles.point_a_r->Draw("hist same");
    profiles.point_b_r->Draw("hist same");
    draw_zero_line();

    canvas.cd(2);
    profiles.midpoint_phi->SetTitle((title + ";R [cm];voxel - point #Delta#phi [rad]").c_str());
    profiles.midpoint_phi->Draw("hist");
    profiles.point_a_phi->Draw("hist same");
    profiles.point_b_phi->Draw("hist same");
    draw_zero_line();

    canvas.cd(3);
    profiles.midpoint_z->SetTitle((title + ";R [cm];voxel - point #Deltaz [cm]").c_str());
    profiles.midpoint_z->Draw("hist");
    profiles.point_a_z->Draw("hist same");
    profiles.point_b_z->Draw("hist same");
    draw_zero_line();

    auto* legend = new TLegend(0.14, 0.74, 0.88, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(profiles.midpoint_r, "midpoint", "l");
    legend->AddEntry(profiles.point_a_r, "positive-track PoCA point", "l");
    legend->AddEntry(profiles.point_b_r, "negative-track PoCA point", "l");
    canvas.cd(1);
    legend->Draw();

    canvas.SaveAs(output_pdf.c_str());
  }

  void draw_dca_scan_component(
      const std::vector<ThresholdProfiles>& side_profiles,
      const int component,
      TH1* b3_hist,
      const std::string& title,
      const std::string& output_pdf)
  {
    const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 1};

    TCanvas canvas("c_closure_dca_scan", title.c_str(), 900, 700);
    canvas.SetGridx();
    canvas.SetGridy();

    bool drew = false;
    TLegend legend(0.14, 0.72, 0.88, 0.9);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);

    for (std::size_t ithr = 0; ithr < side_profiles.size(); ++ithr)
    {
      TProfile* profile = nullptr;
      if (component == 0)
      {
        profile = side_profiles[ithr].midpoint_r;
      }
      else if (component == 1)
      {
        profile = side_profiles[ithr].midpoint_phi;
      }
      else
      {
        profile = side_profiles[ithr].midpoint_z;
      }
      set_line_style(profile, colors[ithr % colors.size()], 1, ithr == 0 ? 3 : 2);
      profile->SetTitle(title.c_str());
      profile->Draw(drew ? "hist same" : "hist");
      drew = true;
      legend.AddEntry(profile, Form("B1 pairs DCA <= %.2f cm", side_profiles[ithr].dca_threshold), "l");
    }

    if (b3_hist)
    {
      set_hist_style(b3_hist, kGray + 2, 2, 3);
      b3_hist->Draw(drew ? "hist same" : "hist");
      legend.AddEntry(b3_hist, "B3 original map", "l");
    }

    draw_zero_line();
    legend.Draw();
    canvas.SaveAs(output_pdf.c_str());
  }

  void draw_entries_scan(
      const std::vector<ThresholdProfiles>& side_profiles,
      const std::string& title,
      const std::string& output_pdf)
  {
    const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 1};
    TCanvas canvas("c_closure_entries", title.c_str(), 900, 700);
    canvas.SetGridx();
    canvas.SetGridy();
    canvas.SetLogy();

    bool drew = false;
    TLegend legend(0.14, 0.72, 0.88, 0.9);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);

    for (std::size_t ithr = 0; ithr < side_profiles.size(); ++ithr)
    {
      auto* hist = side_profiles[ithr].entries;
      set_hist_style(hist, colors[ithr % colors.size()], 1, ithr == 0 ? 3 : 2);
      hist->SetTitle(title.c_str());
      hist->Draw(drew ? "hist same" : "hist");
      drew = true;
      legend.AddEntry(hist, Form("DCA <= %.2f cm", side_profiles[ithr].dca_threshold), "l");
    }

    legend.Draw();
    canvas.SaveAs(output_pdf.c_str());
  }

  void write_profiles(
      TFile& output,
      std::vector<ThresholdProfiles>& profiles)
  {
    output.cd();
    for (auto& profile : profiles)
    {
      profile.point_a_r->Write();
      profile.point_a_phi->Write();
      profile.point_a_z->Write();
      profile.point_b_r->Write();
      profile.point_b_phi->Write();
      profile.point_b_z->Write();
      profile.midpoint_r->Write();
      profile.midpoint_phi->Write();
      profile.midpoint_z->Write();
      profile.dca->Write();
      profile.entries->Write();
    }
  }
}

bool CPM_QA_DrawClosureSlice(
    const std::string& qa_dir = "output/sim_genfit_unweighted_qa",
    const std::string& prefix = "sim_genfit_unweighted",
    const std::string& plot_dir = "",
    const double select_phi = 4.7,
    const double select_z = 10.0,
    const std::string& dca_thresholds_csv = "2.0,1.0,0.5,0.2",
    const std::string& b1_input = "",
    const bool b1_input_is_list = false,
    const std::string& b3_file = "")
{
  using namespace CPMClosureSlice;

  gStyle->SetOptStat(0);

  const SourceList sources = read_source_list(qa_dir, prefix, b1_input, b1_input_is_list);
  if (sources.files.empty())
  {
    std::cout << "CPM_QA_DrawClosureSlice - no B1 pair source found for "
              << qa_dir << " prefix " << prefix << std::endl;
    return false;
  }

  Grid grid;
  if (!load_grid_from_b1_summary(sources.files.front(), grid))
  {
    std::cout << "CPM_QA_DrawClosureSlice - could not read cpm_b1_summary grid from "
              << sources.files.front() << "; using default grid" << std::endl;
    grid.valid = true;
  }

  const double normalized_phi = normalize_phi(select_phi);
  const int target_iphi = find_bin_index(normalized_phi, grid.phi_min, grid.phi_max, grid.phi_bins);
  const int target_iz_pos = find_bin_index(std::abs(select_z), grid.z_min, grid.z_max, grid.z_bins);
  const int target_iz_neg = find_bin_index(-std::abs(select_z), grid.z_min, grid.z_max, grid.z_bins);
  if (target_iphi < 0 || target_iz_pos < 0 || target_iz_neg < 0)
  {
    std::cout << "CPM_QA_DrawClosureSlice - requested phi/z is outside the grid" << std::endl;
    return false;
  }

  const double phi_center = bin_center(target_iphi, grid.phi_min, grid.phi_max, grid.phi_bins);
  const double z_center_pos = bin_center(target_iz_pos, grid.z_min, grid.z_max, grid.z_bins);
  const double z_center_neg = bin_center(target_iz_neg, grid.z_min, grid.z_max, grid.z_bins);

  TChain pairs("cpm_poca_pairs");
  for (const auto& file : sources.files)
  {
    pairs.Add(file.c_str());
  }
  if (pairs.GetEntries() <= 0)
  {
    std::cout << "CPM_QA_DrawClosureSlice - cpm_poca_pairs is empty. "
              << "Run QA with --write-pair-tree first." << std::endl;
    return false;
  }

  Branches branches;
  if (!setup_branches(pairs, branches))
  {
    return false;
  }

  const std::vector<double> dca_thresholds = parse_dca_thresholds(dca_thresholds_csv);
  std::vector<ThresholdProfiles> pos_profiles;
  std::vector<ThresholdProfiles> neg_profiles;
  for (const double threshold : dca_thresholds)
  {
    pos_profiles.push_back(make_threshold_profiles(prefix, "posz", threshold, grid));
    neg_profiles.push_back(make_threshold_profiles(prefix, "negz", threshold, grid));
  }

  const Long64_t entries = pairs.GetEntries();
  unsigned long long selected_phi_z_pairs = 0;
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    pairs.GetEntry(entry);
    if (branches.iphi != target_iphi)
    {
      continue;
    }

    std::vector<ThresholdProfiles>* target_profiles = nullptr;
    if (branches.iz == target_iz_pos)
    {
      target_profiles = &pos_profiles;
    }
    else if (branches.iz == target_iz_neg)
    {
      target_profiles = &neg_profiles;
    }
    else
    {
      continue;
    }

    if (!std::isfinite(branches.dca))
    {
      continue;
    }

    const double voxel_r = std::hypot(branches.voxel_center_x, branches.voxel_center_y);
    const Delta point_a = delta_to_voxel(
        branches.voxel_center_x,
        branches.voxel_center_y,
        branches.voxel_center_z,
        branches.point_ax,
        branches.point_ay,
        branches.point_az);
    const Delta point_b = delta_to_voxel(
        branches.voxel_center_x,
        branches.voxel_center_y,
        branches.voxel_center_z,
        branches.point_bx,
        branches.point_by,
        branches.point_bz);
    const Delta midpoint = delta_to_voxel(
        branches.voxel_center_x,
        branches.voxel_center_y,
        branches.voxel_center_z,
        branches.midpoint_x,
        branches.midpoint_y,
        branches.midpoint_z);

    ++selected_phi_z_pairs;
    for (auto& profile : *target_profiles)
    {
      if (branches.dca > profile.dca_threshold)
      {
        continue;
      }

      profile.point_a_r->Fill(voxel_r, point_a.r);
      profile.point_a_phi->Fill(voxel_r, point_a.phi);
      profile.point_a_z->Fill(voxel_r, point_a.z);
      profile.point_b_r->Fill(voxel_r, point_b.r);
      profile.point_b_phi->Fill(voxel_r, point_b.phi);
      profile.point_b_z->Fill(voxel_r, point_b.z);
      profile.midpoint_r->Fill(voxel_r, midpoint.r);
      profile.midpoint_phi->Fill(voxel_r, midpoint.phi);
      profile.midpoint_z->Fill(voxel_r, midpoint.z);
      profile.dca->Fill(voxel_r, branches.dca);
      profile.entries->Fill(voxel_r, 1.0);
      ++profile.accepted_pairs;
    }
  }

  const std::string resolved_plot_dir =
      plot_dir.empty() ?
      qa_dir + "/closure_slice_phi" + sanitize_number(normalized_phi) + "_z" + sanitize_number(std::abs(select_z)) :
      plot_dir;
  gSystem->mkdir(resolved_plot_dir.c_str(), true);

  const std::string resolved_b3_file =
      b3_file.empty() ?
      qa_dir + "/" + prefix + "_B3_average_correction_histograms.root" :
      b3_file;

  std::unique_ptr<TH1> b3_pos_r(extract_b3_slice(
      resolved_b3_file,
      "hIntDistortionR_posz",
      normalized_phi,
      std::abs(select_z),
      prefix + "_b3_posz_delta_r"));
  std::unique_ptr<TH1> b3_pos_phi(extract_b3_slice(
      resolved_b3_file,
      "hIntDistortionP_posz",
      normalized_phi,
      std::abs(select_z),
      prefix + "_b3_posz_delta_phi"));
  std::unique_ptr<TH1> b3_pos_z(extract_b3_slice(
      resolved_b3_file,
      "hIntDistortionZ_posz",
      normalized_phi,
      std::abs(select_z),
      prefix + "_b3_posz_delta_z"));
  std::unique_ptr<TH1> b3_neg_r(extract_b3_slice(
      resolved_b3_file,
      "hIntDistortionR_negz",
      normalized_phi,
      -std::abs(select_z),
      prefix + "_b3_negz_delta_r"));
  std::unique_ptr<TH1> b3_neg_phi(extract_b3_slice(
      resolved_b3_file,
      "hIntDistortionP_negz",
      normalized_phi,
      -std::abs(select_z),
      prefix + "_b3_negz_delta_phi"));
  std::unique_ptr<TH1> b3_neg_z(extract_b3_slice(
      resolved_b3_file,
      "hIntDistortionZ_negz",
      normalized_phi,
      -std::abs(select_z),
      prefix + "_b3_negz_delta_z"));

  const std::string root_output = resolved_plot_dir + "/" + prefix + "_closure_slice_profiles.root";
  TFile output(root_output.c_str(), "RECREATE");
  if (output.IsZombie())
  {
    std::cout << "CPM_QA_DrawClosureSlice - could not create " << root_output << std::endl;
    return false;
  }

  write_profiles(output, pos_profiles);
  write_profiles(output, neg_profiles);
  if (b3_pos_r)
  {
    b3_pos_r->Write();
    b3_pos_phi->Write();
    b3_pos_z->Write();
    b3_neg_r->Write();
    b3_neg_phi->Write();
    b3_neg_z->Write();
  }

  TTree summary("cpm_closure_slice_summary", "CPM closure slice plotting summary");
  std::string source_description = sources.description;
  std::string output_plot_dir = resolved_plot_dir;
  std::string output_root_file = root_output;
  std::string b3_source = resolved_b3_file;
  int source_files = static_cast<int>(sources.files.size());
  int target_phi_bin = target_iphi;
  int target_pos_z_bin = target_iz_pos;
  int target_neg_z_bin = target_iz_neg;
  double requested_phi = normalized_phi;
  double requested_z = std::abs(select_z);
  double selected_phi_center = phi_center;
  double selected_pos_z_center = z_center_pos;
  double selected_neg_z_center = z_center_neg;
  unsigned long long scanned_pair_rows = static_cast<unsigned long long>(entries);
  unsigned long long selected_pair_rows = selected_phi_z_pairs;
  summary.Branch("source_description", &source_description);
  summary.Branch("source_files", &source_files);
  summary.Branch("b3_source", &b3_source);
  summary.Branch("plot_dir", &output_plot_dir);
  summary.Branch("root_file", &output_root_file);
  summary.Branch("requested_phi", &requested_phi);
  summary.Branch("requested_z", &requested_z);
  summary.Branch("target_phi_bin", &target_phi_bin);
  summary.Branch("target_pos_z_bin", &target_pos_z_bin);
  summary.Branch("target_neg_z_bin", &target_neg_z_bin);
  summary.Branch("selected_phi_center", &selected_phi_center);
  summary.Branch("selected_pos_z_center", &selected_pos_z_center);
  summary.Branch("selected_neg_z_center", &selected_neg_z_center);
  summary.Branch("scanned_pair_rows", &scanned_pair_rows);
  summary.Branch("selected_pair_rows", &selected_pair_rows);
  summary.Fill();
  summary.Write();
  output.Close();

  const std::string slice_label =
      Form("#phi bin center %.3f rad, z = +%.2f / %.2f cm",
           selected_phi_center,
           selected_pos_z_center,
           selected_neg_z_center);

  draw_point_components(
      pos_profiles.front(),
      prefix + " +Z " + slice_label + Form(", DCA <= %.2f cm", pos_profiles.front().dca_threshold),
      resolved_plot_dir + "/" + prefix + "_closure_points_posz_dca" + sanitize_number(pos_profiles.front().dca_threshold) + ".pdf");
  draw_point_components(
      neg_profiles.front(),
      prefix + " -Z " + slice_label + Form(", DCA <= %.2f cm", neg_profiles.front().dca_threshold),
      resolved_plot_dir + "/" + prefix + "_closure_points_negz_dca" + sanitize_number(neg_profiles.front().dca_threshold) + ".pdf");

  draw_dca_scan_component(
      pos_profiles,
      0,
      b3_pos_r.get(),
      prefix + " +Z midpoint #Deltar DCA scan;R [cm];#Deltar [cm]",
      resolved_plot_dir + "/" + prefix + "_closure_midpoint_delta_r_posz_dca_scan.pdf");
  draw_dca_scan_component(
      pos_profiles,
      1,
      b3_pos_phi.get(),
      prefix + " +Z midpoint #Delta#phi DCA scan;R [cm];#Delta#phi [rad]",
      resolved_plot_dir + "/" + prefix + "_closure_midpoint_delta_phi_posz_dca_scan.pdf");
  draw_dca_scan_component(
      pos_profiles,
      2,
      b3_pos_z.get(),
      prefix + " +Z midpoint #Deltaz DCA scan;R [cm];#Deltaz [cm]",
      resolved_plot_dir + "/" + prefix + "_closure_midpoint_delta_z_posz_dca_scan.pdf");
  draw_dca_scan_component(
      neg_profiles,
      0,
      b3_neg_r.get(),
      prefix + " -Z midpoint #Deltar DCA scan;R [cm];#Deltar [cm]",
      resolved_plot_dir + "/" + prefix + "_closure_midpoint_delta_r_negz_dca_scan.pdf");
  draw_dca_scan_component(
      neg_profiles,
      1,
      b3_neg_phi.get(),
      prefix + " -Z midpoint #Delta#phi DCA scan;R [cm];#Delta#phi [rad]",
      resolved_plot_dir + "/" + prefix + "_closure_midpoint_delta_phi_negz_dca_scan.pdf");
  draw_dca_scan_component(
      neg_profiles,
      2,
      b3_neg_z.get(),
      prefix + " -Z midpoint #Deltaz DCA scan;R [cm];#Deltaz [cm]",
      resolved_plot_dir + "/" + prefix + "_closure_midpoint_delta_z_negz_dca_scan.pdf");

  draw_entries_scan(
      pos_profiles,
      prefix + " +Z accepted pairs after DCA refilter;R [cm];pairs",
      resolved_plot_dir + "/" + prefix + "_closure_entries_posz_dca_scan.pdf");
  draw_entries_scan(
      neg_profiles,
      prefix + " -Z accepted pairs after DCA refilter;R [cm];pairs",
      resolved_plot_dir + "/" + prefix + "_closure_entries_negz_dca_scan.pdf");

  std::cout << "CPM_QA_DrawClosureSlice - B1 source: " << sources.description << std::endl;
  std::cout << "CPM_QA_DrawClosureSlice - source files: " << sources.files.size() << std::endl;
  std::cout << "CPM_QA_DrawClosureSlice - scanned pair rows: " << entries << std::endl;
  std::cout << "CPM_QA_DrawClosureSlice - selected phi/z pair rows before DCA scan: "
            << selected_phi_z_pairs << std::endl;
  std::cout << "CPM_QA_DrawClosureSlice - selected phi bin: " << target_iphi
            << " center " << selected_phi_center << std::endl;
  std::cout << "CPM_QA_DrawClosureSlice - selected z bins: + " << target_iz_pos
            << " center " << selected_pos_z_center
            << " / - " << target_iz_neg
            << " center " << selected_neg_z_center << std::endl;
  for (const auto& profile : pos_profiles)
  {
    std::cout << "CPM_QA_DrawClosureSlice - +Z DCA <= "
              << profile.dca_threshold
              << " accepted pairs: " << profile.accepted_pairs << std::endl;
  }
  for (const auto& profile : neg_profiles)
  {
    std::cout << "CPM_QA_DrawClosureSlice - -Z DCA <= "
              << profile.dca_threshold
              << " accepted pairs: " << profile.accepted_pairs << std::endl;
  }
  std::cout << "CPM_QA_DrawClosureSlice - output directory: " << resolved_plot_dir << std::endl;
  std::cout << "CPM_QA_DrawClosureSlice - profile root: " << root_output << std::endl;
  return true;
}
