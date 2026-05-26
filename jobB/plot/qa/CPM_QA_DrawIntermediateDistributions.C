/*
 * Draw 1D QA distributions from the CPM_QA_RunOfflineDiagnostics outputs.
 *
 * Input files are expected to follow the standard QA naming convention:
 *   PREFIX_QA_B0_event_index.root
 *   PREFIX_QA_B1_poca.root
 *   PREFIX_QA_B2_voxel_corrections.root
 *   PREFIX_B3_average_correction_histograms.root
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH3.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace CPMQAPlot
{
  struct TreePlot
  {
    std::string tree_name;
    std::string branch_name;
    std::string hist_name;
    std::string title;
    int bins = 100;
    double min = 0.0;
    double max = 1.0;
    std::string selection;
    bool logy = false;
  };

  struct Hist3Plot
  {
    std::string source_name;
    std::string entries_name;
    std::string hist_name;
    std::string title;
    int bins = 100;
    double min = 0.0;
    double max = 1.0;
    bool require_entries = true;
    bool logy = false;
  };

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

  bool file_exists(const std::string& filename)
  {
    return !filename.empty() && (!gSystem || !gSystem->AccessPathName(filename.c_str()));
  }

  void draw_histogram(
      TH1D& histogram,
      const std::string& plot_dir,
      const bool save_pdf,
      const bool logy)
  {
    TCanvas canvas(
        TString::Format("c_%s", histogram.GetName()),
        histogram.GetTitle(),
        900,
        700);
    if (logy)
    {
      canvas.SetLogy();
    }
    histogram.SetLineWidth(2);
    histogram.Draw("hist");
    if (save_pdf)
    {
      canvas.SaveAs(join_path(plot_dir, std::string(histogram.GetName()) + ".pdf").c_str());
    }
  }

  bool draw_tree_plot(
      TFile& input,
      TFile& output,
      const TreePlot& plot,
      const std::string& plot_dir,
      const bool save_pdf)
  {
    auto* tree = dynamic_cast<TTree*>(input.Get(plot.tree_name.c_str()));
    if (!tree)
    {
      return false;
    }
    if (!tree->GetBranch(plot.branch_name.c_str()))
    {
      std::cout << "CPM_QA_DrawIntermediateDistributions - missing branch "
                << plot.tree_name << "." << plot.branch_name << std::endl;
      return false;
    }

    output.cd();
    TH1D histogram(
        plot.hist_name.c_str(),
        plot.title.c_str(),
        plot.bins,
        plot.min,
        plot.max);
    histogram.Sumw2();
    histogram.SetDirectory(&output);

    const TString draw_expression =
        TString::Format("%s>>%s", plot.branch_name.c_str(), plot.hist_name.c_str());
    tree->Draw(draw_expression, plot.selection.c_str(), "goff");
    histogram.Write();
    draw_histogram(histogram, plot_dir, save_pdf, plot.logy);
    return true;
  }

  bool draw_th3_distribution(
      TFile& input,
      TFile& output,
      const Hist3Plot& plot,
      const std::string& plot_dir,
      const bool save_pdf)
  {
    auto* source = dynamic_cast<TH3*>(input.Get(plot.source_name.c_str()));
    if (!source)
    {
      return false;
    }

    TH3* entries = nullptr;
    if (!plot.entries_name.empty())
    {
      entries = dynamic_cast<TH3*>(input.Get(plot.entries_name.c_str()));
    }

    output.cd();
    TH1D histogram(
        plot.hist_name.c_str(),
        plot.title.c_str(),
        plot.bins,
        plot.min,
        plot.max);
    histogram.Sumw2();
    histogram.SetDirectory(&output);

    for (int ix = 1; ix <= source->GetNbinsX(); ++ix)
    {
      for (int iy = 1; iy <= source->GetNbinsY(); ++iy)
      {
        for (int iz = 1; iz <= source->GetNbinsZ(); ++iz)
        {
          if (plot.require_entries && entries && entries->GetBinContent(ix, iy, iz) <= 0.0)
          {
            continue;
          }
          const double value = source->GetBinContent(ix, iy, iz);
          if (std::isfinite(value))
          {
            histogram.Fill(value);
          }
        }
      }
    }

    histogram.Write();
    draw_histogram(histogram, plot_dir, save_pdf, plot.logy);
    return true;
  }

  unsigned int draw_tree_plots(
      const std::string& filename,
      TFile& output,
      const std::vector<TreePlot>& plots,
      const std::string& plot_dir,
      const bool save_pdf)
  {
    if (!file_exists(filename))
    {
      std::cout << "CPM_QA_DrawIntermediateDistributions - skipping missing file: "
                << filename << std::endl;
      return 0;
    }

    std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
    if (!input || input->IsZombie())
    {
      std::cout << "CPM_QA_DrawIntermediateDistributions - could not open: "
                << filename << std::endl;
      return 0;
    }

    unsigned int drawn = 0;
    for (const auto& plot : plots)
    {
      if (draw_tree_plot(*input, output, plot, plot_dir, save_pdf))
      {
        ++drawn;
      }
    }
    return drawn;
  }

  unsigned int draw_hist3_plots(
      const std::string& filename,
      TFile& output,
      const std::vector<Hist3Plot>& plots,
      const std::string& plot_dir,
      const bool save_pdf)
  {
    if (!file_exists(filename))
    {
      std::cout << "CPM_QA_DrawIntermediateDistributions - skipping missing file: "
                << filename << std::endl;
      return 0;
    }

    std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
    if (!input || input->IsZombie())
    {
      std::cout << "CPM_QA_DrawIntermediateDistributions - could not open: "
                << filename << std::endl;
      return 0;
    }

    unsigned int drawn = 0;
    for (const auto& plot : plots)
    {
      if (draw_th3_distribution(*input, output, plot, plot_dir, save_pdf))
      {
        ++drawn;
      }
    }
    return drawn;
  }
}

bool CPM_QA_DrawIntermediateDistributions(
    const std::string& qa_output_dir = ".",
    const std::string& prefix = "CPM_QA",
    const std::string& plot_output_dir = "",
    const bool save_pdf = true)
{
  if (gStyle)
  {
    gStyle->SetOptStat(1110);
  }

  const std::string output_dir = plot_output_dir.empty() ?
      CPMQAPlot::join_path(qa_output_dir, "qa_plots") :
      plot_output_dir;
  if (gSystem)
  {
    gSystem->mkdir(output_dir.c_str(), true);
  }

  const std::string output_file =
      CPMQAPlot::join_path(output_dir, prefix + "_QA_1d_distributions.root");
  TFile output(output_file.c_str(), "RECREATE");
  if (output.IsZombie())
  {
    std::cout << "CPM_QA_DrawIntermediateDistributions - could not create: "
              << output_file << std::endl;
    return false;
  }

  const std::string b0_file =
      CPMQAPlot::join_path(qa_output_dir, prefix + "_QA_B0_event_index.root");
  const std::string b1_file =
      CPMQAPlot::join_path(qa_output_dir, prefix + "_QA_B1_poca.root");
  const std::string b2_file =
      CPMQAPlot::join_path(qa_output_dir, prefix + "_QA_B2_voxel_corrections.root");
  const std::string b3_file =
      CPMQAPlot::join_path(qa_output_dir, prefix + "_B3_average_correction_histograms.root");

  unsigned int drawn = 0;

  drawn += CPMQAPlot::draw_tree_plots(
      b0_file,
      output,
      {
          {"cpm_event_requests", "object_count", "h_b0_event_object_count", "B0 event object requests;objects/event;events", 100, 0.0, 1000.0, "", true},
          {"cpm_object_requests", "iphi", "h_b0_object_iphi", "B0 requested object iphi;iphi;objects", 36, -0.5, 35.5, "", false},
          {"cpm_object_requests", "ir", "h_b0_object_ir", "B0 requested object ir;ir;objects", 40, -0.5, 39.5, "", false},
          {"cpm_object_requests", "iz", "h_b0_object_iz", "B0 requested object iz;iz;objects", 100, -0.5, 99.5, "", false},
      },
      output_dir,
      save_pdf);

  drawn += CPMQAPlot::draw_tree_plots(
      b1_file,
      output,
      {
          {"cpm_poca_pairs", "dca", "h_b1_pair_dca", "B1 accepted pair DCA;DCA [cm];pairs", 100, 0.0, 5.0, "", true},
          {"cpm_poca_pairs", "pair_weight", "h_b1_pair_weight", "B1 pair weight;weight;pairs", 100, 0.0, 4.0, "", true},
          {"cpm_poca_pairs", "pt_a", "h_b1_pair_pt_a", "B1 positive-track pT;pT [GeV];pairs", 100, 0.0, 10.0, "", true},
          {"cpm_poca_pairs", "pt_b", "h_b1_pair_pt_b", "B1 negative-track pT;pT [GeV];pairs", 100, 0.0, 10.0, "", true},
          {"cpm_poca_pairs", "delta_r", "h_b1_pair_delta_r", "B1 pair #Deltar;#Deltar [cm];pairs", 120, -5.0, 5.0, "", false},
          {"cpm_poca_pairs", "delta_rphi", "h_b1_pair_delta_rphi", "B1 pair r#Delta#phi;r#Delta#phi [cm];pairs", 120, -5.0, 5.0, "", false},
          {"cpm_poca_pairs", "delta_phi", "h_b1_pair_delta_phi", "B1 pair #Delta#phi;#Delta#phi [rad];pairs", 120, -0.05, 0.05, "", false},
          {"cpm_poca_pairs", "delta_z", "h_b1_pair_delta_z", "B1 pair #Deltaz;#Deltaz [cm];pairs", 120, -5.0, 5.0, "", false},
          {"cpm_poca_pairs", "midpoint_z", "h_b1_pair_midpoint_z", "B1 crossing midpoint z;z [cm];pairs", 120, -110.0, 110.0, "", false},
          {"cpm_poca_pairs", "used_line_fallback", "h_b1_pair_line_fallback", "B1 helix line fallback flag;used fallback;pairs", 2, -0.5, 1.5, "", false},
          {"cpm_b1_batch_corrections", "positive_records", "h_b1_batch_positive_records", "B1 batch positive records;records;batch entries", 25, -0.5, 24.5, "", false},
          {"cpm_b1_batch_corrections", "negative_records", "h_b1_batch_negative_records", "B1 batch negative records;records;batch entries", 25, -0.5, 24.5, "", false},
          {"cpm_b1_batch_corrections", "candidate_pairs", "h_b1_batch_candidate_pairs", "B1 batch candidate pairs;pairs;batch entries", 100, 0.0, 1000.0, "", true},
          {"cpm_b1_batch_corrections", "accepted_pairs", "h_b1_batch_accepted_pairs", "B1 batch accepted pairs;pairs;batch entries", 100, 0.0, 1000.0, "", true},
          {"cpm_b1_batch_corrections", "mean_dca", "h_b1_batch_mean_dca", "B1 batch mean DCA;DCA [cm];batch entries", 100, 0.0, 5.0, "", false},
          {"cpm_b1_batch_corrections", "effective_pair_entries", "h_b1_batch_effective_entries", "B1 batch effective pair entries;entries;batch entries", 100, 0.0, 1000.0, "", true},
          {"cpm_b1_batch_corrections", "weighted_mean_delta_r", "h_b1_batch_weighted_mean_delta_r", "B1 batch weighted mean #Deltar;#Deltar [cm];batch entries", 120, -5.0, 5.0, "", false},
          {"cpm_b1_batch_corrections", "weighted_mean_delta_rphi", "h_b1_batch_weighted_mean_delta_rphi", "B1 batch weighted mean r#Delta#phi;r#Delta#phi [cm];batch entries", 120, -5.0, 5.0, "", false},
          {"cpm_b1_batch_corrections", "weighted_mean_delta_phi", "h_b1_batch_weighted_mean_delta_phi", "B1 batch weighted mean #Delta#phi;#Delta#phi [rad];batch entries", 120, -0.05, 0.05, "", false},
          {"cpm_b1_batch_corrections", "weighted_mean_delta_z", "h_b1_batch_weighted_mean_delta_z", "B1 batch weighted mean #Deltaz;#Deltaz [cm];batch entries", 120, -5.0, 5.0, "", false},
          {"cpm_b1_voxel_summary", "records", "h_b1_voxel_records", "B1 voxel input records;records/voxel;voxels", 100, 0.0, 200.0, "", true},
          {"cpm_b1_voxel_summary", "selected_records", "h_b1_voxel_selected_records", "B1 voxel selected records;records/voxel;voxels", 100, 0.0, 100.0, "", true},
          {"cpm_b1_voxel_summary", "positive_records", "h_b1_voxel_positive_records", "B1 voxel positive records;records/voxel;voxels", 100, 0.0, 100.0, "", true},
          {"cpm_b1_voxel_summary", "negative_records", "h_b1_voxel_negative_records", "B1 voxel negative records;records/voxel;voxels", 100, 0.0, 100.0, "", true},
          {"cpm_b1_voxel_summary", "accepted_pairs", "h_b1_voxel_accepted_pairs", "B1 voxel accepted pairs;pairs/voxel;voxels", 100, 0.0, 1000.0, "", true},
          {"cpm_b1_voxel_summary", "batches", "h_b1_voxel_batches", "B1 voxel batches;batches/voxel;voxels", 50, -0.5, 49.5, "", true},
          {"cpm_b1_voxel_summary", "duplicate_dropped_records", "h_b1_voxel_duplicate_dropped_records", "B1 duplicate dropped records;records/voxel;voxels", 100, 0.0, 100.0, "", true},
          {"cpm_b1_voxel_summary", "status", "h_b1_voxel_status", "B1 voxel status;status code;voxels", 10, -0.5, 9.5, "", false},
      },
      output_dir,
      save_pdf);

  drawn += CPMQAPlot::draw_tree_plots(
      b2_file,
      output,
      {
          {"cpm_voxel_corrections", "entries", "h_b2_voxel_entries", "B2 voxel entries;accepted pairs/voxel;voxels", 100, 0.0, 1000.0, "", true},
          {"cpm_voxel_corrections", "effective_pair_entries", "h_b2_voxel_effective_pair_entries", "B2 effective pair entries;effective entries;voxels", 100, 0.0, 1000.0, "", true},
          {"cpm_voxel_corrections", "mean_pair_weight", "h_b2_voxel_mean_pair_weight", "B2 mean pair weight;weight;voxels", 100, 0.0, 4.0, "", true},
          {"cpm_voxel_corrections", "mean_delta_r", "h_b2_voxel_mean_delta_r", "B2 mean #Deltar;#Deltar [cm];voxels", 120, -5.0, 5.0, "", false},
          {"cpm_voxel_corrections", "mean_delta_rphi", "h_b2_voxel_mean_delta_rphi", "B2 mean r#Delta#phi;r#Delta#phi [cm];voxels", 120, -5.0, 5.0, "", false},
          {"cpm_voxel_corrections", "mean_delta_phi", "h_b2_voxel_mean_delta_phi", "B2 mean #Delta#phi;#Delta#phi [rad];voxels", 120, -0.05, 0.05, "", false},
          {"cpm_voxel_corrections", "mean_delta_z", "h_b2_voxel_mean_delta_z", "B2 mean #Deltaz;#Deltaz [cm];voxels", 120, -5.0, 5.0, "", false},
          {"cpm_voxel_corrections", "rms_delta_r", "h_b2_voxel_rms_delta_r", "B2 RMS #Deltar;RMS [cm];voxels", 100, 0.0, 5.0, "", false},
          {"cpm_voxel_corrections", "rms_delta_rphi", "h_b2_voxel_rms_delta_rphi", "B2 RMS r#Delta#phi;RMS [cm];voxels", 100, 0.0, 5.0, "", false},
          {"cpm_voxel_corrections", "rms_delta_z", "h_b2_voxel_rms_delta_z", "B2 RMS #Deltaz;RMS [cm];voxels", 100, 0.0, 5.0, "", false},
          {"cpm_voxel_corrections", "mean_dca", "h_b2_voxel_mean_dca", "B2 mean DCA;DCA [cm];voxels", 100, 0.0, 5.0, "", false},
          {"cpm_voxel_corrections", "rms_dca", "h_b2_voxel_rms_dca", "B2 RMS DCA;DCA [cm];voxels", 100, 0.0, 5.0, "", false},
      },
      output_dir,
      save_pdf);

  drawn += CPMQAPlot::draw_hist3_plots(
      b3_file,
      output,
      {
          {"hentries_rec", "", "h_b3_entries_rec", "B3 filled voxel entries;accepted pairs/voxel;voxels", 100, 0.0, 1000.0, false, true},
          {"hDistortionR_rec", "hentries_rec", "h_b3_distortion_r_rec", "B3 hDistortionR_rec bin contents;#Deltar [cm];voxels", 120, -5.0, 5.0, true, false},
          {"hDistortionP_rec", "hentries_rec", "h_b3_distortion_p_rec", "B3 hDistortionP_rec bin contents;#Delta#phi [rad];voxels", 120, -0.05, 0.05, true, false},
          {"hDistortionZ_rec", "hentries_rec", "h_b3_distortion_z_rec", "B3 hDistortionZ_rec bin contents;#Deltaz [cm];voxels", 120, -5.0, 5.0, true, false},
      },
      output_dir,
      save_pdf);

  output.Close();

  std::cout << "CPM_QA_DrawIntermediateDistributions - output ROOT: "
            << output_file << std::endl;
  std::cout << "CPM_QA_DrawIntermediateDistributions - histograms drawn: "
            << drawn << std::endl;
  if (save_pdf)
  {
    std::cout << "CPM_QA_DrawIntermediateDistributions - PDF directory: "
              << output_dir << std::endl;
  }

  return drawn > 0;
}
