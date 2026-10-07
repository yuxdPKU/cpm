#include "utilities.h"

void prepare_1D_plot(int run = 79516, bool useWeighted = false)
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TString phiType[3] = { "central", "west", "east" };

  const int nfiles = 2;
  oneDHist hists[nfiles]{};
  twoDHist hists_rz[nfiles]{}, hists_pz[nfiles]{}, hists_pr[nfiles]{};
  const TString weightType = useWeighted ? "weighted" : "unweighted";
  const TString outputDir = Form("Rootfiles/%s", weightType.Data());
  if (gSystem->mkdir(outputDir, true) != 0 && gSystem->AccessPathName(outputDir))
  {
    std::cerr << "Failed to create output directory " << outputDir << std::endl;
    return;
  }
  std::cout << "Use CPM " << weightType << " correction maps" << std::endl;
  TString names[nfiles] = 
  {
    "data_reco_FITacts_USEMMSfalse",
    "data_reco_FITacts_USEMMStrue"
  };

  for (int itype = 0; itype < 3; itype++)
  {
    double selectPhi = phi_type2val(phiType[itype]);
    if (std::isnan(selectPhi)) return;

    std::cout << "Use phiType = " << phiType[itype]
          << ", phi = " << selectPhi
          << std::endl;

    // plot distortion vs R, at fixed Z and phi
    std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
    for (const auto& selectZ : selectZs)
    {
      for (int i=0; i<nfiles; i++)
      {
        draw1Dmap_atZ(
          Form("../../output/%s/run%d_CPM_%s_maxdca0p2_B3_average_correction_histograms.root",
               names[i].Data(), run, weightType.Data()),
          Form("CPM_%s", names[i].Data()),
          selectZ, selectPhi, hists[i], kBlack);
      }

      TFile* ofile = new TFile(Form("%s/hist_1Dr_Z%d_%s.root", outputDir.Data(), (int)selectZ, phiType[itype].Data()), "recreate");
      ofile->cd();
      for (int i=0; i<nfiles; i++)
      {
        hists[i].h_N_pos->Write(); hists[i].h_N_neg->Write();
        hists[i].h_R_pos->Write(); hists[i].h_R_neg->Write();
        hists[i].h_P_pos->Write(); hists[i].h_P_neg->Write();
        hists[i].h_Z_pos->Write(); hists[i].h_Z_neg->Write();
      }
      ofile->Write();
      ofile->Close();
      delete ofile;

    }
    
    // plot distortion vs Z, at fixed R and phi
    std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
    for (const auto& selectR : selectRs)
    {
      for (int i=0; i<nfiles; i++)
      {
        draw1Dmap_atR(
          Form("../../output/%s/run%d_CPM_%s_maxdca0p2_B3_average_correction_histograms.root",
               names[i].Data(), run, weightType.Data()),
          Form("CPM_%s", names[i].Data()),
          selectR, selectPhi, hists[i], kBlack);
      }

      TFile* ofile = new TFile(Form("%s/hist_1Dz_R%d_%s.root", outputDir.Data(), (int)selectR, phiType[itype].Data()), "recreate");
      ofile->cd();
      for (int i=0; i<nfiles; i++)
      {
        hists[i].h_N_pos->Write(); hists[i].h_N_neg->Write();
        hists[i].h_R_pos->Write(); hists[i].h_R_neg->Write();
        hists[i].h_P_pos->Write(); hists[i].h_P_neg->Write();
        hists[i].h_Z_pos->Write(); hists[i].h_Z_neg->Write();
      }
      ofile->Write();
      ofile->Close();
      delete ofile;
    }

    // plot distortion at R-Z plane, at fixed phi
    for (int i=0; i<nfiles; i++)
    {

      draw2Dmap_atP(
        Form("../../output/%s/run%d_CPM_%s_maxdca0p2_B3_average_correction_histograms.root",
              names[i].Data(), run, weightType.Data()),
        Form("CPM_%s", names[i].Data()),
        selectPhi, hists_rz[i]);
    }
    TFile* ofile = new TFile(Form("%s/hist_2Drz_%s.root", outputDir.Data(), phiType[itype].Data()),"recreate");
    ofile->cd();
    for (int i=0; i<nfiles; i++)
    {
      hists_rz[i].h_N_pos->Write(); hists_rz[i].h_N_neg->Write();
      hists_rz[i].h_R_pos->Write(); hists_rz[i].h_R_neg->Write();
      hists_rz[i].h_P_pos->Write(); hists_rz[i].h_P_neg->Write();
      hists_rz[i].h_Z_pos->Write(); hists_rz[i].h_Z_neg->Write();
    }
    ofile->Write();
  }

  // plot distortion vs phi, at fixed R and Z
  std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
  std::vector<float> selectZs={5, 10, 15, 30, 60};
  for (const auto& selectR : selectRs)
  {
    for (const auto& selectZ : selectZs)
    {
      for (int i=0; i<nfiles; i++)
      {
        draw1Dmap_atRZ(
          Form("../../output/%s/run%d_CPM_%s_maxdca0p2_B3_average_correction_histograms.root",
              names[i].Data(), run, weightType.Data()),
          Form("CPM_%s", names[i].Data()),
          selectR, selectZ, hists[i], kBlack);
      }

      TFile* ofile = new TFile(Form("%s/hist_1Dp_R%d_Z%d.root",outputDir.Data(),(int)selectR,(int)selectZ),"recreate");
      ofile->cd();
      for (int i=0; i<nfiles; i++)
      {
        hists[i].h_N_pos->Write(); hists[i].h_N_neg->Write();
        hists[i].h_R_pos->Write(); hists[i].h_R_neg->Write();
        hists[i].h_P_pos->Write(); hists[i].h_P_neg->Write();
        hists[i].h_Z_pos->Write(); hists[i].h_Z_neg->Write();
      }
      ofile->Write();
    }
  }

  // plot distortion at phi-R plane, at fixed Z
  for (const auto& selectZ : selectZs)
  {
    for (int i=0; i<nfiles; i++)
    {
      draw2Dmap_atZ(
        Form("../../output/%s/run%d_CPM_%s_maxdca0p2_B3_average_correction_histograms.root",
            names[i].Data(), run, weightType.Data()),
        Form("CPM_%s", names[i].Data()),
        selectZ, hists_pr[i]);
    }

    TFile* ofile = new TFile(Form("%s/hist_2Dpr_Z%d.root",outputDir.Data(),(int)selectZ),"recreate");
    ofile->cd();
    for (int i=0; i<nfiles; i++)
    {
      hists_pr[i].h_N_pos->Write(); hists_pr[i].h_N_neg->Write();
      hists_pr[i].h_R_pos->Write(); hists_pr[i].h_R_neg->Write();
      hists_pr[i].h_P_pos->Write(); hists_pr[i].h_P_neg->Write();
      hists_pr[i].h_Z_pos->Write(); hists_pr[i].h_Z_neg->Write();
    }
    ofile->Write();
  }

  // plot distortion at phi-z plane, at fixed R
  for (const auto& selectR : selectRs)
  {
    for (int i=0; i<nfiles; i++)
    {
      draw2Dmap_atR(
        Form("../../output/%s/run%d_CPM_%s_maxdca0p2_B3_average_correction_histograms.root",
            names[i].Data(), run, weightType.Data()),
        Form("CPM_%s", names[i].Data()),
        selectR, hists_pz[i]);
    }

    TFile* ofile = new TFile(Form("%s/hist_2Dpz_R%d.root",outputDir.Data(),(int)selectR),"recreate");
    ofile->cd();
    for (int i=0; i<nfiles; i++)
    {
      hists_pz[i].h_N_pos->Write(); hists_pz[i].h_N_neg->Write();
      hists_pz[i].h_R_pos->Write(); hists_pz[i].h_R_neg->Write();
      hists_pz[i].h_P_pos->Write(); hists_pz[i].h_P_neg->Write();
      hists_pz[i].h_Z_pos->Write(); hists_pz[i].h_Z_neg->Write();
    }
    ofile->Write();

  }

}
