#include "utilities.h"

void draw1D_z_fitter(bool useWeighted = false)
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
  const TString weightType = useWeighted ? "weighted" : "unweighted";
  const TString inputDir = Form("Rootfiles/%s", weightType.Data());
  const TString figureDir = Form("figure/%s", weightType.Data());
  gSystem->mkdir(figureDir, true);

  int run = 29;
  std::vector<TString> phiTypes = {"central", "west", "east"};
  for (const auto& phiType : phiTypes)
  {
    std::pair<double, double> thetarange_NCO = {0,0};
    std::pair<double, double> thetarange_NCI = {0,0};
    std::pair<double, double> thetarange_SCI = {0,0};
    std::pair<double, double> thetarange_SCO = {0,0};
    bool drawOuterThetaRange = false;

    if (phiType.CompareTo("central")==0)
    {
      thetarange_SCO = {-0.784379,-0.608235};
      thetarange_SCI = {-0.56699,-0.0301199};
      thetarange_NCI = {0.0264637,0.56564};
      thetarange_NCO = {0.606308,0.915119};
      drawOuterThetaRange = true;
    }
    else if (phiType.CompareTo("east")==0)
    {
      thetarange_SCI = {-0.637513,-0.134286};
      thetarange_NCI = {0.135323,0.639332};
    }
    else if (phiType.CompareTo("west")==0)
    {
      thetarange_SCI = {-0.642148,-0.138535};
      thetarange_NCI = {0.123448,0.633285};
    }

    oneDHist hists_acts{};
    oneDHist hists_genfit{};
    oneDHist hists_truth{};

    std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
    for (const auto& selectR : selectRs)
    {
      get1Dhist(Form("%s/hist_1Dz_R%d_%s.root", inputDir.Data(), (int)selectR, phiType.Data()), "CPM_sim_reco_FITacts_EXTRAbidirectional_CLUSERRraw_DISTORTIONINPUTfalse", hists_acts, kRed-7);
      get1Dhist(Form("%s/hist_1Dz_R%d_%s.root", inputDir.Data(), (int)selectR, phiType.Data()), "CPM_sim_reco_FITgenfit_DISTORTIONINPUTfalse", hists_genfit, kGreen+3);
      get1Dhist(Form("%s/hist_1Dz_R%d_%s.root", inputDir.Data(), (int)selectR, phiType.Data()), "CPM_sim_reco_FITtruth_DISTORTIONINPUTfalse", hists_truth, kAzure-2);

      if (!hists_acts.h_N_pos || !hists_genfit.h_N_pos || !hists_truth.h_N_pos) {
        std::cerr << "Error: some histograms are null for Z = " << selectR << std::endl;
        continue;
      }

      std::vector<oneDHist> allHists = {
        hists_acts,
        hists_genfit,
        hists_truth
      };

      std::vector<TH1*> hists_P; hists_P.clear();
      std::vector<TH1*> hists_R; hists_R.clear();
      std::vector<TH1*> hists_Z; hists_Z.clear();

      AddOneDHistVectorToVectors(allHists, hists_P, hists_R, hists_Z);

      std::pair<double,double> yrange_P = SetCommonYRange(hists_P);
      std::pair<double,double> yrange_R = SetCommonYRange(hists_R);
      std::pair<double,double> yrange_Z = SetCommonYRange(hists_Z);

      TLine *l_nco_o_R = nullptr, *l_nco_i_R = nullptr, *l_nci_o_R = nullptr, *l_nci_i_R = nullptr;
      TLine *l_nco_o_P = nullptr, *l_nco_i_P = nullptr, *l_nci_o_P = nullptr, *l_nci_i_P = nullptr;
      TLine *l_nco_o_Z = nullptr, *l_nco_i_Z = nullptr, *l_nci_o_Z = nullptr, *l_nci_i_Z = nullptr;
      TLine *l_sci_i_R = nullptr, *l_sci_o_R = nullptr, *l_sco_i_R = nullptr, *l_sco_o_R = nullptr;
      TLine *l_sci_i_P = nullptr, *l_sci_o_P = nullptr, *l_sco_i_P = nullptr, *l_sco_o_P = nullptr;
      TLine *l_sci_i_Z = nullptr, *l_sci_o_Z = nullptr, *l_sco_i_Z = nullptr, *l_sco_o_Z = nullptr;

      if (drawOuterThetaRange)
      {
        l_nco_o_R = new TLine(selectR*tan(thetarange_NCO.second),yrange_R.first,selectR*tan(thetarange_NCO.second),yrange_R.second);
        l_nco_o_P = new TLine(selectR*tan(thetarange_NCO.second),yrange_P.first,selectR*tan(thetarange_NCO.second),yrange_P.second);
        l_nco_o_Z = new TLine(selectR*tan(thetarange_NCO.second),yrange_Z.first,selectR*tan(thetarange_NCO.second),yrange_Z.second);

        l_nco_i_R = new TLine(selectR*tan(thetarange_NCO.first),yrange_R.first,selectR*tan(thetarange_NCO.first),yrange_R.second);
        l_nco_i_P = new TLine(selectR*tan(thetarange_NCO.first),yrange_P.first,selectR*tan(thetarange_NCO.first),yrange_P.second);
        l_nco_i_Z = new TLine(selectR*tan(thetarange_NCO.first),yrange_Z.first,selectR*tan(thetarange_NCO.first),yrange_Z.second);

        l_sco_i_R = new TLine(selectR*tan(thetarange_SCO.second),yrange_R.first,selectR*tan(thetarange_SCO.second),yrange_R.second);
        l_sco_i_P = new TLine(selectR*tan(thetarange_SCO.second),yrange_P.first,selectR*tan(thetarange_SCO.second),yrange_P.second);
        l_sco_i_Z = new TLine(selectR*tan(thetarange_SCO.second),yrange_Z.first,selectR*tan(thetarange_SCO.second),yrange_Z.second);

        l_sco_o_R = new TLine(selectR*tan(thetarange_SCO.first),yrange_R.first,selectR*tan(thetarange_SCO.first),yrange_R.second);
        l_sco_o_P = new TLine(selectR*tan(thetarange_SCO.first),yrange_P.first,selectR*tan(thetarange_SCO.first),yrange_P.second);
        l_sco_o_Z = new TLine(selectR*tan(thetarange_SCO.first),yrange_Z.first,selectR*tan(thetarange_SCO.first),yrange_Z.second);
      }

      l_nci_o_R = new TLine(selectR*tan(thetarange_NCI.second),yrange_R.first,selectR*tan(thetarange_NCI.second),yrange_R.second);
      l_nci_o_P = new TLine(selectR*tan(thetarange_NCI.second),yrange_P.first,selectR*tan(thetarange_NCI.second),yrange_P.second);
      l_nci_o_Z = new TLine(selectR*tan(thetarange_NCI.second),yrange_Z.first,selectR*tan(thetarange_NCI.second),yrange_Z.second);

      l_nci_i_R = new TLine(selectR*tan(thetarange_NCI.first),yrange_R.first,selectR*tan(thetarange_NCI.first),yrange_R.second);
      l_nci_i_P = new TLine(selectR*tan(thetarange_NCI.first),yrange_P.first,selectR*tan(thetarange_NCI.first),yrange_P.second);
      l_nci_i_Z = new TLine(selectR*tan(thetarange_NCI.first),yrange_Z.first,selectR*tan(thetarange_NCI.first),yrange_Z.second);

      l_sci_i_R = new TLine(selectR*tan(thetarange_SCI.second),yrange_R.first,selectR*tan(thetarange_SCI.second),yrange_R.second);
      l_sci_i_P = new TLine(selectR*tan(thetarange_SCI.second),yrange_P.first,selectR*tan(thetarange_SCI.second),yrange_P.second);
      l_sci_i_Z = new TLine(selectR*tan(thetarange_SCI.second),yrange_Z.first,selectR*tan(thetarange_SCI.second),yrange_Z.second);

      l_sci_o_R = new TLine(selectR*tan(thetarange_SCI.first),yrange_R.first,selectR*tan(thetarange_SCI.first),yrange_R.second);
      l_sci_o_P = new TLine(selectR*tan(thetarange_SCI.first),yrange_P.first,selectR*tan(thetarange_SCI.first),yrange_P.second);
      l_sci_o_Z = new TLine(selectR*tan(thetarange_SCI.first),yrange_Z.first,selectR*tan(thetarange_SCI.first),yrange_Z.second);

      TCanvas* can = new TCanvas("can","",2400,1200);
      can->Divide(4,2);
      can->cd(1);
      gPad->SetLogy(1);
      hists_acts.h_N_pos->Draw("hist,same");
      hists_genfit.h_N_pos->Draw("hist,same");
      hists_truth.h_N_pos->Draw("hist,same");
      can->cd(2);
      hists_acts.h_P_pos->SetMinimum(-0.01); hists_acts.h_P_pos->SetMaximum(0.03);
      hists_acts.h_P_pos->Draw("hist,same");
      hists_genfit.h_P_pos->Draw("hist,same");
      hists_truth.h_P_pos->Draw("hist,same");
      if (l_nco_i_P) l_nco_i_P->Draw();
      if (l_nco_o_P) l_nco_o_P->Draw();
      l_nci_i_P->Draw();
      l_nci_o_P->Draw();
      can->cd(3);
      gPad->SetLogy(0);
      hists_acts.h_R_pos->SetMinimum(-0.5); hists_acts.h_R_pos->SetMaximum(0.5);
      hists_acts.h_R_pos->Draw("hist,same");
      hists_genfit.h_R_pos->Draw("hist,same");
      hists_truth.h_R_pos->Draw("hist,same");
      if (l_nco_i_R) l_nco_i_R->Draw();
      if (l_nco_o_R) l_nco_o_R->Draw();
      l_nci_i_R->Draw();
      l_nci_o_R->Draw();
      can->cd(4);
      gPad->SetLogy(0);
      hists_acts.h_Z_pos->SetMinimum(-1); hists_acts.h_Z_pos->SetMaximum(1);
      hists_acts.h_Z_pos->Draw("hist,same");
      hists_genfit.h_Z_pos->Draw("hist,same");
      hists_truth.h_Z_pos->Draw("hist,same");
      if (l_nco_i_Z) l_nco_i_Z->Draw();
      if (l_nco_o_Z) l_nco_o_Z->Draw();
      l_nci_i_Z->Draw();
      l_nci_o_Z->Draw();
      can->cd(5);
      gPad->SetLogy(1);
      hists_acts.h_N_neg->Draw("hist,same");
      hists_genfit.h_N_neg->Draw("hist,same");
      hists_truth.h_N_neg->Draw("hist,same");
      can->cd(6);
      hists_acts.h_P_neg->SetMinimum(-0.01); hists_acts.h_P_neg->SetMaximum(0.03);
      hists_acts.h_P_neg->Draw("hist,same");
      hists_genfit.h_P_neg->Draw("hist,same");
      hists_truth.h_P_neg->Draw("hist,same");
      if (l_sco_i_P) l_sco_i_P->Draw();
      if (l_sco_o_P) l_sco_o_P->Draw();
      l_sci_i_P->Draw();
      l_sci_o_P->Draw();
      can->cd(7);
      gPad->SetLogy(0);
      hists_acts.h_R_neg->SetMinimum(-0.5); hists_acts.h_R_neg->SetMaximum(0.5);
      hists_acts.h_R_neg->Draw("hist,same");
      hists_genfit.h_R_neg->Draw("hist,same");
      hists_truth.h_R_neg->Draw("hist,same");
      if (l_sco_i_R) l_sco_i_R->Draw();
      if (l_sco_o_R) l_sco_o_R->Draw();
      l_sci_i_R->Draw();
      l_sci_o_R->Draw();
      can->cd(8);
      gPad->SetLogy(0);
      hists_acts.h_Z_neg->SetMinimum(-1); hists_acts.h_Z_neg->SetMaximum(1);
      hists_acts.h_Z_neg->Draw("hist,same");
      hists_genfit.h_Z_neg->Draw("hist,same");
      hists_truth.h_Z_neg->Draw("hist,same");
      if (l_sco_i_Z) l_sco_i_Z->Draw();
      if (l_sco_o_Z) l_sco_o_Z->Draw();
      l_sci_i_Z->Draw();
      l_sci_o_Z->Draw();

      gPad->RedrawAxis();

      can->Update();
      can->SaveAs(Form("%s/resid_vsZ_from3D_atR%d_%s.pdf", figureDir.Data(), (int)selectR, phiType.Data()));

      TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
      legend->SetHeader("CPM + TRUTH Seeding, No distortion input");
      legend->AddEntry(hists_acts.h_P_pos, Form("Acts"), "F");
      legend->AddEntry(hists_genfit.h_P_pos, Form("Genfit"), "F");
      legend->AddEntry(hists_truth.h_P_pos, Form("Truth fitter"), "F");
      legend->Draw();
      TCanvas* can_leg = new TCanvas("can_leg","",2000,1200);
      legend->Draw();
      can_leg->SaveAs(Form("%s/resid_vsZ_from3D_leg_%s.pdf", figureDir.Data(), phiType.Data()));

      delete can;
      delete can_leg;

    }
  }

}
