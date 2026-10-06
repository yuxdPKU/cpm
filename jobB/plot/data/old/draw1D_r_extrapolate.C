#include "utilities.h"

void draw1D_r_extrapolate(bool useWeighted = false)
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
  oneDHist hists_default;
  oneDHist hists_forward;
  oneDHist hists_backward;
  oneDHist hists_bidirectional;

  std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
  for (const auto& selectZ : selectZs)
  {
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "CPM_data_reco_FITacts_EXTRAdefault", hists_default, kRed-7);
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "CPM_data_reco_FITacts_EXTRAforward", hists_forward, kGreen+3);
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "CPM_data_reco_FITacts_EXTRAbackward", hists_backward, kAzure-2);
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "CPM_data_reco_FITacts_EXTRAbidirectional", hists_bidirectional, kBlack);

    if (!hists_default.h_N_pos || !hists_forward.h_N_pos || 
        !hists_backward.h_N_pos || !hists_bidirectional.h_N_pos) {
      std::cerr << "Error: some histograms are null for Z = " << selectZ << std::endl;
      continue;
    }

    std::vector<oneDHist> allHists = {
      hists_default,
      hists_forward,
      hists_backward,
      hists_bidirectional
    };

    std::vector<TH1*> hists_P; hists_P.clear();
    std::vector<TH1*> hists_R; hists_R.clear();
    std::vector<TH1*> hists_Z; hists_Z.clear();

    AddOneDHistVectorToVectors(allHists, hists_P, hists_R, hists_Z);

    std::pair<double,double> yrange_P = SetCommonYRange(hists_P);
    std::pair<double,double> yrange_R = SetCommonYRange(hists_R);
    std::pair<double,double> yrange_Z = SetCommonYRange(hists_Z);

    TCanvas* can = new TCanvas("can","",2400,1200);
    can->Divide(4,2);
    can->cd(1);
    gPad->SetLogy(1);
    hists_default.h_N_pos->Draw("histsame");
    hists_forward.h_N_pos->Draw("histsame");
    hists_backward.h_N_pos->Draw("histsame");
    hists_bidirectional.h_N_pos->Draw("histsame");
    can->cd(2);
    hists_default.h_P_pos->Draw("histsame");
    hists_forward.h_P_pos->Draw("histsame");
    hists_backward.h_P_pos->Draw("histsame");
    hists_bidirectional.h_P_pos->Draw("histsame");
    can->cd(3);
    gPad->SetLogy(0);
    hists_default.h_R_pos->Draw("histsame");
    hists_forward.h_R_pos->Draw("histsame");
    hists_backward.h_R_pos->Draw("histsame");
    hists_bidirectional.h_R_pos->Draw("histsame");
    can->cd(4);
    gPad->SetLogy(0);
    hists_default.h_Z_pos->Draw("histsame");
    hists_forward.h_Z_pos->Draw("histsame");
    hists_backward.h_Z_pos->Draw("histsame");
    hists_bidirectional.h_Z_pos->Draw("histsame");
    can->cd(5);
    gPad->SetLogy(1);
    hists_default.h_N_neg->Draw("histsame");
    hists_forward.h_N_neg->Draw("histsame");
    hists_backward.h_N_neg->Draw("histsame");
    hists_bidirectional.h_N_neg->Draw("histsame");
    can->cd(6);
    hists_default.h_P_neg->Draw("histsame");
    hists_forward.h_P_neg->Draw("histsame");
    hists_backward.h_P_neg->Draw("histsame");
    hists_bidirectional.h_P_neg->Draw("histsame");
    can->cd(7);
    gPad->SetLogy(0);
    hists_default.h_R_neg->Draw("histsame");
    hists_forward.h_R_neg->Draw("histsame");
    hists_backward.h_R_neg->Draw("histsame");
    hists_bidirectional.h_R_neg->Draw("histsame");
    can->cd(8);
    gPad->SetLogy(0);
    hists_default.h_Z_neg->Draw("histsame");
    hists_forward.h_Z_neg->Draw("histsame");
    hists_backward.h_Z_neg->Draw("histsame");
    hists_bidirectional.h_Z_neg->Draw("histsame");

    gPad->RedrawAxis();

    can->Update();
    can->SaveAs(Form("%s/resid_vsR_from3D_atZ%d_%s.pdf", figureDir.Data(), (int)selectZ, phiType.Data()));

    TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend->SetHeader("Run 79516, CPM + Acts + Cluster error data");
    legend->AddEntry(hists_default.h_P_pos, Form("Default extrapolation"), "F");
    legend->AddEntry(hists_forward.h_P_pos, Form("Forward extrapolation"), "F");
    legend->AddEntry(hists_backward.h_P_pos, Form("Backward extrapolation"), "F");
    legend->AddEntry(hists_bidirectional.h_P_pos, Form("Bidirectional extrapolation"), "F");
    legend->Draw();
    TCanvas* can_leg = new TCanvas("can_leg","",2000,1200);
    legend->Draw();
    can_leg->SaveAs(Form("%s/resid_vsR_from3D_leg_%s.pdf", figureDir.Data(), phiType.Data()));

    delete can;
    delete can_leg;

  }
  }

}
