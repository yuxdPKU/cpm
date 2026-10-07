#include "utilities.h"

void draw2D_pr_fitter(bool useWeighted = false)
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
  const TString weightType = useWeighted ? "weighted" : "unweighted";
  const TString inputDir = Form("Rootfiles/%s", weightType.Data());
  const TString figureDir = Form("figure/%s", weightType.Data());
  gSystem->mkdir(figureDir, true);

  std::vector<float> selectZs={5, 10, 15, 30, 60};
  for (const auto& selectZ : selectZs)
  {

    twoDHist hists{};
    get2Dhist(Form("%s/hist_2Dpr_Z%d.root", inputDir.Data(), (int)selectZ), "CPM_data_reco_FITacts_USEMMSfalse", hists);

    if (!hists.h_R_pos) {
      std::cerr << "Error: some histograms are null for R = " << selectZ << std::endl;
      continue;
    }

    std::vector<twoDHist> allHists = {
      hists
    };

    std::vector<TH2*> hists_N; hists_N.clear();
    std::vector<TH2*> hists_P; hists_P.clear();
    std::vector<TH2*> hists_R; hists_R.clear();
    std::vector<TH2*> hists_Z; hists_Z.clear();

    AddTwoDHistVectorToVectors(allHists, hists_N, hists_P, hists_R, hists_Z);

    //std::pair<double,double> zrange_N = SetCommonZRange(hists_N);
    std::pair<double,double> zrange_P = SetCommonZRange(hists_P);
    std::pair<double,double> zrange_R = SetCommonZRange(hists_R);
    std::pair<double,double> zrange_Z = SetCommonZRange(hists_Z);

    TCanvas* can = new TCanvas("can","",6000,2400);
    can->Divide(4,2);
    can->cd(1);
    gPad->SetRightMargin(0.15);
    hists.h_N_pos->Draw("colz");
    can->cd(2);
    gPad->SetRightMargin(0.15);
    hists.h_P_pos->Draw("colz");
    can->cd(3);
    gPad->SetRightMargin(0.15);
    hists.h_R_pos->Draw("colz");
    can->cd(4);
    gPad->SetRightMargin(0.15);
    hists.h_Z_pos->Draw("colz");
    can->cd(5);
    gPad->SetRightMargin(0.15);
    hists.h_N_neg->Draw("colz");
    can->cd(6);
    gPad->SetRightMargin(0.15);
    hists.h_P_neg->Draw("colz");
    can->cd(7);
    gPad->SetRightMargin(0.15);
    hists.h_R_neg->Draw("colz");
    can->cd(8);
    gPad->SetRightMargin(0.15);
    hists.h_Z_neg->Draw("colz");

    gPad->RedrawAxis();

    can->Update();
    can->SaveAs(Form("%s/resid_vsPR_from3D_Z%d.pdf",figureDir.Data(),(int)selectZ));

    TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend->SetHeader("ana573_2026p003_v001, Run 79516, poly seed");
    legend->AddEntry(hists.h_P_pos, Form("No TPOT"), "F");
    legend->Draw();
    TCanvas* can_leg = new TCanvas("can_leg","",2000,1200);
    legend->Draw();
    can_leg->SaveAs(Form("%s/resid_vsPR_from3D_leg.pdf",figureDir.Data()));

    delete can;
    delete can_leg;

  }

}
