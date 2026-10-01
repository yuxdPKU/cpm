#include "utilities.h"

void draw1D_r_method_acts(bool useWeighted = false)
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
  const TString weightType = useWeighted ? "weighted" : "unweighted";
  const TString inputDir = Form("Rootfiles/%s", weightType.Data());
  const TString figureDir = Form("figure/%s", weightType.Data());
  gSystem->mkdir(figureDir, true);

  int run = 79516;
  std::vector<TString> phiTypes = {"central", "west", "east"};
  for (const auto& phiType : phiTypes)
  {
  oneDHist hists_mi{};
  oneDHist hists_cpm{};
  oneDHist hists_lamination{};

  std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
  for (const auto& selectZ : selectZs)
  {
    get1Dhist(Form("/gpfs/mnt/gpfs02/sphenix/user/xyu3/TPCdistortion/Si_TPOT_fit/run3pp_extrapolate/jobB/Rootfiles/hist_1Dr_Z%d_%s.root", (int)selectZ, phiType.Data()), "CPM_data_reco_FITacts_EXTRAbidirectional", hists_mi, kRed - 7);
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "CPM_data_reco_FITacts_EXTRAbidirectional", hists_cpm, kGreen + 3);
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "Lamination", hists_lamination, kAzure - 2, 1, false);

    if (!hists_mi.h_N_pos || !hists_cpm.h_N_pos || !hists_lamination.h_P_pos) {
      std::cerr << "Error: some histograms are null for Z = " << selectZ << std::endl;
      continue;
    }

    std::vector<oneDHist> allHists = {
      hists_mi,
      hists_cpm,
      hists_lamination
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
    hists_mi.h_N_pos->Draw("hist,same");
    hists_cpm.h_N_pos->Draw("hist,same");
    can->cd(2);
    hists_mi.h_P_pos->Draw("hist,same");
    hists_cpm.h_P_pos->Draw("hist,same");
    hists_lamination.h_P_pos->Draw("hist,same");
    can->cd(3);
    gPad->SetLogy(0);
    hists_mi.h_R_pos->Draw("hist,same");
    hists_cpm.h_R_pos->Draw("hist,same");
    hists_lamination.h_R_pos->Draw("hist,same");
    can->cd(4);
    gPad->SetLogy(0);
    hists_mi.h_Z_pos->Draw("hist,same");
    hists_cpm.h_Z_pos->Draw("hist,same");
    hists_lamination.h_Z_pos->Draw("hist,same");
    can->cd(5);
    gPad->SetLogy(1);
    hists_mi.h_N_neg->Draw("hist,same");
    hists_cpm.h_N_neg->Draw("hist,same");
    can->cd(6);
    hists_mi.h_P_neg->Draw("hist,same");
    hists_cpm.h_P_neg->Draw("hist,same");
    hists_lamination.h_P_neg->Draw("hist,same");
    can->cd(7);
    gPad->SetLogy(0);
    hists_mi.h_R_neg->Draw("hist,same");
    hists_cpm.h_R_neg->Draw("hist,same");
    hists_lamination.h_R_neg->Draw("hist,same");
    can->cd(8);
    gPad->SetLogy(0);
    hists_mi.h_Z_neg->Draw("hist,same");
    hists_cpm.h_Z_neg->Draw("hist,same");
    hists_lamination.h_Z_neg->Draw("hist,same");

    gPad->RedrawAxis();

    can->Update();
    can->SaveAs(Form("%s/resid_vsR_from3D_atZ%d_%s.pdf", figureDir.Data(), (int)selectZ, phiType.Data()));

    TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend->SetHeader("Run 79516, Acts (bidirectional, cluster error data)");
    legend->AddEntry(hists_mi.h_P_pos, Form("MI"), "F");
    legend->AddEntry(hists_cpm.h_P_pos, Form("CPM"), "F");
    legend->AddEntry(hists_lamination.h_P_pos, Form("Lamination"), "F");
    legend->Draw();
    TCanvas* can_leg = new TCanvas("can_leg","",2000,1200);
    legend->Draw();
    can_leg->SaveAs(Form("%s/resid_vsR_from3D_leg_%s.pdf", figureDir.Data(), phiType.Data()));

    delete can;
    delete can_leg;

  }
  }

}
