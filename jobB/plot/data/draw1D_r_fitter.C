#include "utilities.h"

void draw1D_r_fitter(bool useWeighted = false)
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
  oneDHist hists_acts{};
  oneDHist hists_genfit{};
  oneDHist hists_CM_Devon_75073{};
  oneDHist hists_CM_Devon_83319{};
  oneDHist hists_PHGARFIELD_Hugo{};
  oneDHist hists_PHGARFIELD_Hugo_new{};

  std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
  for (const auto& selectZ : selectZs)
  {
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "CPM_data_reco_FITacts", hists_acts, kBlack);
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "CPM_data_reco_FITgenfit", hists_genfit, kRed - 7);
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "CM_Devon_reference75073", hists_CM_Devon_75073, kGreen + 3, 1, false);
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "CM_Devon_reference83319", hists_CM_Devon_83319, kGreen + 3, 2, false);
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "PHGARFIELD_Hugo", hists_PHGARFIELD_Hugo, kAzure - 2, 1, false);
    get1Dhist(Form("%s/hist_1Dr_Z%d_%s.root", inputDir.Data(), (int)selectZ, phiType.Data()), "PHGARFIELD_Hugo_new", hists_PHGARFIELD_Hugo_new, kAzure - 2, 2, false);

    if (!hists_acts.h_R_pos || !hists_genfit.h_R_pos || !hists_CM_Devon_75073.h_R_pos || !hists_CM_Devon_83319.h_R_pos || !hists_PHGARFIELD_Hugo.h_R_pos || !hists_PHGARFIELD_Hugo_new.h_R_pos) {
      std::cerr << "Error: some histograms are null for Z = " << selectZ << std::endl;
      continue;
    }

    std::vector<oneDHist> allHists = {
      hists_acts,
      hists_genfit,
      hists_CM_Devon_75073,
      hists_CM_Devon_83319,
      hists_PHGARFIELD_Hugo,
      hists_PHGARFIELD_Hugo_new
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
    hists_acts.h_N_pos->Draw("hist,same");
    hists_genfit.h_N_pos->Draw("hist,same");
    //hists_CM_Devon_75073.h_N_pos->Draw("hist,same");
    //hists_CM_Devon_83319.h_N_pos->Draw("hist,same");
    //hists_PHGARFIELD_Hugo.h_N_pos->Draw("hist,same");
    //hists_PHGARFIELD_Hugo_new.h_N_pos->Draw("hist,same");
    can->cd(2);
    hists_acts.h_P_pos->SetMaximum(15e-3);
    hists_acts.h_P_pos->Draw("hist,same");
    hists_genfit.h_P_pos->Draw("hist,same");
    hists_CM_Devon_75073.h_P_pos->Draw("hist,same");
    hists_CM_Devon_83319.h_P_pos->Draw("hist,same");
    hists_PHGARFIELD_Hugo.h_P_pos->Draw("hist,same");
    hists_PHGARFIELD_Hugo_new.h_P_pos->Draw("hist,same");
    can->cd(3);
    gPad->SetLogy(0);
    hists_acts.h_R_pos->Draw("hist,same");
    hists_genfit.h_R_pos->Draw("hist,same");
    hists_CM_Devon_75073.h_R_pos->Draw("hist,same");
    hists_CM_Devon_83319.h_R_pos->Draw("hist,same");
    hists_PHGARFIELD_Hugo.h_R_pos->Draw("hist,same");
    hists_PHGARFIELD_Hugo_new.h_R_pos->Draw("hist,same");
    can->cd(4);
    gPad->SetLogy(0);
    hists_acts.h_Z_pos->Draw("hist,same");
    hists_genfit.h_Z_pos->Draw("hist,same");
    hists_CM_Devon_75073.h_Z_pos->Draw("hist,same");
    hists_CM_Devon_83319.h_Z_pos->Draw("hist,same");
    hists_PHGARFIELD_Hugo.h_Z_pos->Draw("hist,same");
    hists_PHGARFIELD_Hugo_new.h_Z_pos->Draw("hist,same");
    can->cd(5);
    gPad->SetLogy(1);
    hists_acts.h_N_neg->Draw("hist,same");
    hists_genfit.h_N_neg->Draw("hist,same");
    //hists_CM_Devon_75073.h_N_neg->Draw("hist,same");
    //hists_CM_Devon_83319.h_N_neg->Draw("hist,same");
    //hists_PHGARFIELD_Hugo.h_N_neg->Draw("hist,same");
    //hists_PHGARFIELD_Hugo_new.h_N_neg->Draw("hist,same");
    can->cd(6);
    hists_acts.h_P_neg->SetMaximum(15e-3);
    hists_acts.h_P_neg->Draw("hist,same");
    hists_genfit.h_P_neg->Draw("hist,same");
    hists_CM_Devon_75073.h_P_neg->Draw("hist,same");
    hists_CM_Devon_83319.h_P_neg->Draw("hist,same");
    hists_PHGARFIELD_Hugo.h_P_neg->Draw("hist,same");
    hists_PHGARFIELD_Hugo_new.h_P_neg->Draw("hist,same");
    can->cd(7);
    gPad->SetLogy(0);
    hists_acts.h_R_neg->Draw("hist,same");
    hists_genfit.h_R_neg->Draw("hist,same");
    hists_CM_Devon_75073.h_R_neg->Draw("hist,same");
    hists_CM_Devon_83319.h_R_neg->Draw("hist,same");
    hists_PHGARFIELD_Hugo.h_R_neg->Draw("hist,same");
    hists_PHGARFIELD_Hugo_new.h_R_neg->Draw("hist,same");
    can->cd(8);
    gPad->SetLogy(0);
    hists_acts.h_Z_neg->Draw("hist,same");
    hists_genfit.h_Z_neg->Draw("hist,same");
    hists_CM_Devon_75073.h_Z_neg->Draw("hist,same");
    hists_CM_Devon_83319.h_Z_neg->Draw("hist,same");
    hists_PHGARFIELD_Hugo.h_Z_neg->Draw("hist,same");
    hists_PHGARFIELD_Hugo_new.h_Z_neg->Draw("hist,same");

    gPad->RedrawAxis();

    can->Update();
    can->SaveAs(Form("%s/resid_vsR_from3D_atZ%d_%s.pdf", figureDir.Data(), (int)selectZ,phiType.Data()));

    TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend->SetHeader("Run 79516");
    legend->AddEntry(hists_acts.h_P_pos, Form("CPM + Acts"), "F");
    legend->AddEntry(hists_genfit.h_P_pos, Form("CPM + Genfit"), "F");
    legend->AddEntry(hists_CM_Devon_75073.h_P_pos, Form("Devon CM reference 75073"), "F");
    legend->AddEntry(hists_CM_Devon_83319.h_P_pos, Form("Devon CM reference 83319"), "F");
    legend->AddEntry(hists_PHGARFIELD_Hugo.h_P_pos, Form("Hugo PHGARFIELD"), "F");
    legend->AddEntry(hists_PHGARFIELD_Hugo_new.h_P_pos, Form("Hugo PHGARFIELD new"), "F");
    legend->Draw();
    TCanvas* can_leg = new TCanvas("can_leg","",2000,1200);
    legend->Draw();
    can_leg->SaveAs(Form("%s/resid_vsR_from3D_leg_%s.pdf", figureDir.Data(), phiType.Data()));

    delete can;
    delete can_leg;

  }
  }

}
