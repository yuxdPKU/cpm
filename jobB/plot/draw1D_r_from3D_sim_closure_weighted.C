#include <TGaxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>

#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

struct threeDHist {
    TH3* h_N_prz_pos;
    TH3* h_R_prz_pos;
    TH3* h_P_prz_pos;
    TH3* h_Z_prz_pos;
    TH3* h_N_prz_neg;
    TH3* h_R_prz_neg;
    TH3* h_P_prz_neg;
    TH3* h_Z_prz_neg;
};

struct oneDHist {
    TH1* h_N_pos;
    TH1* h_R_pos;
    TH1* h_P_pos;
    TH1* h_Z_pos;
    TH1* h_N_neg;
    TH1* h_R_neg;
    TH1* h_P_neg;
    TH1* h_Z_neg;
};

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);
void draw1Dmap(TString filename, TString tag, double selectZ, float selectPhi, oneDHist& hists1D, int color=1, bool convert_RP_2_P=false);
void draw1Dmap_2D(TString filename, TString tag, TH1*& h_R_pos, TH1*& h_R_neg, TH1*& h_P_pos, TH1*& h_P_neg, TH1*& h_Z_pos, TH1*& h_Z_neg, int color);

void plot1D_Zbin(TH3* h3, TH1* h1, float z, float selectPhi, bool convert_RP_2_P=false)
{
  int xbin = h3->GetXaxis()->FindBin(selectPhi);
  int zbin = h3->GetZaxis()->FindBin(z);
  for (int i = 1; i <= h3->GetNbinsY(); i++)
  {
      double value = h3->GetBinContent(xbin, i, zbin);
      double error = h3->GetBinError(xbin, i, zbin);
      if (convert_RP_2_P)
      {
        double rcenter = h3->GetYaxis()->GetBinCenter(i);
        value /= rcenter;
        error /= rcenter;
      }
      h1->SetBinContent(i, value);
      h1->SetBinError(i, error);
  }
}

void plot1D(TH2* h2, TH1* h1)
{
  int xbin = h2->GetXaxis()->FindBin(((-1.72932)+(-1.42713))/2.+2*TMath::Pi());
  for (int i = 1; i <= h2->GetNbinsY(); i++)
  {
      int xminbin = h2->GetXaxis()->FindBin((-1.72932)+2*TMath::Pi());
      int xmaxbin = h2->GetXaxis()->FindBin((-1.42713)+2*TMath::Pi());
      float bincontent=0, binerror=0;
      int nbin=0;
      for (int j = xminbin; j <= xmaxbin; j++)
      {
        bincontent+=h2->GetBinContent(j, i);
        binerror+=h2->GetBinError(j, i);
	      nbin++;
      }
      bincontent/=nbin;
      binerror/=nbin;
      h1->SetBinContent(i, bincontent);
      h1->SetBinError(i, binerror);
  }
}

void draw1D_r_from3D_sim_closure_weighted(TString phiType="central")
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

int run = 29;
phiType.ToLower();
float selectPhi = 0;
if (phiType.CompareTo("central")==0)
{
  selectPhi = ((-1.72932)+(-1.42713))/2.+2*TMath::Pi();
}
else if (phiType.CompareTo("west")==0)
{
  selectPhi = ((-1.20634)+(-0.902555))/2.+2*TMath::Pi();
}
else if (phiType.CompareTo("east")==0)
{
  selectPhi = ((-2.25883)+(-1.9569))/2.+2*TMath::Pi();
}
else
{
  std::cerr << "Error: unknown phiType = " << phiType << std::endl;
  std::cerr << "Please use phiType = central, west, or east." << std::endl;
  return;
}

std::cout << "Use phiType = " << phiType
          << ", phi = " << selectPhi
          << std::endl;

std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
for (const auto& selectZ : selectZs)
{

oneDHist hists_CPM_sim_acts, hists_CPM_sim_genfit, hists_CPM_sim_truth;
draw1Dmap(Form("../output/sim_run%d_acts_weighted/sim_run%d_acts_weighted_maxdca0p2_B3_average_correction_histograms.root", run, run), Form("CPM_sim_reco_acts"), selectZ, selectPhi, hists_CPM_sim_acts, 2);
draw1Dmap(Form("../output/sim_run%d_genfit_weighted/sim_run%d_genfit_weighted_maxdca0p2_B3_average_correction_histograms.root", run, run), Form("CPM_sim_reco_genfit"), selectZ, selectPhi, hists_CPM_sim_genfit, 3);
draw1Dmap(Form("../output/sim_run%d_truth_weighted/sim_run%d_truth_weighted_maxdca0p2_B3_average_correction_histograms.root", run, run), Form("CPM_sim_reco_truth"), selectZ, selectPhi, hists_CPM_sim_truth, 4);

std::vector<TH1*> hists_P; hists_P.clear();
std::vector<TH1*> hists_R; hists_R.clear();
std::vector<TH1*> hists_Z; hists_Z.clear();
hists_P.push_back(hists_CPM_sim_acts.h_P_neg);
hists_P.push_back(hists_CPM_sim_acts.h_P_pos);
hists_P.push_back(hists_CPM_sim_genfit.h_P_neg);
hists_P.push_back(hists_CPM_sim_genfit.h_P_pos);
hists_P.push_back(hists_CPM_sim_truth.h_P_neg);
hists_P.push_back(hists_CPM_sim_truth.h_P_pos);
hists_R.push_back(hists_CPM_sim_acts.h_R_neg);
hists_R.push_back(hists_CPM_sim_acts.h_R_pos);
hists_R.push_back(hists_CPM_sim_genfit.h_R_neg);
hists_R.push_back(hists_CPM_sim_genfit.h_R_pos);
hists_R.push_back(hists_CPM_sim_truth.h_R_neg);
hists_R.push_back(hists_CPM_sim_truth.h_R_pos);
hists_Z.push_back(hists_CPM_sim_acts.h_Z_neg);
hists_Z.push_back(hists_CPM_sim_acts.h_Z_pos);
hists_Z.push_back(hists_CPM_sim_genfit.h_Z_neg);
hists_Z.push_back(hists_CPM_sim_genfit.h_Z_pos);
hists_Z.push_back(hists_CPM_sim_truth.h_Z_neg);
hists_Z.push_back(hists_CPM_sim_truth.h_Z_pos);
std::pair<double,double> yrange_P = SetCommonYRange(hists_P);
std::pair<double,double> yrange_R = SetCommonYRange(hists_R);
std::pair<double,double> yrange_Z = SetCommonYRange(hists_Z);

TFile* ofile_N = new TFile(Form("Rootfiles/hist_n_1Dr_Z%d_weighted_%s.root",(int)selectZ,phiType.Data()),"recreate");
ofile_N->cd();
hists_CPM_sim_acts.h_N_pos->Write();
hists_CPM_sim_acts.h_N_neg->Write();
hists_CPM_sim_genfit.h_N_pos->Write();
hists_CPM_sim_genfit.h_N_neg->Write();
hists_CPM_sim_truth.h_N_pos->Write();
hists_CPM_sim_truth.h_N_neg->Write();
ofile_N->Write();

TFile* ofile_R = new TFile(Form("Rootfiles/hist_dr_1Dr_Z%d_weighted_%s.root",(int)selectZ,phiType.Data()),"recreate");
ofile_R->cd();
hists_CPM_sim_acts.h_R_pos->Write();
hists_CPM_sim_acts.h_R_neg->Write();
hists_CPM_sim_genfit.h_R_pos->Write();
hists_CPM_sim_genfit.h_R_neg->Write();
hists_CPM_sim_truth.h_R_pos->Write();
hists_CPM_sim_truth.h_R_neg->Write();
ofile_R->Write();

TFile* ofile_P = new TFile(Form("Rootfiles/hist_dphi_1Dr_Z%d_weighted_%s.root",(int)selectZ,phiType.Data()),"recreate");
ofile_P->cd();
hists_CPM_sim_acts.h_P_pos->Write();
hists_CPM_sim_acts.h_P_neg->Write();
hists_CPM_sim_genfit.h_P_pos->Write();
hists_CPM_sim_genfit.h_P_neg->Write();
hists_CPM_sim_truth.h_P_pos->Write();
hists_CPM_sim_truth.h_P_neg->Write();
ofile_P->Write();

TFile* ofile_Z = new TFile(Form("Rootfiles/hist_dz_1Dr_Z%d_weighted_%s.root",(int)selectZ,phiType.Data()),"recreate");
ofile_Z->cd();
hists_CPM_sim_acts.h_Z_pos->Write();
hists_CPM_sim_acts.h_Z_neg->Write();
hists_CPM_sim_genfit.h_Z_pos->Write();
hists_CPM_sim_genfit.h_Z_neg->Write();
hists_CPM_sim_truth.h_Z_pos->Write();
hists_CPM_sim_truth.h_Z_neg->Write();
ofile_Z->Write();

TCanvas* can = new TCanvas("can","",2400,1200);
can->Divide(4,2);
can->cd(1);
gPad->SetLogy(1);
hists_CPM_sim_acts.h_N_pos->Draw("hist,same");
hists_CPM_sim_genfit.h_N_pos->Draw("hist,same");
hists_CPM_sim_truth.h_N_pos->Draw("hist,same");
can->cd(2);
hists_CPM_sim_acts.h_P_pos->SetMinimum(-0.01); hists_CPM_sim_acts.h_P_pos->SetMaximum(0.03);
hists_CPM_sim_acts.h_P_pos->Draw("hist,same");
hists_CPM_sim_genfit.h_P_pos->Draw("hist,same");
hists_CPM_sim_truth.h_P_pos->Draw("hist,same");
can->cd(3);
gPad->SetLogy(0);
hists_CPM_sim_acts.h_R_pos->SetMinimum(-0.5); hists_CPM_sim_acts.h_R_pos->SetMaximum(0.5);
hists_CPM_sim_acts.h_R_pos->Draw("hist,same");
hists_CPM_sim_genfit.h_R_pos->Draw("hist,same");
hists_CPM_sim_truth.h_R_pos->Draw("hist,same");
can->cd(4);
gPad->SetLogy(0);
hists_CPM_sim_acts.h_Z_pos->SetMinimum(-0.1); hists_CPM_sim_acts.h_Z_pos->SetMaximum(0.1);
hists_CPM_sim_acts.h_Z_pos->Draw("hist,same");
hists_CPM_sim_genfit.h_Z_pos->Draw("hist,same");
hists_CPM_sim_truth.h_Z_pos->Draw("hist,same");
can->cd(5);
gPad->SetLogy(1);
hists_CPM_sim_acts.h_N_neg->Draw("hist,same");
hists_CPM_sim_genfit.h_N_neg->Draw("hist,same");
hists_CPM_sim_truth.h_N_neg->Draw("hist,same");
can->cd(6);
hists_CPM_sim_acts.h_P_neg->SetMinimum(-0.01); hists_CPM_sim_acts.h_P_neg->SetMaximum(0.03);
hists_CPM_sim_acts.h_P_neg->Draw("hist,same");
hists_CPM_sim_genfit.h_P_neg->Draw("hist,same");
hists_CPM_sim_truth.h_P_neg->Draw("hist,same");
can->cd(7);
gPad->SetLogy(0);
hists_CPM_sim_acts.h_R_neg->SetMinimum(-0.5); hists_CPM_sim_acts.h_R_neg->SetMaximum(0.5);
hists_CPM_sim_acts.h_R_neg->Draw("hist,same");
hists_CPM_sim_genfit.h_R_neg->Draw("hist,same");
hists_CPM_sim_truth.h_R_neg->Draw("hist,same");
can->cd(8);
gPad->SetLogy(0);
hists_CPM_sim_acts.h_Z_neg->SetMinimum(-0.1); hists_CPM_sim_acts.h_Z_neg->SetMaximum(0.1);
hists_CPM_sim_acts.h_Z_neg->Draw("hist,same");
hists_CPM_sim_genfit.h_Z_neg->Draw("hist,same");
hists_CPM_sim_truth.h_Z_neg->Draw("hist,same");

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_weighted_%s.pdf",(int)selectZ,phiType.Data()));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->SetHeader("CPM weighted");
legend->AddEntry(hists_CPM_sim_acts.h_P_pos, Form("ACTS"), "F");
legend->AddEntry(hists_CPM_sim_genfit.h_P_pos, Form("GENFIT"), "F");
legend->AddEntry(hists_CPM_sim_truth.h_P_pos, Form("TRUTH"), "F");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs(Form("figure/resid_vsR_from3D_leg_weighted_%s.pdf",phiType.Data()));

delete can;
delete can_leg;

}

}

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms)
{
    Double_t yMin = TMath::Infinity();
    Double_t yMax = -TMath::Infinity();

    for (TH1* h : histograms) {
        if (!h) continue;
        yMin = TMath::Min(yMin, h->GetMinimum());
        yMax = TMath::Max(yMax, h->GetMaximum());
    }

    if (yMin == TMath::Infinity()) yMin = 0;
    if (yMax == -TMath::Infinity()) yMax = 1;

    for (TH1* h : histograms) {
        if (h) {
            h->SetMinimum(yMin);
            h->SetMaximum(yMax);
        }
    }
    return {yMin, yMax};
}

TH1* SubtractHistograms(const TH1* h1, const TH1* h2, const char* name = "diff")
{
    if (!h1 || !h2) return nullptr;
    if (h1->GetNbinsX() != h2->GetNbinsX()) {
        std::cerr << "Error: histograms have different number of bins." << std::endl;
        return nullptr;
    }
    TH1* hdiff = (TH1*)h1->Clone(name);
    hdiff->Reset();
    for (int i = 1; i <= h1->GetNbinsX(); ++i) {
        double val1 = h1->GetBinContent(i);
        double err1 = h1->GetBinError(i);
        double val2 = h2->GetBinContent(i);
        double err2 = h2->GetBinError(i);
        double diff = val1 - val2;
        double err = sqrt(err1*err1 + err2*err2);
        hdiff->SetBinContent(i, diff);
        hdiff->SetBinError(i, err);
    }
    return hdiff;
}

void draw1Dmap(TString filename, TString tag, double selectZ, float selectPhi,
  oneDHist& hists1D, int color, bool convert_RP_2_P)
{
  TFile* file_3D_map = new TFile(filename,"");
  if (!file_3D_map || file_3D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  threeDHist hists3D;
  hists3D.h_N_prz_pos = (TH3*) file_3D_map->Get("hentries_posz");
  hists3D.h_R_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  hists3D.h_P_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  hists3D.h_Z_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  hists3D.h_N_prz_neg = (TH3*) file_3D_map->Get("hentries_negz");
  hists3D.h_R_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  hists3D.h_P_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  hists3D.h_Z_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");

  if (!hists3D.h_N_prz_pos || !hists3D.h_R_prz_pos || !hists3D.h_P_prz_pos || !hists3D.h_Z_prz_pos ||
      !hists3D.h_N_prz_neg || !hists3D.h_R_prz_neg || !hists3D.h_P_prz_neg || !hists3D.h_Z_prz_neg) {
      std::cerr << "Error: missing 3D histograms in file" << std::endl;
      delete file_3D_map;
      return;
  }

  hists1D.h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),hists3D.h_N_prz_pos->GetYaxis()->GetNbins(),hists3D.h_N_prz_pos->GetYaxis()->GetXmin(),hists3D.h_N_prz_pos->GetYaxis()->GetXmax());
  hists1D.h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),hists3D.h_R_prz_pos->GetYaxis()->GetNbins(),hists3D.h_R_prz_pos->GetYaxis()->GetXmin(),hists3D.h_R_prz_pos->GetYaxis()->GetXmax());
  hists1D.h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),hists3D.h_P_prz_pos->GetYaxis()->GetNbins(),hists3D.h_P_prz_pos->GetYaxis()->GetXmin(),hists3D.h_P_prz_pos->GetYaxis()->GetXmax());
  hists1D.h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),hists3D.h_Z_prz_pos->GetYaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetYaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetYaxis()->GetXmax());
  hists1D.h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),hists3D.h_N_prz_neg->GetYaxis()->GetNbins(),hists3D.h_N_prz_neg->GetYaxis()->GetXmin(),hists3D.h_N_prz_neg->GetYaxis()->GetXmax());
  hists1D.h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),hists3D.h_R_prz_neg->GetYaxis()->GetNbins(),hists3D.h_R_prz_neg->GetYaxis()->GetXmin(),hists3D.h_R_prz_neg->GetYaxis()->GetXmax());
  hists1D.h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ Z=-%d cm;R (cm);dphi (rad)",(int)selectZ),hists3D.h_P_prz_neg->GetYaxis()->GetNbins(),hists3D.h_P_prz_neg->GetYaxis()->GetXmin(),hists3D.h_P_prz_neg->GetYaxis()->GetXmax());
  hists1D.h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),hists3D.h_Z_prz_neg->GetYaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetYaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetYaxis()->GetXmax());

  // do not save in the file_3D_map
  // directly saved in the stack memory
  hists1D.h_N_pos->SetDirectory(0);
  hists1D.h_R_pos->SetDirectory(0);
  hists1D.h_P_pos->SetDirectory(0);
  hists1D.h_Z_pos->SetDirectory(0);
  hists1D.h_N_neg->SetDirectory(0);
  hists1D.h_R_neg->SetDirectory(0);
  hists1D.h_P_neg->SetDirectory(0);
  hists1D.h_Z_neg->SetDirectory(0);

  plot1D_Zbin(hists3D.h_N_prz_pos,hists1D.h_N_pos,selectZ, selectPhi, false);
  plot1D_Zbin(hists3D.h_R_prz_pos,hists1D.h_R_pos,selectZ, selectPhi, false);
  plot1D_Zbin(hists3D.h_P_prz_pos,hists1D.h_P_pos,selectZ, selectPhi, convert_RP_2_P);
  plot1D_Zbin(hists3D.h_Z_prz_pos,hists1D.h_Z_pos,selectZ, selectPhi, false);
  plot1D_Zbin(hists3D.h_N_prz_neg,hists1D.h_N_neg,-selectZ, selectPhi, false);
  plot1D_Zbin(hists3D.h_R_prz_neg,hists1D.h_R_neg,-selectZ, selectPhi, false);
  plot1D_Zbin(hists3D.h_P_prz_neg,hists1D.h_P_neg,-selectZ, selectPhi, convert_RP_2_P);
  plot1D_Zbin(hists3D.h_Z_prz_neg,hists1D.h_Z_neg,-selectZ, selectPhi, false);

  hists1D.h_N_pos->SetLineColor(color); hists1D.h_N_pos->SetLineWidth(1); hists1D.h_N_pos->SetFillColor(0); hists1D.h_N_pos->SetMarkerColor(color);
  hists1D.h_R_pos->SetLineColor(color); hists1D.h_R_pos->SetLineWidth(1); hists1D.h_R_pos->SetFillColor(0); hists1D.h_R_pos->SetMarkerColor(color);
  hists1D.h_P_pos->SetLineColor(color); hists1D.h_P_pos->SetLineWidth(1); hists1D.h_P_pos->SetFillColor(0); hists1D.h_P_pos->SetMarkerColor(color);
  hists1D.h_Z_pos->SetLineColor(color); hists1D.h_Z_pos->SetLineWidth(1); hists1D.h_Z_pos->SetFillColor(0); hists1D.h_Z_pos->SetMarkerColor(color);
  hists1D.h_N_neg->SetLineColor(color); hists1D.h_N_neg->SetLineWidth(1); hists1D.h_N_neg->SetFillColor(0); hists1D.h_N_neg->SetMarkerColor(color);
  hists1D.h_R_neg->SetLineColor(color); hists1D.h_R_neg->SetLineWidth(1); hists1D.h_R_neg->SetFillColor(0); hists1D.h_R_neg->SetMarkerColor(color);
  hists1D.h_P_neg->SetLineColor(color); hists1D.h_P_neg->SetLineWidth(1); hists1D.h_P_neg->SetFillColor(0); hists1D.h_P_neg->SetMarkerColor(color);
  hists1D.h_Z_neg->SetLineColor(color); hists1D.h_Z_neg->SetLineWidth(1); hists1D.h_Z_neg->SetFillColor(0); hists1D.h_Z_neg->SetMarkerColor(color);

  delete file_3D_map;

  return;
}

void draw1Dmap_2D(TString filename, TString tag,
  TH1*& h_R_pos, TH1*& h_R_neg,
  TH1*& h_P_pos, TH1*& h_P_neg,
  TH1*& h_Z_pos, TH1*& h_Z_neg,
  int color=1)
{
  TFile* file_2D_map = new TFile(filename,"");
  if (!file_2D_map || file_2D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  TH2* h_R_pr_pos = (TH2*) file_2D_map->Get("hIntDistortionR_posz");
  TH2* h_P_pr_pos = (TH2*) file_2D_map->Get("hIntDistortionP_posz");
  TH2* h_Z_pr_pos = (TH2*) file_2D_map->Get("hIntDistortionZ_posz");
  //TH2* h_N_pr_pos = (TH2*) file_2D_map->Get("hentries_posz");
  TH2* h_R_pr_neg = (TH2*) file_2D_map->Get("hIntDistortionR_negz");
  TH2* h_P_pr_neg = (TH2*) file_2D_map->Get("hIntDistortionP_negz");
  TH2* h_Z_pr_neg = (TH2*) file_2D_map->Get("hIntDistortionZ_negz");
  //TH2* h_N_pr_neg = (TH2*) file_2D_map->Get("hentries_negz");

  if (!h_R_pr_pos || !h_P_pr_pos || !h_Z_pr_pos ||
      !h_R_pr_neg || !h_P_pr_neg || !h_Z_pr_neg) {
      std::cerr << "Error: missing 2D histograms in file" << std::endl;
      delete file_2D_map;
      return;
  }

  //h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N;R (cm);N"),h_N_pr_pos->GetYaxis()->GetNbins(),h_N_pr_pos->GetYaxis()->GetXmin(),h_N_pr_pos->GetYaxis()->GetXmax());
  h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR;R (cm);dR (cm)"),h_R_pr_pos->GetYaxis()->GetNbins(),h_R_pr_pos->GetYaxis()->GetXmin(),h_R_pr_pos->GetYaxis()->GetXmax());
  h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi;R (cm);dphi (rad)"),h_P_pr_pos->GetYaxis()->GetNbins(),h_P_pr_pos->GetYaxis()->GetXmin(),h_P_pr_pos->GetYaxis()->GetXmax());
  h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz;R (cm);dz (cm)"),h_Z_pr_pos->GetYaxis()->GetNbins(),h_Z_pr_pos->GetYaxis()->GetXmin(),h_Z_pr_pos->GetYaxis()->GetXmax());
  //h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N;R (cm);N"),h_N_pr_neg->GetYaxis()->GetNbins(),h_N_pr_neg->GetYaxis()->GetXmin(),h_N_pr_neg->GetYaxis()->GetXmax());
  h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR;R (cm);dR (cm)"),h_R_pr_neg->GetYaxis()->GetNbins(),h_R_pr_neg->GetYaxis()->GetXmin(),h_R_pr_neg->GetYaxis()->GetXmax());
  h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi;R (cm);dphi (rad)"),h_P_pr_neg->GetYaxis()->GetNbins(),h_P_pr_neg->GetYaxis()->GetXmin(),h_P_pr_neg->GetYaxis()->GetXmax());
  h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz;R (cm);dz (cm)"),h_Z_pr_neg->GetYaxis()->GetNbins(),h_Z_pr_neg->GetYaxis()->GetXmin(),h_Z_pr_neg->GetYaxis()->GetXmax());

  // do not save in the file_2D_map
  // directly saved in the stack memory
  //h_N_pos->SetDirectory(0);
  h_R_pos->SetDirectory(0);
  h_P_pos->SetDirectory(0);
  h_Z_pos->SetDirectory(0);
  //h_N_neg->SetDirectory(0);
  h_R_neg->SetDirectory(0);
  h_P_neg->SetDirectory(0);
  h_Z_neg->SetDirectory(0);

  //plot1D(h_N_pr_pos,h_N_pos);
  plot1D(h_R_pr_pos,h_R_pos);
  plot1D(h_P_pr_pos,h_P_pos);
  plot1D(h_Z_pr_pos,h_Z_pos);
  //plot1D(h_N_pr_neg,h_N_neg);
  plot1D(h_R_pr_neg,h_R_neg);
  plot1D(h_P_pr_neg,h_P_neg);
  plot1D(h_Z_pr_neg,h_Z_neg);

  //h_N_pos->SetLineColor(color); h_N_pos->SetLineWidth(1); h_N_pos->SetFillColor(0); h_N_pos->SetMarkerColor(color);
  h_R_pos->SetLineColor(color); h_R_pos->SetLineWidth(1); h_R_pos->SetFillColor(0); h_R_pos->SetMarkerColor(color);
  h_P_pos->SetLineColor(color); h_P_pos->SetLineWidth(1); h_P_pos->SetFillColor(0); h_P_pos->SetMarkerColor(color);
  h_Z_pos->SetLineColor(color); h_Z_pos->SetLineWidth(1); h_Z_pos->SetFillColor(0); h_Z_pos->SetMarkerColor(color);
  //h_N_neg->SetLineColor(color); h_N_neg->SetLineWidth(1); h_N_neg->SetFillColor(0); h_R_neg->SetMarkerColor(color);
  h_R_neg->SetLineColor(color); h_R_neg->SetLineWidth(1); h_R_neg->SetFillColor(0); h_R_neg->SetMarkerColor(color);
  h_P_neg->SetLineColor(color); h_P_neg->SetLineWidth(1); h_P_neg->SetFillColor(0); h_P_neg->SetMarkerColor(color);
  h_Z_neg->SetLineColor(color); h_Z_neg->SetLineWidth(1); h_Z_neg->SetFillColor(0); h_Z_neg->SetMarkerColor(color);

  delete file_2D_map;

  return;
}
