#include <TGaxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TLine.h>
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
void draw1Dmap(TString filename, TString tag, double selectR, float selectPhi, oneDHist& hists1D, int color=1, bool convert_RP_2_P=false);
void draw1Dmap_2D(TString filename, TString tag, TH1*& h_R_pos, TH1*& h_R_neg, TH1*& h_P_pos, TH1*& h_P_neg, TH1*& h_Z_pos, TH1*& h_Z_neg, int color);
void draw1Dmap_2D(TString filename, TString tag, double selectR, TLine*& l_R_pos, TLine*& l_R_neg, TLine*& l_P_pos, TLine*& l_P_neg, TLine*& l_Z_pos, TLine*& l_Z_neg, int color=1);

void plot1D_Rbin(TH3* h3, TH1* h1, float r, float selectPhi, bool convert_RP_2_P=false)
{
  int xbin = h3->GetXaxis()->FindBin(selectPhi);
  int rbin = h3->GetYaxis()->FindBin(r);
  double rcenter = h3->GetYaxis()->GetBinCenter(rbin);
  for (int i = 1; i <= h3->GetNbinsZ(); i++)
  {
      double value = h3->GetBinContent(xbin, rbin, i);
      double error = h3->GetBinError(xbin, rbin, i);
      if (convert_RP_2_P)
      {
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

void plot1D_line(TH2* h2, TLine*& line, float y, bool is_north_or_south)
{
  int ybin = h2->GetYaxis()->FindBin(y);
  int xminbin = h2->GetXaxis()->FindBin((-1.72932)+2*TMath::Pi());
  int xmaxbin = h2->GetXaxis()->FindBin((-1.42713)+2*TMath::Pi());
  float bincontent=0, binerror=0;
  int nbin=0;
  for (int j = xminbin; j <= xmaxbin; j++)
  {
    bincontent+=h2->GetBinContent(j, ybin);
    binerror+=h2->GetBinError(j, ybin);
    nbin++;
  }
  bincontent/=nbin;
  binerror/=nbin;
  if (is_north_or_south)
  {
    line = new TLine(0,bincontent,100,bincontent);
  }
  else
  {
    line = new TLine(-100,bincontent,0,bincontent);
  }
}

void draw1D_z_from3D_sim_closure_weighted(TString phiType="central")
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

std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
for (const auto& selectR : selectRs)
{

oneDHist hists_CPM_sim_acts, hists_CPM_sim_genfit, hists_CPM_sim_truth;
draw1Dmap(Form("../output/sim_run%d_acts_weighted/sim_run%d_acts_weighted_maxdca0p2_B3_average_correction_histograms.root", run, run), Form("CPM_sim_reco_acts"), selectR, selectPhi, hists_CPM_sim_acts, 2);
draw1Dmap(Form("../output/sim_run%d_genfit_weighted/sim_run%d_genfit_weighted_maxdca0p2_B3_average_correction_histograms.root", run, run), Form("CPM_sim_reco_genfit"), selectR, selectPhi, hists_CPM_sim_genfit, 3);
draw1Dmap(Form("../output/sim_run%d_truth_weighted/sim_run%d_truth_weighted_maxdca0p2_B3_average_correction_histograms.root", run, run), Form("CPM_sim_reco_truth"), selectR, selectPhi, hists_CPM_sim_truth, 4);

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
std::pair<double,double> yrange_P = {-0.01, 0.01}; //SetCommonYRange(hists_P);
std::pair<double,double> yrange_R = {-0.5, 0.5}; //SetCommonYRange(hists_R);
std::pair<double,double> yrange_Z = {-1, 1}; //SetCommonYRange(hists_Z);

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

TFile* ofile_N = new TFile(Form("Rootfiles/hist_n_1Dz_R%d_weighted_%s.root",(int)selectR,phiType.Data()),"recreate");
ofile_N->cd();
hists_CPM_sim_acts.h_N_pos->Write();
hists_CPM_sim_acts.h_N_neg->Write();
hists_CPM_sim_genfit.h_N_pos->Write();
hists_CPM_sim_genfit.h_N_neg->Write();
hists_CPM_sim_truth.h_N_pos->Write();
hists_CPM_sim_truth.h_N_neg->Write();
ofile_N->Write();

TFile* ofile_R = new TFile(Form("Rootfiles/hist_dr_1Dz_R%d_weighted_%s.root",(int)selectR,phiType.Data()),"recreate");
ofile_R->cd();
hists_CPM_sim_acts.h_R_pos->Write();
hists_CPM_sim_acts.h_R_neg->Write();
hists_CPM_sim_genfit.h_R_pos->Write();
hists_CPM_sim_genfit.h_R_neg->Write();
hists_CPM_sim_truth.h_R_pos->Write();
hists_CPM_sim_truth.h_R_neg->Write();
ofile_R->Write();

TFile* ofile_P = new TFile(Form("Rootfiles/hist_dphi_1Dz_R%d_weighted_%s.root",(int)selectR,phiType.Data()),"recreate");
ofile_P->cd();
hists_CPM_sim_acts.h_P_pos->Write();
hists_CPM_sim_acts.h_P_neg->Write();
hists_CPM_sim_genfit.h_P_pos->Write();
hists_CPM_sim_genfit.h_P_neg->Write();
hists_CPM_sim_truth.h_P_pos->Write();
hists_CPM_sim_truth.h_P_neg->Write();
ofile_P->Write();

TFile* ofile_Z = new TFile(Form("Rootfiles/hist_dz_1Dz_R%d_weighted_%s.root",(int)selectR,phiType.Data()),"recreate");
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
hists_CPM_sim_acts.h_P_pos->SetMinimum(-0.01); hists_CPM_sim_acts.h_P_pos->SetMaximum(0.01);
hists_CPM_sim_acts.h_P_pos->Draw("hist,same");
hists_CPM_sim_genfit.h_P_pos->Draw("hist,same");
hists_CPM_sim_truth.h_P_pos->Draw("hist,same");
if (l_nco_i_P) l_nco_i_P->Draw();
if (l_nco_o_P) l_nco_o_P->Draw();
l_nci_i_P->Draw();
l_nci_o_P->Draw();
can->cd(3);
gPad->SetLogy(0);
hists_CPM_sim_acts.h_R_pos->SetMinimum(-0.5); hists_CPM_sim_acts.h_R_pos->SetMaximum(0.5);
hists_CPM_sim_acts.h_R_pos->Draw("hist,same");
hists_CPM_sim_genfit.h_R_pos->Draw("hist,same");
hists_CPM_sim_truth.h_R_pos->Draw("hist,same");
if (l_nco_i_R) l_nco_i_R->Draw();
if (l_nco_o_R) l_nco_o_R->Draw();
l_nci_i_R->Draw();
l_nci_o_R->Draw();
can->cd(4);
gPad->SetLogy(0);
hists_CPM_sim_acts.h_Z_pos->SetMinimum(-1); hists_CPM_sim_acts.h_Z_pos->SetMaximum(1);
hists_CPM_sim_acts.h_Z_pos->Draw("hist,same");
hists_CPM_sim_genfit.h_Z_pos->Draw("hist,same");
hists_CPM_sim_truth.h_Z_pos->Draw("hist,same");
if (l_nco_i_Z) l_nco_i_Z->Draw();
if (l_nco_o_Z) l_nco_o_Z->Draw();
l_nci_i_Z->Draw();
l_nci_o_Z->Draw();
can->cd(5);
gPad->SetLogy(1);
hists_CPM_sim_acts.h_N_neg->Draw("hist,same");
hists_CPM_sim_genfit.h_N_neg->Draw("hist,same");
hists_CPM_sim_truth.h_N_neg->Draw("hist,same");
can->cd(6);
hists_CPM_sim_acts.h_P_neg->SetMinimum(-0.01); hists_CPM_sim_acts.h_P_neg->SetMaximum(0.01);
hists_CPM_sim_acts.h_P_neg->Draw("hist,same");
hists_CPM_sim_genfit.h_P_neg->Draw("hist,same");
hists_CPM_sim_truth.h_P_neg->Draw("hist,same");
if (l_sco_i_P) l_sco_i_P->Draw();
if (l_sco_o_P) l_sco_o_P->Draw();
l_sci_i_P->Draw();
l_sci_o_P->Draw();
can->cd(7);
gPad->SetLogy(0);
hists_CPM_sim_acts.h_R_neg->SetMinimum(-0.5); hists_CPM_sim_acts.h_R_neg->SetMaximum(0.5);
hists_CPM_sim_acts.h_R_neg->Draw("hist,same");
hists_CPM_sim_genfit.h_R_neg->Draw("hist,same");
hists_CPM_sim_truth.h_R_neg->Draw("hist,same");
if (l_sco_i_R) l_sco_i_R->Draw();
if (l_sco_o_R) l_sco_o_R->Draw();
l_sci_i_R->Draw();
l_sci_o_R->Draw();
can->cd(8);
gPad->SetLogy(0);
hists_CPM_sim_acts.h_Z_neg->SetMinimum(-1); hists_CPM_sim_acts.h_Z_neg->SetMaximum(1);
hists_CPM_sim_acts.h_Z_neg->Draw("hist,same");
hists_CPM_sim_genfit.h_Z_neg->Draw("hist,same");
hists_CPM_sim_truth.h_Z_neg->Draw("hist,same");
if (l_sco_i_Z) l_sco_i_Z->Draw();
if (l_sco_o_Z) l_sco_o_Z->Draw();
l_sci_i_Z->Draw();
l_sci_o_Z->Draw();

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsZ_from3D_atR%d_weighted_%s.pdf",(int)selectR,phiType.Data()));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->SetHeader("CPM weighted");
legend->AddEntry(hists_CPM_sim_acts.h_P_pos, Form("ACTS"), "F");
legend->AddEntry(hists_CPM_sim_genfit.h_P_pos, Form("GENFIT"), "F");
legend->AddEntry(hists_CPM_sim_truth.h_P_pos, Form("TRUTH"), "F");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs(Form("figure/resid_vsZ_from3D_leg_weighted_%s.pdf",phiType.Data()));

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

void draw1Dmap(TString filename, TString tag, double selectR, float selectPhi,
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

  hists1D.h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N @ R=%d cm;Z (cm);N",(int)selectR),hists3D.h_N_prz_pos->GetZaxis()->GetNbins(),hists3D.h_N_prz_pos->GetZaxis()->GetXmin(),hists3D.h_N_prz_pos->GetZaxis()->GetXmax());
  hists1D.h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),hists3D.h_R_prz_pos->GetZaxis()->GetNbins(),hists3D.h_R_prz_pos->GetZaxis()->GetXmin(),hists3D.h_R_prz_pos->GetZaxis()->GetXmax());
  hists1D.h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),hists3D.h_P_prz_pos->GetZaxis()->GetNbins(),hists3D.h_P_prz_pos->GetZaxis()->GetXmin(),hists3D.h_P_prz_pos->GetZaxis()->GetXmax());
  hists1D.h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),hists3D.h_Z_prz_pos->GetZaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetZaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetZaxis()->GetXmax());
  hists1D.h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N @ R=%d cm;Z (cm);N",(int)selectR),hists3D.h_N_prz_neg->GetZaxis()->GetNbins(),hists3D.h_N_prz_neg->GetZaxis()->GetXmin(),hists3D.h_N_prz_neg->GetZaxis()->GetXmax());
  hists1D.h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),hists3D.h_R_prz_neg->GetZaxis()->GetNbins(),hists3D.h_R_prz_neg->GetZaxis()->GetXmin(),hists3D.h_R_prz_neg->GetZaxis()->GetXmax());
  hists1D.h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),hists3D.h_P_prz_neg->GetZaxis()->GetNbins(),hists3D.h_P_prz_neg->GetZaxis()->GetXmin(),hists3D.h_P_prz_neg->GetZaxis()->GetXmax());
  hists1D.h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),hists3D.h_Z_prz_neg->GetZaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetZaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetZaxis()->GetXmax());

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

  plot1D_Rbin(hists3D.h_N_prz_pos,hists1D.h_N_pos,selectR, selectPhi, false);
  plot1D_Rbin(hists3D.h_R_prz_pos,hists1D.h_R_pos,selectR, selectPhi, false);
  plot1D_Rbin(hists3D.h_P_prz_pos,hists1D.h_P_pos,selectR, selectPhi, convert_RP_2_P);
  plot1D_Rbin(hists3D.h_Z_prz_pos,hists1D.h_Z_pos,selectR, selectPhi, false);
  plot1D_Rbin(hists3D.h_N_prz_neg,hists1D.h_N_neg,selectR, selectPhi, false);
  plot1D_Rbin(hists3D.h_R_prz_neg,hists1D.h_R_neg,selectR, selectPhi, false);
  plot1D_Rbin(hists3D.h_P_prz_neg,hists1D.h_P_neg,selectR, selectPhi, convert_RP_2_P);
  plot1D_Rbin(hists3D.h_Z_prz_neg,hists1D.h_Z_neg,selectR, selectPhi, false);

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
  int color)
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

void draw1Dmap_2D(TString filename, TString tag, double selectR,
  TLine*& l_R_pos, TLine*& l_R_neg,
  TLine*& l_P_pos, TLine*& l_P_neg,
  TLine*& l_Z_pos, TLine*& l_Z_neg,
  int color)
{
  TFile* file_2D_map = new TFile(filename,"");
  if (!file_2D_map || file_2D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  TH2* h_R_pr_pos = (TH2*) file_2D_map->Get("hIntDistortionR_posz");
  TH2* h_P_pr_pos = (TH2*) file_2D_map->Get("hIntDistortionP_posz");
  TH2* h_Z_pr_pos = (TH2*) file_2D_map->Get("hIntDistortionZ_posz");
  TH2* h_R_pr_neg = (TH2*) file_2D_map->Get("hIntDistortionR_negz");
  TH2* h_P_pr_neg = (TH2*) file_2D_map->Get("hIntDistortionP_negz");
  TH2* h_Z_pr_neg = (TH2*) file_2D_map->Get("hIntDistortionZ_negz");

  if (!h_R_pr_pos || !h_P_pr_pos || !h_Z_pr_pos ||
      !h_R_pr_neg || !h_P_pr_neg || !h_Z_pr_neg) {
      std::cerr << "Error: missing 2D histograms in file" << std::endl;
      delete file_2D_map;
      return;
  }

  l_R_pos = new TLine();
  l_P_pos = new TLine();
  l_Z_pos = new TLine();
  l_R_neg = new TLine();
  l_P_neg = new TLine();
  l_Z_neg = new TLine();

  plot1D_line(h_R_pr_pos,l_R_pos,selectR,true);
  plot1D_line(h_P_pr_pos,l_P_pos,selectR,true);
  plot1D_line(h_Z_pr_pos,l_Z_pos,selectR,true);
  plot1D_line(h_R_pr_neg,l_R_neg,selectR,false);
  plot1D_line(h_P_pr_neg,l_P_neg,selectR,false);
  plot1D_line(h_Z_pr_neg,l_Z_neg,selectR,false);

  l_R_pos->SetLineColor(color); l_R_pos->SetLineWidth(1);
  l_P_pos->SetLineColor(color); l_P_pos->SetLineWidth(1);
  l_Z_pos->SetLineColor(color); l_Z_pos->SetLineWidth(1);
  l_R_neg->SetLineColor(color); l_R_neg->SetLineWidth(1);
  l_P_neg->SetLineColor(color); l_P_neg->SetLineWidth(1);
  l_Z_neg->SetLineColor(color); l_Z_neg->SetLineWidth(1);

  delete file_2D_map;

  return;
}
