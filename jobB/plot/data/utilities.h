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
    TH3* h_N_prz_pos; TH3* h_R_prz_pos; TH3* h_P_prz_pos; TH3* h_Z_prz_pos;
    TH3* h_N_prz_neg; TH3* h_R_prz_neg; TH3* h_P_prz_neg; TH3* h_Z_prz_neg;
};

struct twoDHist {
    TH2* h_N_pos; TH2* h_R_pos; TH2* h_P_pos; TH2* h_Z_pos;
    TH2* h_N_neg; TH2* h_R_neg; TH2* h_P_neg; TH2* h_Z_neg;
};

struct oneDHist {
    TH1* h_N_pos; TH1* h_R_pos; TH1* h_P_pos; TH1* h_Z_pos;
    TH1* h_N_neg; TH1* h_R_neg; TH1* h_P_neg; TH1* h_Z_neg;
};

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);
std::pair<double,double> SetCommonZRange(const std::vector<TH2*>& histograms);
void draw1Dmap_atZ(TString filename, TString tag, double selectZ, float selectPhi, oneDHist& hists1D, int color=1, bool convert_RP_2_P=false, int style=1, bool read_N=true);
void draw1Dmap_atR(TString filename, TString tag, double selectR, float selectPhi, oneDHist& hists1D, int color=1, bool convert_RP_2_P=false, int style=1, bool read_N=true);
void draw1Dmap_atRZ(TString filename, TString tag, double selectR, double selectZ, oneDHist& hists1D, int color=1, bool convert_RP_2_P=false, int style=1, bool read_N=true);
void draw2Dmap_atP(TString filename, TString tag, double selectPhi, twoDHist& hists2D, bool convert_RP_2_P=false, bool read_N=true);
void draw2Dmap_atZ(TString filename, TString tag, double selectZ, twoDHist& hists2D, bool convert_RP_2_P=false, bool read_N=true);
void draw2Dmap_atR(TString filename, TString tag, double selectR, twoDHist& hists2D, bool convert_RP_2_P=false, bool read_N=true);
void draw1Dmap_2D(TString filename, TString tag, double selectR, float selectPhi, oneDHist& hists1D, int color=1, bool convert_RP_2_P=false, int style=1, bool read_N=true);
void plot1D_Zbin(TH3* h3, TH1* h1, float z, float selectPhi, bool convert_RP_2_P=false);
void plot1D_Rbin(TH3* h3, TH1* h1, float r, float selectPhi, bool convert_RP_2_P=false);
void plot1D_RZbin(TH3* h3, TH1* h1, float r, float z, bool convert_RP_2_P=false);
void plot2D_RZplane(TH3* h3, TH2* h2, float phi, bool convert_RP_2_P=false);
void plot2D_PRplane(TH3* h3, TH2* h2, float z, bool convert_RP_2_P=false);
void plot2D_PZplane(TH3* h3, TH2* h2, float r, bool convert_RP_2_P=false);
void plot1D_line(TH2* h2, TLine*& line, float y, bool is_north_or_south);
void plot1D(TH2* h2, TH1* h1, float selectPhi);
double phi_type2val(TString phiType);
void get1Dhist(TString filename, TString tag, oneDHist& hists1D, int color = 1, int style = 1, bool read_N = true);
void get2Dhist(TString filename, TString tag, twoDHist& hists2D, bool read_N = true);
void AddOneDHistVectorToVectors(const std::vector<oneDHist>& histVector, std::vector<TH1*>& hists_P, std::vector<TH1*>& hists_R, std::vector<TH1*>& hists_Z);
void AddTwoDHistVectorToVectors(const std::vector<twoDHist>& histVector, std::vector<TH2*>& hists_P, std::vector<TH2*>& hists_R, std::vector<TH2*>& hists_Z);

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

    if (yMin<0) yMin *= 1.1;
    else if (yMin>0) yMin *= 0.9;

    if (yMax<0) yMax *= 0.9;
    else if (yMax>0) yMax *= 1.1;

    for (TH1* h : histograms) {
        if (h) {
            h->SetMinimum(yMin);
            h->SetMaximum(yMax);
        }
    }
    return {yMin, yMax};
}

std::pair<double,double> SetCommonZRange(const std::vector<TH2*>& histograms)
{
    Double_t zMin = TMath::Infinity();
    Double_t zMax = -TMath::Infinity();

    for (TH1* h : histograms) {
        if (!h) continue;
        zMin = TMath::Min(zMin, h->GetMinimum());
        zMax = TMath::Max(zMax, h->GetMaximum());
    }

    if (zMin == TMath::Infinity()) zMin = 0;
    if (zMax == -TMath::Infinity()) zMax = 1;

    if (zMin<0) zMin *= 1.1;
    else if (zMin>0) zMin *= 0.9;

    if (zMax<0) zMax *= 0.9;
    else if (zMax>0) zMax *= 1.1;

    for (TH2* h : histograms) {
        if (h) {
            h->SetMinimum(zMin);
            h->SetMaximum(zMax);
        }
    }
    return {zMin, zMax};
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

double phi_type2val(TString phiType)
{
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
    return std::numeric_limits<double>::quiet_NaN();
  }
  return selectPhi;
}

void plot1D_Zbin(TH3* h3, TH1* h1, float z, float selectPhi, bool convert_RP_2_P)
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

void plot1D_Rbin(TH3* h3, TH1* h1, float r, float selectPhi, bool convert_RP_2_P)
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

void plot1D_RZbin(TH3* h3, TH1* h1, float r, float z, bool convert_RP_2_P)
{
  int rbin = h3->GetYaxis()->FindBin(r);
  int zbin = h3->GetZaxis()->FindBin(z);
  double rcenter = h3->GetYaxis()->GetBinCenter(rbin);
  for (int i = 1; i <= h3->GetNbinsX(); i++)
  {
      double value = h3->GetBinContent(i, rbin, zbin);
      double error = h3->GetBinError(i, rbin, zbin);
      if (convert_RP_2_P)
      {
        value /= rcenter;
        error /= rcenter;
      }
      h1->SetBinContent(i, value);
      h1->SetBinError(i, error);
  }
}

void plot2D_RZplane(TH3* h3, TH2* h2, float phi, bool convert_RP_2_P)
{
  int xbin = h3->GetXaxis()->FindBin(phi); // p
  for (int ybin = 1; ybin <= h3->GetNbinsY(); ybin++) // r
  {
    for (int zbin = 1; zbin <= h3->GetNbinsZ(); zbin++) // z
    {
      double rcenter = h3->GetYaxis()->GetBinCenter(ybin);

      double value = h3->GetBinContent(xbin, ybin, zbin);
      double error = h3->GetBinError(xbin, ybin, zbin);
      if (convert_RP_2_P)
      {
        value /= rcenter;
        error /= rcenter;
      }
      h2->SetBinContent(ybin, zbin, value);
      h2->SetBinError(ybin, zbin, error);
    }
  }
}

void plot2D_PRplane(TH3* h3, TH2* h2, float z, bool convert_RP_2_P)
{
  int zbin = h3->GetZaxis()->FindBin(z); // z
  for (int xbin = 1; xbin <= h3->GetNbinsX(); xbin++) // p
  {
    for (int ybin = 1; ybin <= h3->GetNbinsY(); ybin++) // r
    {
      double rcenter = h3->GetYaxis()->GetBinCenter(ybin);

      double value = h3->GetBinContent(xbin, ybin, zbin);
      double error = h3->GetBinError(xbin, ybin, zbin);
      if (convert_RP_2_P)
      {
        value /= rcenter;
        error /= rcenter;
      }
      h2->SetBinContent(xbin, ybin, value);
      h2->SetBinError(xbin, ybin, error);
    }
  }
}

void plot2D_PZplane(TH3* h3, TH2* h2, float r, bool convert_RP_2_P)
{
  int ybin = h3->GetYaxis()->FindBin(r); // r
  double rcenter = h3->GetYaxis()->GetBinCenter(ybin);
  for (int xbin = 1; xbin <= h3->GetNbinsX(); xbin++) // p
  {
    for (int zbin = 1; zbin <= h3->GetNbinsZ(); zbin++) // z
    {
      double value = h3->GetBinContent(xbin, ybin, zbin);
      double error = h3->GetBinError(xbin, ybin, zbin);
      if (convert_RP_2_P)
      {
        value /= rcenter;
        error /= rcenter;
      }
      h2->SetBinContent(xbin, zbin, value);
      h2->SetBinError(xbin, zbin, error);
    }
  }
}

void plot1D(TH2* h2, TH1* h1, float selectPhi)
{
  int xbin = h2->GetXaxis()->FindBin(selectPhi);
  for (int i = 1; i <= h2->GetNbinsY(); i++)
  {
      float bincontent=h2->GetBinContent(xbin, i);
      float binerror=h2->GetBinError(xbin, i);
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

void draw1Dmap_atZ(
  TString filename, TString tag,
  double selectZ, float selectPhi,
  oneDHist& hists1D, int color,
  bool convert_RP_2_P, int style,
  bool read_N)
{
  TFile* file_3D_map = new TFile(filename,"");
  if (!file_3D_map || file_3D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  threeDHist hists3D;
  if (read_N) hists3D.h_N_prz_pos = (TH3*) file_3D_map->Get("hentries_posz");
  hists3D.h_R_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  hists3D.h_P_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  hists3D.h_Z_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  if (read_N) hists3D.h_N_prz_neg = (TH3*) file_3D_map->Get("hentries_negz");
  hists3D.h_R_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  hists3D.h_P_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  hists3D.h_Z_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");

  if ((read_N && !hists3D.h_N_prz_pos) || !hists3D.h_R_prz_pos || !hists3D.h_P_prz_pos || !hists3D.h_Z_prz_pos ||
      (read_N && !hists3D.h_N_prz_neg) || !hists3D.h_R_prz_neg || !hists3D.h_P_prz_neg || !hists3D.h_Z_prz_neg) {
      std::cerr << "Error: missing 3D histograms in file" << std::endl;
      delete file_3D_map;
      return;
  }

  if (read_N) hists1D.h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),hists3D.h_N_prz_pos->GetYaxis()->GetNbins(),hists3D.h_N_prz_pos->GetYaxis()->GetXmin(),hists3D.h_N_prz_pos->GetYaxis()->GetXmax());
  hists1D.h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),hists3D.h_R_prz_pos->GetYaxis()->GetNbins(),hists3D.h_R_prz_pos->GetYaxis()->GetXmin(),hists3D.h_R_prz_pos->GetYaxis()->GetXmax());
  hists1D.h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),hists3D.h_P_prz_pos->GetYaxis()->GetNbins(),hists3D.h_P_prz_pos->GetYaxis()->GetXmin(),hists3D.h_P_prz_pos->GetYaxis()->GetXmax());
  hists1D.h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),hists3D.h_Z_prz_pos->GetYaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetYaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetYaxis()->GetXmax());
  if (read_N) hists1D.h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),hists3D.h_N_prz_neg->GetYaxis()->GetNbins(),hists3D.h_N_prz_neg->GetYaxis()->GetXmin(),hists3D.h_N_prz_neg->GetYaxis()->GetXmax());
  hists1D.h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),hists3D.h_R_prz_neg->GetYaxis()->GetNbins(),hists3D.h_R_prz_neg->GetYaxis()->GetXmin(),hists3D.h_R_prz_neg->GetYaxis()->GetXmax());
  hists1D.h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ Z=-%d cm;R (cm);dphi (rad)",(int)selectZ),hists3D.h_P_prz_neg->GetYaxis()->GetNbins(),hists3D.h_P_prz_neg->GetYaxis()->GetXmin(),hists3D.h_P_prz_neg->GetYaxis()->GetXmax());
  hists1D.h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),hists3D.h_Z_prz_neg->GetYaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetYaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetYaxis()->GetXmax());

  // do not save in the file_3D_map
  // directly saved in the stack memory
  if (read_N) hists1D.h_N_pos->SetDirectory(0);
  hists1D.h_R_pos->SetDirectory(0);
  hists1D.h_P_pos->SetDirectory(0);
  hists1D.h_Z_pos->SetDirectory(0);
  if (read_N) hists1D.h_N_neg->SetDirectory(0);
  hists1D.h_R_neg->SetDirectory(0);
  hists1D.h_P_neg->SetDirectory(0);
  hists1D.h_Z_neg->SetDirectory(0);

  if (read_N) plot1D_Zbin(hists3D.h_N_prz_pos,hists1D.h_N_pos,selectZ, selectPhi, false);
  plot1D_Zbin(hists3D.h_R_prz_pos,hists1D.h_R_pos,selectZ, selectPhi, false);
  plot1D_Zbin(hists3D.h_P_prz_pos,hists1D.h_P_pos,selectZ, selectPhi, convert_RP_2_P);
  plot1D_Zbin(hists3D.h_Z_prz_pos,hists1D.h_Z_pos,selectZ, selectPhi, false);
  if (read_N) plot1D_Zbin(hists3D.h_N_prz_neg,hists1D.h_N_neg,-selectZ, selectPhi, false);
  plot1D_Zbin(hists3D.h_R_prz_neg,hists1D.h_R_neg,-selectZ, selectPhi, false);
  plot1D_Zbin(hists3D.h_P_prz_neg,hists1D.h_P_neg,-selectZ, selectPhi, convert_RP_2_P);
  plot1D_Zbin(hists3D.h_Z_prz_neg,hists1D.h_Z_neg,-selectZ, selectPhi, false);

  if (read_N) {hists1D.h_N_pos->SetLineColor(color); hists1D.h_N_pos->SetLineWidth(1); hists1D.h_N_pos->SetFillColor(0); hists1D.h_N_pos->SetMarkerColor(color); hists1D.h_N_pos->SetLineStyle(style);}
  hists1D.h_R_pos->SetLineColor(color); hists1D.h_R_pos->SetLineWidth(1); hists1D.h_R_pos->SetFillColor(0); hists1D.h_R_pos->SetMarkerColor(color); hists1D.h_R_pos->SetLineStyle(style);
  hists1D.h_P_pos->SetLineColor(color); hists1D.h_P_pos->SetLineWidth(1); hists1D.h_P_pos->SetFillColor(0); hists1D.h_P_pos->SetMarkerColor(color); hists1D.h_P_pos->SetLineStyle(style);
  hists1D.h_Z_pos->SetLineColor(color); hists1D.h_Z_pos->SetLineWidth(1); hists1D.h_Z_pos->SetFillColor(0); hists1D.h_Z_pos->SetMarkerColor(color); hists1D.h_Z_pos->SetLineStyle(style);
  if (read_N) {hists1D.h_N_neg->SetLineColor(color); hists1D.h_N_neg->SetLineWidth(1); hists1D.h_N_neg->SetFillColor(0); hists1D.h_N_neg->SetMarkerColor(color); hists1D.h_N_neg->SetLineStyle(style);}
  hists1D.h_R_neg->SetLineColor(color); hists1D.h_R_neg->SetLineWidth(1); hists1D.h_R_neg->SetFillColor(0); hists1D.h_R_neg->SetMarkerColor(color); hists1D.h_R_neg->SetLineStyle(style);
  hists1D.h_P_neg->SetLineColor(color); hists1D.h_P_neg->SetLineWidth(1); hists1D.h_P_neg->SetFillColor(0); hists1D.h_P_neg->SetMarkerColor(color); hists1D.h_P_neg->SetLineStyle(style);
  hists1D.h_Z_neg->SetLineColor(color); hists1D.h_Z_neg->SetLineWidth(1); hists1D.h_Z_neg->SetFillColor(0); hists1D.h_Z_neg->SetMarkerColor(color); hists1D.h_Z_neg->SetLineStyle(style);

  delete file_3D_map;

  return;
}

void draw1Dmap_atR(
  TString filename, TString tag,
  double selectR, float selectPhi,
  oneDHist& hists1D, int color,
  bool convert_RP_2_P, int style,
  bool read_N)
{
  TFile* file_3D_map = new TFile(filename,"");
  if (!file_3D_map || file_3D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  threeDHist hists3D;
  if (read_N) hists3D.h_N_prz_pos = (TH3*) file_3D_map->Get("hentries_posz");
  hists3D.h_R_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  hists3D.h_P_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  hists3D.h_Z_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  if (read_N) hists3D.h_N_prz_neg = (TH3*) file_3D_map->Get("hentries_negz");
  hists3D.h_R_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  hists3D.h_P_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  hists3D.h_Z_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");

  if ((read_N && !hists3D.h_N_prz_pos) || !hists3D.h_R_prz_pos || !hists3D.h_P_prz_pos || !hists3D.h_Z_prz_pos ||
      (read_N && !hists3D.h_N_prz_neg) || !hists3D.h_R_prz_neg || !hists3D.h_P_prz_neg || !hists3D.h_Z_prz_neg) {
      std::cerr << "Error: missing 3D histograms in file" << std::endl;
      delete file_3D_map;
      return;
  }

  if (read_N) hists1D.h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N @ R=%d cm;Z (cm);N",(int)selectR),hists3D.h_N_prz_pos->GetZaxis()->GetNbins(),hists3D.h_N_prz_pos->GetZaxis()->GetXmin(),hists3D.h_N_prz_pos->GetZaxis()->GetXmax());
  hists1D.h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),hists3D.h_R_prz_pos->GetZaxis()->GetNbins(),hists3D.h_R_prz_pos->GetZaxis()->GetXmin(),hists3D.h_R_prz_pos->GetZaxis()->GetXmax());
  hists1D.h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),hists3D.h_P_prz_pos->GetZaxis()->GetNbins(),hists3D.h_P_prz_pos->GetZaxis()->GetXmin(),hists3D.h_P_prz_pos->GetZaxis()->GetXmax());
  hists1D.h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),hists3D.h_Z_prz_pos->GetZaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetZaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetZaxis()->GetXmax());
  if (read_N) hists1D.h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N @ R=%d cm;Z (cm);N",(int)selectR),hists3D.h_N_prz_neg->GetZaxis()->GetNbins(),hists3D.h_N_prz_neg->GetZaxis()->GetXmin(),hists3D.h_N_prz_neg->GetZaxis()->GetXmax());
  hists1D.h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),hists3D.h_R_prz_neg->GetZaxis()->GetNbins(),hists3D.h_R_prz_neg->GetZaxis()->GetXmin(),hists3D.h_R_prz_neg->GetZaxis()->GetXmax());
  hists1D.h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),hists3D.h_P_prz_neg->GetZaxis()->GetNbins(),hists3D.h_P_prz_neg->GetZaxis()->GetXmin(),hists3D.h_P_prz_neg->GetZaxis()->GetXmax());
  hists1D.h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),hists3D.h_Z_prz_neg->GetZaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetZaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetZaxis()->GetXmax());

  // do not save in the file_3D_map
  // directly saved in the stack memory
  if (read_N) hists1D.h_N_pos->SetDirectory(0);
  hists1D.h_R_pos->SetDirectory(0);
  hists1D.h_P_pos->SetDirectory(0);
  hists1D.h_Z_pos->SetDirectory(0);
  if (read_N) hists1D.h_N_neg->SetDirectory(0);
  hists1D.h_R_neg->SetDirectory(0);
  hists1D.h_P_neg->SetDirectory(0);
  hists1D.h_Z_neg->SetDirectory(0);

  if (read_N) plot1D_Rbin(hists3D.h_N_prz_pos,hists1D.h_N_pos,selectR, selectPhi, false);
  plot1D_Rbin(hists3D.h_R_prz_pos,hists1D.h_R_pos,selectR, selectPhi, false);
  plot1D_Rbin(hists3D.h_P_prz_pos,hists1D.h_P_pos,selectR, selectPhi, convert_RP_2_P);
  plot1D_Rbin(hists3D.h_Z_prz_pos,hists1D.h_Z_pos,selectR, selectPhi, false);
  if (read_N) plot1D_Rbin(hists3D.h_N_prz_neg,hists1D.h_N_neg,selectR, selectPhi, false);
  plot1D_Rbin(hists3D.h_R_prz_neg,hists1D.h_R_neg,selectR, selectPhi, false);
  plot1D_Rbin(hists3D.h_P_prz_neg,hists1D.h_P_neg,selectR, selectPhi, convert_RP_2_P);
  plot1D_Rbin(hists3D.h_Z_prz_neg,hists1D.h_Z_neg,selectR, selectPhi, false);

  if (read_N) {hists1D.h_N_pos->SetLineColor(color); hists1D.h_N_pos->SetLineWidth(1); hists1D.h_N_pos->SetFillColor(0); hists1D.h_N_pos->SetMarkerColor(color); hists1D.h_N_pos->SetLineStyle(style);}
  hists1D.h_R_pos->SetLineColor(color); hists1D.h_R_pos->SetLineWidth(1); hists1D.h_R_pos->SetFillColor(0); hists1D.h_R_pos->SetMarkerColor(color); hists1D.h_R_pos->SetLineStyle(style);
  hists1D.h_P_pos->SetLineColor(color); hists1D.h_P_pos->SetLineWidth(1); hists1D.h_P_pos->SetFillColor(0); hists1D.h_P_pos->SetMarkerColor(color); hists1D.h_P_pos->SetLineStyle(style);
  hists1D.h_Z_pos->SetLineColor(color); hists1D.h_Z_pos->SetLineWidth(1); hists1D.h_Z_pos->SetFillColor(0); hists1D.h_Z_pos->SetMarkerColor(color); hists1D.h_Z_pos->SetLineStyle(style);
  if (read_N) {hists1D.h_N_neg->SetLineColor(color); hists1D.h_N_neg->SetLineWidth(1); hists1D.h_N_neg->SetFillColor(0); hists1D.h_N_neg->SetMarkerColor(color); hists1D.h_N_neg->SetLineStyle(style);}
  hists1D.h_R_neg->SetLineColor(color); hists1D.h_R_neg->SetLineWidth(1); hists1D.h_R_neg->SetFillColor(0); hists1D.h_R_neg->SetMarkerColor(color); hists1D.h_R_neg->SetLineStyle(style);
  hists1D.h_P_neg->SetLineColor(color); hists1D.h_P_neg->SetLineWidth(1); hists1D.h_P_neg->SetFillColor(0); hists1D.h_P_neg->SetMarkerColor(color); hists1D.h_P_neg->SetLineStyle(style);
  hists1D.h_Z_neg->SetLineColor(color); hists1D.h_Z_neg->SetLineWidth(1); hists1D.h_Z_neg->SetFillColor(0); hists1D.h_Z_neg->SetMarkerColor(color); hists1D.h_Z_neg->SetLineStyle(style);

  delete file_3D_map;

  return;
}

void draw1Dmap_atRZ(
  TString filename, TString tag,
  double selectR, double selectZ,
  oneDHist& hists1D, int color,
  bool convert_RP_2_P, int style,
  bool read_N)
{
  TFile* file_3D_map = new TFile(filename,"");
  if (!file_3D_map || file_3D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  threeDHist hists3D;
  if (read_N) hists3D.h_N_prz_pos = (TH3*) file_3D_map->Get("hentries_posz");
  hists3D.h_R_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  hists3D.h_P_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  hists3D.h_Z_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  if (read_N) hists3D.h_N_prz_neg = (TH3*) file_3D_map->Get("hentries_negz");
  hists3D.h_R_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  hists3D.h_P_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  hists3D.h_Z_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");

  if ((read_N && !hists3D.h_N_prz_pos) || !hists3D.h_R_prz_pos || !hists3D.h_P_prz_pos || !hists3D.h_Z_prz_pos ||
      (read_N && !hists3D.h_N_prz_neg) || !hists3D.h_R_prz_neg || !hists3D.h_P_prz_neg || !hists3D.h_Z_prz_neg) {
      std::cerr << "Error: missing 3D histograms in file" << std::endl;
      delete file_3D_map;
      return;
  }

  if (read_N) hists1D.h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N @ R=%d cm, Z=%d cm;#phi (rad);N",(int)selectR, (int)selectZ),hists3D.h_N_prz_pos->GetXaxis()->GetNbins(),hists3D.h_N_prz_pos->GetXaxis()->GetXmin(),hists3D.h_N_prz_pos->GetXaxis()->GetXmax());
  hists1D.h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ R=%d cm, Z=%d cm;#phi (rad);dR (cm)",(int)selectR, (int)selectZ),hists3D.h_R_prz_pos->GetXaxis()->GetNbins(),hists3D.h_R_prz_pos->GetXaxis()->GetXmin(),hists3D.h_R_prz_pos->GetXaxis()->GetXmax());
  hists1D.h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ R=%d cm, Z=%d cm;#phi (rad);dphi (rad)",(int)selectR, (int)selectZ),hists3D.h_P_prz_pos->GetXaxis()->GetNbins(),hists3D.h_P_prz_pos->GetXaxis()->GetXmin(),hists3D.h_P_prz_pos->GetXaxis()->GetXmax());
  hists1D.h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz @ R=%d cm, Z=%d cm;#phi (rad);dz (cm)",(int)selectR, (int)selectZ),hists3D.h_Z_prz_pos->GetXaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetXaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetXaxis()->GetXmax());
  if (read_N) hists1D.h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N @ R=%d cm, Z=%d cm;#phi (rad);N",(int)selectR, (int)selectZ),hists3D.h_N_prz_neg->GetXaxis()->GetNbins(),hists3D.h_N_prz_neg->GetXaxis()->GetXmin(),hists3D.h_N_prz_neg->GetXaxis()->GetXmax());
  hists1D.h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ R=%d cm, Z=%d cm;#phi (rad);dR (cm)",(int)selectR, (int)selectZ),hists3D.h_R_prz_neg->GetXaxis()->GetNbins(),hists3D.h_R_prz_neg->GetXaxis()->GetXmin(),hists3D.h_R_prz_neg->GetXaxis()->GetXmax());
  hists1D.h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ R=%d cm, Z=%d cm;#phi (rad);dphi (rad)",(int)selectR, (int)selectZ),hists3D.h_P_prz_neg->GetXaxis()->GetNbins(),hists3D.h_P_prz_neg->GetXaxis()->GetXmin(),hists3D.h_P_prz_neg->GetXaxis()->GetXmax());
  hists1D.h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz @ R=%d cm, Z=%d cm;#phi (rad);dz (cm)",(int)selectR, (int)selectZ),hists3D.h_Z_prz_neg->GetXaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetXaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetXaxis()->GetXmax());

  // do not save in the file_3D_map
  // directly saved in the stack memory
  if (read_N) hists1D.h_N_pos->SetDirectory(0);
  hists1D.h_R_pos->SetDirectory(0);
  hists1D.h_P_pos->SetDirectory(0);
  hists1D.h_Z_pos->SetDirectory(0);
  if (read_N) hists1D.h_N_neg->SetDirectory(0);
  hists1D.h_R_neg->SetDirectory(0);
  hists1D.h_P_neg->SetDirectory(0);
  hists1D.h_Z_neg->SetDirectory(0);

  if (read_N) plot1D_RZbin(hists3D.h_N_prz_pos,hists1D.h_N_pos,selectR, selectZ, false);
  plot1D_RZbin(hists3D.h_R_prz_pos,hists1D.h_R_pos,selectR, selectZ, false);
  plot1D_RZbin(hists3D.h_P_prz_pos,hists1D.h_P_pos,selectR, selectZ, convert_RP_2_P);
  plot1D_RZbin(hists3D.h_Z_prz_pos,hists1D.h_Z_pos,selectR, selectZ, false);
  if (read_N) plot1D_RZbin(hists3D.h_N_prz_neg,hists1D.h_N_neg,selectR, -selectZ, false);
  plot1D_RZbin(hists3D.h_R_prz_neg,hists1D.h_R_neg,selectR, -selectZ, false);
  plot1D_RZbin(hists3D.h_P_prz_neg,hists1D.h_P_neg,selectR, -selectZ, convert_RP_2_P);
  plot1D_RZbin(hists3D.h_Z_prz_neg,hists1D.h_Z_neg,selectR, -selectZ, false);

  if (read_N) {hists1D.h_N_pos->SetLineColor(color); hists1D.h_N_pos->SetLineWidth(1); hists1D.h_N_pos->SetFillColor(0); hists1D.h_N_pos->SetMarkerColor(color); hists1D.h_N_pos->SetLineStyle(style);}
  hists1D.h_R_pos->SetLineColor(color); hists1D.h_R_pos->SetLineWidth(1); hists1D.h_R_pos->SetFillColor(0); hists1D.h_R_pos->SetMarkerColor(color); hists1D.h_R_pos->SetLineStyle(style);
  hists1D.h_P_pos->SetLineColor(color); hists1D.h_P_pos->SetLineWidth(1); hists1D.h_P_pos->SetFillColor(0); hists1D.h_P_pos->SetMarkerColor(color); hists1D.h_P_pos->SetLineStyle(style);
  hists1D.h_Z_pos->SetLineColor(color); hists1D.h_Z_pos->SetLineWidth(1); hists1D.h_Z_pos->SetFillColor(0); hists1D.h_Z_pos->SetMarkerColor(color); hists1D.h_Z_pos->SetLineStyle(style);
  if (read_N) {hists1D.h_N_neg->SetLineColor(color); hists1D.h_N_neg->SetLineWidth(1); hists1D.h_N_neg->SetFillColor(0); hists1D.h_N_neg->SetMarkerColor(color); hists1D.h_N_neg->SetLineStyle(style);}
  hists1D.h_R_neg->SetLineColor(color); hists1D.h_R_neg->SetLineWidth(1); hists1D.h_R_neg->SetFillColor(0); hists1D.h_R_neg->SetMarkerColor(color); hists1D.h_R_neg->SetLineStyle(style);
  hists1D.h_P_neg->SetLineColor(color); hists1D.h_P_neg->SetLineWidth(1); hists1D.h_P_neg->SetFillColor(0); hists1D.h_P_neg->SetMarkerColor(color); hists1D.h_P_neg->SetLineStyle(style);
  hists1D.h_Z_neg->SetLineColor(color); hists1D.h_Z_neg->SetLineWidth(1); hists1D.h_Z_neg->SetFillColor(0); hists1D.h_Z_neg->SetMarkerColor(color); hists1D.h_Z_neg->SetLineStyle(style);

  delete file_3D_map;

  return;
}

void draw2Dmap_atP(
  TString filename, TString tag,
  double selectPhi,
  twoDHist& hists2D,
  bool convert_RP_2_P,
  bool read_N)
{
  TFile* file_3D_map = new TFile(filename,"");
  if (!file_3D_map || file_3D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  threeDHist hists3D;
  if (read_N) hists3D.h_N_prz_pos = (TH3*) file_3D_map->Get("hentries_posz");
  hists3D.h_R_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  hists3D.h_P_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  hists3D.h_Z_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  if (read_N) hists3D.h_N_prz_neg = (TH3*) file_3D_map->Get("hentries_negz");
  hists3D.h_R_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  hists3D.h_P_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  hists3D.h_Z_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");

  if ((read_N && !hists3D.h_N_prz_pos) || !hists3D.h_R_prz_pos || !hists3D.h_P_prz_pos || !hists3D.h_Z_prz_pos ||
      (read_N && !hists3D.h_N_prz_neg) || !hists3D.h_R_prz_neg || !hists3D.h_P_prz_neg || !hists3D.h_Z_prz_neg) {
      std::cerr << "Error: missing 3D histograms in file" << std::endl;
      delete file_3D_map;
      return;
  }

  if (read_N) hists2D.h_N_pos = new TH2F(Form("hentries_posz_%s",tag.Data()),Form("N @ phi %.5f rad;R (cm);Z (cm);N", selectPhi),hists3D.h_N_prz_pos->GetYaxis()->GetNbins(),hists3D.h_N_prz_pos->GetYaxis()->GetXmin(),hists3D.h_N_prz_pos->GetYaxis()->GetXmax(),hists3D.h_N_prz_pos->GetZaxis()->GetNbins(),hists3D.h_N_prz_pos->GetZaxis()->GetXmin(),hists3D.h_N_prz_pos->GetZaxis()->GetXmax());
  hists2D.h_R_pos = new TH2F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ phi %.5f rad;R (cm);Z (cm);dR (cm)",selectPhi),hists3D.h_R_prz_pos->GetYaxis()->GetNbins(),hists3D.h_R_prz_pos->GetYaxis()->GetXmin(),hists3D.h_R_prz_pos->GetYaxis()->GetXmax(),hists3D.h_R_prz_pos->GetZaxis()->GetNbins(),hists3D.h_R_prz_pos->GetZaxis()->GetXmin(),hists3D.h_R_prz_pos->GetZaxis()->GetXmax());
  hists2D.h_P_pos = new TH2F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ phi %.5f rad;R (cm);Z (cm);dphi (rad)",selectPhi),hists3D.h_P_prz_pos->GetYaxis()->GetNbins(),hists3D.h_P_prz_pos->GetYaxis()->GetXmin(),hists3D.h_P_prz_pos->GetYaxis()->GetXmax(),hists3D.h_P_prz_pos->GetZaxis()->GetNbins(),hists3D.h_P_prz_pos->GetZaxis()->GetXmin(),hists3D.h_P_prz_pos->GetZaxis()->GetXmax());
  hists2D.h_Z_pos = new TH2F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz @ phi %.5f rad;R (cm);Z (cm);dz (cm)",selectPhi),hists3D.h_Z_prz_pos->GetYaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetYaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetYaxis()->GetXmax(),hists3D.h_Z_prz_pos->GetZaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetZaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetZaxis()->GetXmax());
  if (read_N) hists2D.h_N_neg = new TH2F(Form("hentries_negz_%s",tag.Data()),Form("N @ phi %.5f rad;R (cm);Z (cm);N",selectPhi),hists3D.h_N_prz_neg->GetYaxis()->GetNbins(),hists3D.h_N_prz_neg->GetYaxis()->GetXmin(),hists3D.h_N_prz_neg->GetYaxis()->GetXmax(),hists3D.h_N_prz_neg->GetZaxis()->GetNbins(),hists3D.h_N_prz_neg->GetZaxis()->GetXmin(),hists3D.h_N_prz_neg->GetZaxis()->GetXmax());
  hists2D.h_R_neg = new TH2F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ @ phi %.5f rad;R (cm);Z (cm);dR (cm)",selectPhi),hists3D.h_R_prz_neg->GetYaxis()->GetNbins(),hists3D.h_R_prz_neg->GetYaxis()->GetXmin(),hists3D.h_R_prz_neg->GetYaxis()->GetXmax(),hists3D.h_R_prz_neg->GetZaxis()->GetNbins(),hists3D.h_R_prz_neg->GetZaxis()->GetXmin(),hists3D.h_R_prz_neg->GetZaxis()->GetXmax());
  hists2D.h_P_neg = new TH2F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ phi %.5f rad;R (cm);Z (cm);dphi (rad)",selectPhi),hists3D.h_P_prz_neg->GetYaxis()->GetNbins(),hists3D.h_P_prz_neg->GetYaxis()->GetXmin(),hists3D.h_P_prz_neg->GetYaxis()->GetXmax(),hists3D.h_P_prz_neg->GetZaxis()->GetNbins(),hists3D.h_P_prz_neg->GetZaxis()->GetXmin(),hists3D.h_P_prz_neg->GetZaxis()->GetXmax());
  hists2D.h_Z_neg = new TH2F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz @ phi %.5f rad;R (cm);Z (cm);dz (cm)",selectPhi),hists3D.h_Z_prz_neg->GetYaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetYaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetYaxis()->GetXmax(),hists3D.h_Z_prz_neg->GetZaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetZaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetZaxis()->GetXmax());

  // do not save in the file_3D_map
  // directly saved in the stack memory
  if (read_N) hists2D.h_N_pos->SetDirectory(0);
  hists2D.h_R_pos->SetDirectory(0);
  hists2D.h_P_pos->SetDirectory(0);
  hists2D.h_Z_pos->SetDirectory(0);
  if (read_N) hists2D.h_N_neg->SetDirectory(0);
  hists2D.h_R_neg->SetDirectory(0);
  hists2D.h_P_neg->SetDirectory(0);
  hists2D.h_Z_neg->SetDirectory(0);

  if (read_N) plot2D_RZplane(hists3D.h_N_prz_pos, hists2D.h_N_pos,selectPhi, false);
  plot2D_RZplane(hists3D.h_R_prz_pos, hists2D.h_R_pos, selectPhi, false);
  plot2D_RZplane(hists3D.h_P_prz_pos, hists2D.h_P_pos, selectPhi, convert_RP_2_P);
  plot2D_RZplane(hists3D.h_Z_prz_pos, hists2D.h_Z_pos, selectPhi, false);
  if (read_N) plot2D_RZplane(hists3D.h_N_prz_neg, hists2D.h_N_neg, selectPhi, false);
  plot2D_RZplane(hists3D.h_R_prz_neg, hists2D.h_R_neg, selectPhi, false);
  plot2D_RZplane(hists3D.h_P_prz_neg, hists2D.h_P_neg, selectPhi, convert_RP_2_P);
  plot2D_RZplane(hists3D.h_Z_prz_neg, hists2D.h_Z_neg, selectPhi, false);
  delete file_3D_map;

  return;
}

void draw2Dmap_atZ(
  TString filename, TString tag,
  double selectZ,
  twoDHist& hists2D,
  bool convert_RP_2_P,
  bool read_N)
{
  TFile* file_3D_map = new TFile(filename,"");
  if (!file_3D_map || file_3D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  threeDHist hists3D;
  if (read_N) hists3D.h_N_prz_pos = (TH3*) file_3D_map->Get("hentries_posz");
  hists3D.h_R_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  hists3D.h_P_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  hists3D.h_Z_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  if (read_N) hists3D.h_N_prz_neg = (TH3*) file_3D_map->Get("hentries_negz");
  hists3D.h_R_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  hists3D.h_P_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  hists3D.h_Z_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");

  if ((read_N && !hists3D.h_N_prz_pos) || !hists3D.h_R_prz_pos || !hists3D.h_P_prz_pos || !hists3D.h_Z_prz_pos ||
      (read_N && !hists3D.h_N_prz_neg) || !hists3D.h_R_prz_neg || !hists3D.h_P_prz_neg || !hists3D.h_Z_prz_neg) {
      std::cerr << "Error: missing 3D histograms in file" << std::endl;
      delete file_3D_map;
      return;
  }

  if (read_N) hists2D.h_N_pos = new TH2F(Form("hentries_posz_%s",tag.Data()),Form("N @ Z %d cm;#phi (rad);R (cm);N", (int)selectZ),hists3D.h_N_prz_pos->GetXaxis()->GetNbins(),hists3D.h_N_prz_pos->GetXaxis()->GetXmin(),hists3D.h_N_prz_pos->GetXaxis()->GetXmax(),hists3D.h_N_prz_pos->GetYaxis()->GetNbins(),hists3D.h_N_prz_pos->GetYaxis()->GetXmin(),hists3D.h_N_prz_pos->GetYaxis()->GetXmax());
  hists2D.h_R_pos = new TH2F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ Z %d cm;#phi (rad);R (cm);dR (cm)", (int)selectZ),hists3D.h_R_prz_pos->GetXaxis()->GetNbins(),hists3D.h_R_prz_pos->GetXaxis()->GetXmin(),hists3D.h_R_prz_pos->GetXaxis()->GetXmax(),hists3D.h_R_prz_pos->GetYaxis()->GetNbins(),hists3D.h_R_prz_pos->GetYaxis()->GetXmin(),hists3D.h_R_prz_pos->GetYaxis()->GetXmax());
  hists2D.h_P_pos = new TH2F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ Z %d cm;#phi (rad);R (cm);dphi (rad)", (int)selectZ),hists3D.h_P_prz_pos->GetXaxis()->GetNbins(),hists3D.h_P_prz_pos->GetXaxis()->GetXmin(),hists3D.h_P_prz_pos->GetXaxis()->GetXmax(),hists3D.h_P_prz_pos->GetYaxis()->GetNbins(),hists3D.h_P_prz_pos->GetYaxis()->GetXmin(),hists3D.h_P_prz_pos->GetYaxis()->GetXmax());
  hists2D.h_Z_pos = new TH2F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dZ @ Z %d cm;#phi (rad);R (cm);dz (cm)", (int)selectZ),hists3D.h_Z_prz_pos->GetXaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetXaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetXaxis()->GetXmax(),hists3D.h_Z_prz_pos->GetYaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetYaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetYaxis()->GetXmax());
  if (read_N) hists2D.h_N_neg = new TH2F(Form("hentries_negz_%s",tag.Data()),Form("N @ Z %d cm;#phi (rad);R (cm);N", (int)selectZ),hists3D.h_N_prz_neg->GetXaxis()->GetNbins(),hists3D.h_N_prz_neg->GetXaxis()->GetXmin(),hists3D.h_N_prz_neg->GetXaxis()->GetXmax(),hists3D.h_N_prz_neg->GetYaxis()->GetNbins(),hists3D.h_N_prz_neg->GetYaxis()->GetXmin(),hists3D.h_N_prz_neg->GetYaxis()->GetXmax());
  hists2D.h_R_neg = new TH2F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ Z %d cm;#phi (rad);R (cm);dR (cm)", (int)selectZ),hists3D.h_R_prz_neg->GetXaxis()->GetNbins(),hists3D.h_R_prz_neg->GetXaxis()->GetXmin(),hists3D.h_R_prz_neg->GetXaxis()->GetXmax(),hists3D.h_R_prz_neg->GetYaxis()->GetNbins(),hists3D.h_R_prz_neg->GetYaxis()->GetXmin(),hists3D.h_R_prz_neg->GetYaxis()->GetXmax());
  hists2D.h_P_neg = new TH2F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ Z %d cm;#phi (rad);R (cm);dphi (rad)", (int)selectZ),hists3D.h_P_prz_neg->GetXaxis()->GetNbins(),hists3D.h_P_prz_neg->GetXaxis()->GetXmin(),hists3D.h_P_prz_neg->GetXaxis()->GetXmax(),hists3D.h_P_prz_neg->GetYaxis()->GetNbins(),hists3D.h_P_prz_neg->GetYaxis()->GetXmin(),hists3D.h_P_prz_neg->GetYaxis()->GetXmax());
  hists2D.h_Z_neg = new TH2F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dZ @ Z %d cm;#phi (rad);R (cm);dz (cm)", (int)selectZ),hists3D.h_Z_prz_neg->GetXaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetXaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetXaxis()->GetXmax(),hists3D.h_Z_prz_neg->GetYaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetYaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetYaxis()->GetXmax());

  // do not save in the file_3D_map
  // directly saved in the stack memory
  if (read_N) hists2D.h_N_pos->SetDirectory(0);
  hists2D.h_R_pos->SetDirectory(0);
  hists2D.h_P_pos->SetDirectory(0);
  hists2D.h_Z_pos->SetDirectory(0);
  if (read_N) hists2D.h_N_neg->SetDirectory(0);
  hists2D.h_R_neg->SetDirectory(0);
  hists2D.h_P_neg->SetDirectory(0);
  hists2D.h_Z_neg->SetDirectory(0);

  if (read_N) plot2D_PRplane(hists3D.h_N_prz_pos, hists2D.h_N_pos,selectZ, false);
  plot2D_PRplane(hists3D.h_R_prz_pos, hists2D.h_R_pos, selectZ, false);
  plot2D_PRplane(hists3D.h_P_prz_pos, hists2D.h_P_pos, selectZ, convert_RP_2_P);
  plot2D_PRplane(hists3D.h_Z_prz_pos, hists2D.h_Z_pos, selectZ, false);
  if (read_N) plot2D_PRplane(hists3D.h_N_prz_neg, hists2D.h_N_neg, -selectZ, false);
  plot2D_PRplane(hists3D.h_R_prz_neg, hists2D.h_R_neg, -selectZ, false);
  plot2D_PRplane(hists3D.h_P_prz_neg, hists2D.h_P_neg, -selectZ, convert_RP_2_P);
  plot2D_PRplane(hists3D.h_Z_prz_neg, hists2D.h_Z_neg, -selectZ, false);
  delete file_3D_map;

  return;
}

void draw2Dmap_atR(
  TString filename, TString tag,
  double selectR,
  twoDHist& hists2D,
  bool convert_RP_2_P,
  bool read_N)
{
  TFile* file_3D_map = new TFile(filename,"");
  if (!file_3D_map || file_3D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  threeDHist hists3D;
  if (read_N) hists3D.h_N_prz_pos = (TH3*) file_3D_map->Get("hentries_posz");
  hists3D.h_R_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  hists3D.h_P_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  hists3D.h_Z_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  if (read_N) hists3D.h_N_prz_neg = (TH3*) file_3D_map->Get("hentries_negz");
  hists3D.h_R_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  hists3D.h_P_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  hists3D.h_Z_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");

  if ((read_N && !hists3D.h_N_prz_pos) || !hists3D.h_R_prz_pos || !hists3D.h_P_prz_pos || !hists3D.h_Z_prz_pos ||
      (read_N && !hists3D.h_N_prz_neg) || !hists3D.h_R_prz_neg || !hists3D.h_P_prz_neg || !hists3D.h_Z_prz_neg) {
      std::cerr << "Error: missing 3D histograms in file" << std::endl;
      delete file_3D_map;
      return;
  }

  if (read_N) hists2D.h_N_pos = new TH2F(Form("hentries_posz_%s",tag.Data()),Form("N @ R %d cm;#phi (rad);Z (cm);N", (int)selectR),hists3D.h_N_prz_pos->GetXaxis()->GetNbins(),hists3D.h_N_prz_pos->GetXaxis()->GetXmin(),hists3D.h_N_prz_pos->GetXaxis()->GetXmax(),hists3D.h_N_prz_pos->GetZaxis()->GetNbins(),hists3D.h_N_prz_pos->GetZaxis()->GetXmin(),hists3D.h_N_prz_pos->GetZaxis()->GetXmax());
  hists2D.h_R_pos = new TH2F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ R %d cm;#phi (rad);Z (cm);dR (cm)", (int)selectR),hists3D.h_R_prz_pos->GetXaxis()->GetNbins(),hists3D.h_R_prz_pos->GetXaxis()->GetXmin(),hists3D.h_R_prz_pos->GetXaxis()->GetXmax(),hists3D.h_R_prz_pos->GetZaxis()->GetNbins(),hists3D.h_R_prz_pos->GetZaxis()->GetXmin(),hists3D.h_R_prz_pos->GetZaxis()->GetXmax());
  hists2D.h_P_pos = new TH2F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ R %d cm;#phi (rad);Z (cm);dphi (rad)", (int)selectR),hists3D.h_P_prz_pos->GetXaxis()->GetNbins(),hists3D.h_P_prz_pos->GetXaxis()->GetXmin(),hists3D.h_P_prz_pos->GetXaxis()->GetXmax(),hists3D.h_P_prz_pos->GetZaxis()->GetNbins(),hists3D.h_P_prz_pos->GetZaxis()->GetXmin(),hists3D.h_P_prz_pos->GetZaxis()->GetXmax());
  hists2D.h_Z_pos = new TH2F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dZ @ R %d cm;#phi (rad);Z (cm);dz (cm)", (int)selectR),hists3D.h_Z_prz_pos->GetXaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetXaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetXaxis()->GetXmax(),hists3D.h_Z_prz_pos->GetZaxis()->GetNbins(),hists3D.h_Z_prz_pos->GetZaxis()->GetXmin(),hists3D.h_Z_prz_pos->GetZaxis()->GetXmax());
  if (read_N) hists2D.h_N_neg = new TH2F(Form("hentries_negz_%s",tag.Data()),Form("N @ R %d cm;#phi (rad);Z (cm);N", (int)selectR),hists3D.h_N_prz_neg->GetXaxis()->GetNbins(),hists3D.h_N_prz_neg->GetXaxis()->GetXmin(),hists3D.h_N_prz_neg->GetXaxis()->GetXmax(),hists3D.h_N_prz_neg->GetZaxis()->GetNbins(),hists3D.h_N_prz_neg->GetZaxis()->GetXmin(),hists3D.h_N_prz_neg->GetZaxis()->GetXmax());
  hists2D.h_R_neg = new TH2F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ R %d cm;#phi (rad);Z (cm);dR (cm)", (int)selectR),hists3D.h_R_prz_neg->GetXaxis()->GetNbins(),hists3D.h_R_prz_neg->GetXaxis()->GetXmin(),hists3D.h_R_prz_neg->GetXaxis()->GetXmax(),hists3D.h_R_prz_neg->GetZaxis()->GetNbins(),hists3D.h_R_prz_neg->GetZaxis()->GetXmin(),hists3D.h_R_prz_neg->GetZaxis()->GetXmax());
  hists2D.h_P_neg = new TH2F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ R %d cm;#phi (rad);Z (cm);dphi (rad)", (int)selectR),hists3D.h_P_prz_neg->GetXaxis()->GetNbins(),hists3D.h_P_prz_neg->GetXaxis()->GetXmin(),hists3D.h_P_prz_neg->GetXaxis()->GetXmax(),hists3D.h_P_prz_neg->GetZaxis()->GetNbins(),hists3D.h_P_prz_neg->GetZaxis()->GetXmin(),hists3D.h_P_prz_neg->GetZaxis()->GetXmax());
  hists2D.h_Z_neg = new TH2F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dZ @ R %d cm;#phi (rad);Z (cm);dz (cm)", (int)selectR),hists3D.h_Z_prz_neg->GetXaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetXaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetXaxis()->GetXmax(),hists3D.h_Z_prz_neg->GetZaxis()->GetNbins(),hists3D.h_Z_prz_neg->GetZaxis()->GetXmin(),hists3D.h_Z_prz_neg->GetZaxis()->GetXmax());

  // do not save in the file_3D_map
  // directly saved in the stack memory
  if (read_N) hists2D.h_N_pos->SetDirectory(0);
  hists2D.h_R_pos->SetDirectory(0);
  hists2D.h_P_pos->SetDirectory(0);
  hists2D.h_Z_pos->SetDirectory(0);
  if (read_N) hists2D.h_N_neg->SetDirectory(0);
  hists2D.h_R_neg->SetDirectory(0);
  hists2D.h_P_neg->SetDirectory(0);
  hists2D.h_Z_neg->SetDirectory(0);

  if (read_N) plot2D_PZplane(hists3D.h_N_prz_pos, hists2D.h_N_pos,selectR, false);
  plot2D_PZplane(hists3D.h_R_prz_pos, hists2D.h_R_pos, selectR, false);
  plot2D_PZplane(hists3D.h_P_prz_pos, hists2D.h_P_pos, selectR, convert_RP_2_P);
  plot2D_PZplane(hists3D.h_Z_prz_pos, hists2D.h_Z_pos, selectR, false);
  if (read_N) plot2D_PZplane(hists3D.h_N_prz_neg, hists2D.h_N_neg, selectR, false);
  plot2D_PZplane(hists3D.h_R_prz_neg, hists2D.h_R_neg, selectR, false);
  plot2D_PZplane(hists3D.h_P_prz_neg, hists2D.h_P_neg, selectR, convert_RP_2_P);
  plot2D_PZplane(hists3D.h_Z_prz_neg, hists2D.h_Z_neg, selectR, false);
  delete file_3D_map;

  return;
}

void draw1Dmap_2D(
  TString filename, TString tag,
  double selectR, float selectPhi,
  oneDHist& hists1D, int color,
  bool convert_RP_2_P, int style,
  bool read_N)
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

  //hists1D.h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N;R (cm);N"),h_N_pr_pos->GetYaxis()->GetNbins(),h_N_pr_pos->GetYaxis()->GetXmin(),h_N_pr_pos->GetYaxis()->GetXmax());
  hists1D.h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR;R (cm);dR (cm)"),h_R_pr_pos->GetYaxis()->GetNbins(),h_R_pr_pos->GetYaxis()->GetXmin(),h_R_pr_pos->GetYaxis()->GetXmax());
  hists1D.h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi;R (cm);dphi (rad)"),h_P_pr_pos->GetYaxis()->GetNbins(),h_P_pr_pos->GetYaxis()->GetXmin(),h_P_pr_pos->GetYaxis()->GetXmax());
  hists1D.h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz;R (cm);dz (cm)"),h_Z_pr_pos->GetYaxis()->GetNbins(),h_Z_pr_pos->GetYaxis()->GetXmin(),h_Z_pr_pos->GetYaxis()->GetXmax());
  //hists1D.h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N;R (cm);N"),h_N_pr_neg->GetYaxis()->GetNbins(),h_N_pr_neg->GetYaxis()->GetXmin(),h_N_pr_neg->GetYaxis()->GetXmax());
  hists1D.h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR;R (cm);dR (cm)"),h_R_pr_neg->GetYaxis()->GetNbins(),h_R_pr_neg->GetYaxis()->GetXmin(),h_R_pr_neg->GetYaxis()->GetXmax());
  hists1D.h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi;R (cm);dphi (rad)"),h_P_pr_neg->GetYaxis()->GetNbins(),h_P_pr_neg->GetYaxis()->GetXmin(),h_P_pr_neg->GetYaxis()->GetXmax());
  hists1D.h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz;R (cm);dz (cm)"),h_Z_pr_neg->GetYaxis()->GetNbins(),h_Z_pr_neg->GetYaxis()->GetXmin(),h_Z_pr_neg->GetYaxis()->GetXmax());

  // do not save in the file_2D_map
  // directly saved in the stack memory
  //hists1D.h_N_pos->SetDirectory(0);
  hists1D.h_R_pos->SetDirectory(0);
  hists1D.h_P_pos->SetDirectory(0);
  hists1D.h_Z_pos->SetDirectory(0);
  //hists1D.h_N_neg->SetDirectory(0);
  hists1D.h_R_neg->SetDirectory(0);
  hists1D.h_P_neg->SetDirectory(0);
  hists1D.h_Z_neg->SetDirectory(0);

  //plot1D(h_N_pr_pos, hists1D.h_N_pos, selectPhi);
  plot1D(h_R_pr_pos, hists1D.h_R_pos, selectPhi);
  plot1D(h_P_pr_pos, hists1D.h_P_pos, selectPhi);
  plot1D(h_Z_pr_pos, hists1D.h_Z_pos, selectPhi);
  //plot1D(h_N_pr_neg, hists1D.h_N_neg, selectPhi);
  plot1D(h_R_pr_neg, hists1D.h_R_neg, selectPhi);
  plot1D(h_P_pr_neg, hists1D.h_P_neg, selectPhi);
  plot1D(h_Z_pr_neg, hists1D.h_Z_neg, selectPhi);

  //hists1D.h_N_pos->SetLineColor(color); hists1D.h_N_pos->SetLineWidth(1); hists1D.h_N_pos->SetFillColor(0); hists1D.h_N_pos->SetMarkerColor(color); hists1D.h_N_pos->SetLineStyle(style);
  hists1D.h_R_pos->SetLineColor(color); hists1D.h_R_pos->SetLineWidth(1); hists1D.h_R_pos->SetFillColor(0); hists1D.h_R_pos->SetMarkerColor(color); hists1D.h_R_pos->SetLineStyle(style);
  hists1D.h_P_pos->SetLineColor(color); hists1D.h_P_pos->SetLineWidth(1); hists1D.h_P_pos->SetFillColor(0); hists1D.h_P_pos->SetMarkerColor(color); hists1D.h_P_pos->SetLineStyle(style);
  hists1D.h_Z_pos->SetLineColor(color); hists1D.h_Z_pos->SetLineWidth(1); hists1D.h_Z_pos->SetFillColor(0); hists1D.h_Z_pos->SetMarkerColor(color); hists1D.h_Z_pos->SetLineStyle(style);
  //hists1D.h_N_neg->SetLineColor(color); hists1D.h_N_neg->SetLineWidth(1); hists1D.h_N_neg->SetFillColor(0); hists1D.h_R_neg->SetMarkerColor(color); hists1D.h_N_neg->SetLineStyle(style);
  hists1D.h_R_neg->SetLineColor(color); hists1D.h_R_neg->SetLineWidth(1); hists1D.h_R_neg->SetFillColor(0); hists1D.h_R_neg->SetMarkerColor(color); hists1D.h_R_neg->SetLineStyle(style);
  hists1D.h_P_neg->SetLineColor(color); hists1D.h_P_neg->SetLineWidth(1); hists1D.h_P_neg->SetFillColor(0); hists1D.h_P_neg->SetMarkerColor(color); hists1D.h_P_neg->SetLineStyle(style);
  hists1D.h_Z_neg->SetLineColor(color); hists1D.h_Z_neg->SetLineWidth(1); hists1D.h_Z_neg->SetFillColor(0); hists1D.h_Z_neg->SetMarkerColor(color); hists1D.h_Z_neg->SetLineStyle(style);

  delete file_2D_map;

  return;
}

void get1Dhist(
  TString filename, TString tag,
  oneDHist& hists, int color,
  int style, bool read_N)
{
  TFile* file_map = new TFile(filename,"");
  if (!file_map || file_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  if (read_N) hists.h_N_pos = (TH1*) file_map->Get(Form("hentries_posz_%s",tag.Data()));
  hists.h_R_pos = (TH1*) file_map->Get(Form("hIntDistortionR_posz_%s",tag.Data()));
  hists.h_P_pos = (TH1*) file_map->Get(Form("hIntDistortionP_posz_%s",tag.Data()));
  hists.h_Z_pos = (TH1*) file_map->Get(Form("hIntDistortionZ_posz_%s",tag.Data()));
  if (read_N) hists.h_N_neg = (TH1*) file_map->Get(Form("hentries_negz_%s",tag.Data()));
  hists.h_R_neg = (TH1*) file_map->Get(Form("hIntDistortionR_negz_%s",tag.Data()));
  hists.h_P_neg = (TH1*) file_map->Get(Form("hIntDistortionP_negz_%s",tag.Data()));
  hists.h_Z_neg = (TH1*) file_map->Get(Form("hIntDistortionZ_negz_%s",tag.Data()));

  if ((read_N && !hists.h_N_pos) || !hists.h_R_pos || !hists.h_P_pos || !hists.h_Z_pos ||
      (read_N && !hists.h_N_neg) || !hists.h_R_neg || !hists.h_P_neg || !hists.h_Z_neg) {
      std::cerr << "Error: missing 1D histograms for tag " << tag << " in file " << filename << std::endl;
      delete file_map;
      return;
  }

  // do not save in the file_map
  // directly saved in the stack memory
  if (read_N) hists.h_N_pos->SetDirectory(0);
  hists.h_R_pos->SetDirectory(0);
  hists.h_P_pos->SetDirectory(0);
  hists.h_Z_pos->SetDirectory(0);
  if (read_N) hists.h_N_neg->SetDirectory(0);
  hists.h_R_neg->SetDirectory(0);
  hists.h_P_neg->SetDirectory(0);
  hists.h_Z_neg->SetDirectory(0);

std::vector<TH1*> histograms;
if (read_N)
{
  histograms = {
    hists.h_N_pos, hists.h_R_pos, hists.h_P_pos, hists.h_Z_pos,
    hists.h_N_neg, hists.h_R_neg, hists.h_P_neg, hists.h_Z_neg
  };
}
else
{
  histograms = {
    hists.h_R_pos, hists.h_P_pos, hists.h_Z_pos,
    hists.h_R_neg, hists.h_P_neg, hists.h_Z_neg
  };
}
    
  for (TH1* h : histograms) {
      if (h) {
          h->SetLineColor(color);
          h->SetLineWidth(1);
          h->SetFillColor(0);
          h->SetMarkerColor(color);
          h->SetLineStyle(style);
      }
  }

  delete file_map;

  return;
}

void get2Dhist(
  TString filename, TString tag,
  twoDHist& hists, bool read_N)
{
  TFile* file_map = new TFile(filename,"");
  if (!file_map || file_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  if (read_N) hists.h_N_pos = (TH2*) file_map->Get(Form("hentries_posz_%s",tag.Data()));
  hists.h_R_pos = (TH2*) file_map->Get(Form("hIntDistortionR_posz_%s",tag.Data()));
  hists.h_P_pos = (TH2*) file_map->Get(Form("hIntDistortionP_posz_%s",tag.Data()));
  hists.h_Z_pos = (TH2*) file_map->Get(Form("hIntDistortionZ_posz_%s",tag.Data()));
  if (read_N) hists.h_N_neg = (TH2*) file_map->Get(Form("hentries_negz_%s",tag.Data()));
  hists.h_R_neg = (TH2*) file_map->Get(Form("hIntDistortionR_negz_%s",tag.Data()));
  hists.h_P_neg = (TH2*) file_map->Get(Form("hIntDistortionP_negz_%s",tag.Data()));
  hists.h_Z_neg = (TH2*) file_map->Get(Form("hIntDistortionZ_negz_%s",tag.Data()));

  if ((read_N && !hists.h_N_pos) || !hists.h_R_pos || !hists.h_P_pos || !hists.h_Z_pos ||
      (read_N && !hists.h_N_neg) || !hists.h_R_neg || !hists.h_P_neg || !hists.h_Z_neg) {
      std::cerr << "Error: missing 1D histograms for tag " << tag << " in file " << filename << std::endl;
      delete file_map;
      return;
  }

  // do not save in the file_map
  // directly saved in the stack memory
  if (read_N) hists.h_N_pos->SetDirectory(0);
  hists.h_R_pos->SetDirectory(0);
  hists.h_P_pos->SetDirectory(0);
  hists.h_Z_pos->SetDirectory(0);
  if (read_N) hists.h_N_neg->SetDirectory(0);
  hists.h_R_neg->SetDirectory(0);
  hists.h_P_neg->SetDirectory(0);
  hists.h_Z_neg->SetDirectory(0);

std::vector<TH2*> histograms;
if (read_N)
{
  histograms = {
    hists.h_N_pos, hists.h_R_pos, hists.h_P_pos, hists.h_Z_pos,
    hists.h_N_neg, hists.h_R_neg, hists.h_P_neg, hists.h_Z_neg
  };
}
else
{
  histograms = {
    hists.h_R_pos, hists.h_P_pos, hists.h_Z_pos,
    hists.h_R_neg, hists.h_P_neg, hists.h_Z_neg
  };
}
    
  delete file_map;

  return;
}

void AddOneDHistVectorToVectors(
    const std::vector<oneDHist>& histVector,
    std::vector<TH1*>& hists_N,
    std::vector<TH1*>& hists_P,
    std::vector<TH1*>& hists_R,
    std::vector<TH1*>& hists_Z
) {
    hists_N.clear();
    hists_P.clear();
    hists_R.clear();
    hists_Z.clear();
    
    for (const auto& hist : histVector) {
        if (hist.h_N_neg) hists_N.push_back(hist.h_N_neg);
        if (hist.h_N_pos) hists_N.push_back(hist.h_N_pos);
        if (hist.h_P_neg) hists_P.push_back(hist.h_P_neg);
        if (hist.h_P_pos) hists_P.push_back(hist.h_P_pos);
        if (hist.h_R_neg) hists_R.push_back(hist.h_R_neg);
        if (hist.h_R_pos) hists_R.push_back(hist.h_R_pos);
        if (hist.h_Z_neg) hists_Z.push_back(hist.h_Z_neg);
        if (hist.h_Z_pos) hists_Z.push_back(hist.h_Z_pos);
    }
}

void AddTwoDHistVectorToVectors(
    const std::vector<twoDHist>& histVector,
    std::vector<TH2*>& hists_N,
    std::vector<TH2*>& hists_P,
    std::vector<TH2*>& hists_R,
    std::vector<TH2*>& hists_Z
) {
    hists_N.clear();
    hists_P.clear();
    hists_R.clear();
    hists_Z.clear();
    
    for (const auto& hist : histVector) {
        if (hist.h_N_neg) hists_N.push_back(hist.h_N_neg);
        if (hist.h_N_pos) hists_N.push_back(hist.h_N_pos);
        if (hist.h_P_neg) hists_P.push_back(hist.h_P_neg);
        if (hist.h_P_pos) hists_P.push_back(hist.h_P_pos);
        if (hist.h_R_neg) hists_R.push_back(hist.h_R_neg);
        if (hist.h_R_pos) hists_R.push_back(hist.h_R_pos);
        if (hist.h_Z_neg) hists_Z.push_back(hist.h_Z_neg);
        if (hist.h_Z_pos) hists_Z.push_back(hist.h_Z_pos);
    }
}
