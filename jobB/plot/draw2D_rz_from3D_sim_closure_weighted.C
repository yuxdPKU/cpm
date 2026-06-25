#include <TGaxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPad.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

struct threeDHist {
    TH3* h_N_pos;
    TH3* h_P_pos;
    TH3* h_R_pos;
    TH3* h_Z_pos;
    TH3* h_N_neg;
    TH3* h_P_neg;
    TH3* h_R_neg;
    TH3* h_Z_neg;
};

struct twoDHist {
    TH2* h_N_pos;
    TH2* h_P_pos;
    TH2* h_R_pos;
    TH2* h_Z_pos;
    TH2* h_N_neg;
    TH2* h_P_neg;
    TH2* h_R_neg;
    TH2* h_Z_neg;
};

TH2* MakeRZHist(TH3* h3, TString name, TString title)
{
  return new TH2F(name, title,
                  h3->GetYaxis()->GetNbins(), h3->GetYaxis()->GetXmin(), h3->GetYaxis()->GetXmax(),
                  h3->GetZaxis()->GetNbins(), h3->GetZaxis()->GetXmin(), h3->GetZaxis()->GetXmax());
}

void plot2D_Xbin(TH3* h3, TH2* h2, float phi)
{
  int xbin = h3->GetXaxis()->FindBin(phi);
  for (int i = 1; i <= h3->GetNbinsY(); i++)
  {
    for (int j = 1; j <= h3->GetNbinsZ(); j++)
    {
      double value = h3->GetBinContent(xbin, i, j);
      double error = h3->GetBinError(xbin, i, j);
      if (std::isnan(value))
      {
        value = 0.0;
        error = 0.0;
      }
      h2->SetBinContent(i, j, value);
      h2->SetBinError(i, j, error);
    }
  }
}

int Mask2DWithN(TH2* h_N, TH2* h)
{
  int nMasked = 0;
  for (int i = 1; i <= h->GetNbinsX(); i++)
  {
    for (int j = 1; j <= h->GetNbinsY(); j++)
    {
      if (h_N->GetBinContent(i, j) <= 0)
      {
        h->SetBinContent(i, j, std::numeric_limits<double>::quiet_NaN());
        h->SetBinError(i, j, 0.0);
        nMasked++;
      }
    }
  }
  return nMasked;
}

void SetCommonZRange(const std::vector<TH2*>& histograms, bool positiveOnly=false)
{
  double zmin = 1e30;
  double zmax = -1e30;
  for (TH2* h : histograms)
  {
    if (!h) continue;
    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
      for (int j = 1; j <= h->GetNbinsY(); j++)
      {
        double value = h->GetBinContent(i, j);
        if (!std::isfinite(value)) continue;
        zmin = std::min(zmin, value);
        zmax = std::max(zmax, value);
      }
    }
  }

  if (zmin == 1e30) zmin = 0;
  if (zmax == -1e30) zmax = 1;

  if (positiveOnly)
  {
    zmin = 0;
    if (zmax <= 0) zmax = 1;
  }
  else
  {
    if (zmin < 0) zmin *= 1.1; else zmin = -0.1;
    if (zmax > 0) zmax *= 1.1; else zmax = 0.1;
  }

  for (TH2* h : histograms)
  {
    if (!h) continue;
    h->SetMinimum(zmin);
    h->SetMaximum(zmax);
  }
}

void DrawMaskedBinsWhite(TH2* h_N)
{
  for (int i = 1; i <= h_N->GetNbinsX(); i++)
  {
    double xlow = h_N->GetXaxis()->GetBinLowEdge(i);
    double xup = h_N->GetXaxis()->GetBinUpEdge(i);
    for (int j = 1; j <= h_N->GetNbinsY(); j++)
    {
      if (h_N->GetBinContent(i, j) > 0) continue;
      double ylow = h_N->GetYaxis()->GetBinLowEdge(j);
      double yup = h_N->GetYaxis()->GetBinUpEdge(j);
      TBox* box = new TBox(xlow, ylow, xup, yup);
      box->SetFillColor(kWhite);
      box->SetLineColor(kWhite);
      box->Draw("same");
    }
  }
}

void Draw2D(TH2* h, TH2* h_N=nullptr)
{
  gPad->SetRightMargin(0.15);
  h->Draw(h_N ? "colz0" : "colz");
  if (h_N) DrawMaskedBinsWhite(h_N);
  gPad->RedrawAxis();
}

bool Load3DHists(TFile* file_3D_map, threeDHist& hists3D)
{
  hists3D.h_N_pos = (TH3*) file_3D_map->Get("hentries_posz");
  hists3D.h_P_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  hists3D.h_R_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  hists3D.h_Z_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  hists3D.h_N_neg = (TH3*) file_3D_map->Get("hentries_negz");
  hists3D.h_P_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  hists3D.h_R_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  hists3D.h_Z_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");

  return hists3D.h_N_pos && hists3D.h_P_pos && hists3D.h_R_pos && hists3D.h_Z_pos &&
         hists3D.h_N_neg && hists3D.h_P_neg && hists3D.h_R_neg && hists3D.h_Z_neg;
}

void draw2D_rz_from3D_sim_closure_weighted(int run=29, TString mode="acts", TString phiType="central")
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);
gSystem->mkdir("figure", true);
gSystem->mkdir("Rootfiles", true);

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
  std::cerr << "Please use phiType = central, east, or west." << std::endl;
  return;
}

TString tag = Form("CPM_sim_reco_%s", mode.Data());
TString filename = Form("../output/sim_run%d_%s_weighted/sim_run%d_%s_weighted_maxdca0p2_B3_average_correction_histograms.root", run, mode.Data(), run, mode.Data());
TFile* file_3D_map = new TFile(filename, "");
if (!file_3D_map || file_3D_map->IsZombie())
{
  std::cerr << "Error opening file: " << filename << std::endl;
  return;
}

threeDHist hists3D;
if (!Load3DHists(file_3D_map, hists3D))
{
  std::cerr << "Error: missing N/P/R/Z 3D histograms in " << filename << std::endl;
  delete file_3D_map;
  return;
}

int xbin = hists3D.h_N_pos->GetXaxis()->FindBin(selectPhi);
double phiBinCenter = hists3D.h_N_pos->GetXaxis()->GetBinCenter(xbin);
std::cout << "Use phiType = " << phiType
          << ", phi = " << selectPhi
          << ", xbin = " << xbin
          << ", bin center = " << phiBinCenter
          << std::endl;

twoDHist hists2D;
hists2D.h_N_pos = MakeRZHist(hists3D.h_N_pos, Form("hentries_posz_%s",tag.Data()), Form("N posz @ #phi=%.3f;R (cm);Z (cm);N",selectPhi));
hists2D.h_P_pos = MakeRZHist(hists3D.h_P_pos, Form("hIntDistortionP_posz_%s",tag.Data()), Form("d#phi posz @ #phi=%.3f;R (cm);Z (cm);d#phi (rad)",selectPhi));
hists2D.h_R_pos = MakeRZHist(hists3D.h_R_pos, Form("hIntDistortionR_posz_%s",tag.Data()), Form("dR posz @ #phi=%.3f;R (cm);Z (cm);dR (cm)",selectPhi));
hists2D.h_Z_pos = MakeRZHist(hists3D.h_Z_pos, Form("hIntDistortionZ_posz_%s",tag.Data()), Form("dz posz @ #phi=%.3f;R (cm);Z (cm);dz (cm)",selectPhi));
hists2D.h_N_neg = MakeRZHist(hists3D.h_N_neg, Form("hentries_negz_%s",tag.Data()), Form("N negz @ #phi=%.3f;R (cm);Z (cm);N",selectPhi));
hists2D.h_P_neg = MakeRZHist(hists3D.h_P_neg, Form("hIntDistortionP_negz_%s",tag.Data()), Form("d#phi negz @ #phi=%.3f;R (cm);Z (cm);d#phi (rad)",selectPhi));
hists2D.h_R_neg = MakeRZHist(hists3D.h_R_neg, Form("hIntDistortionR_negz_%s",tag.Data()), Form("dR negz @ #phi=%.3f;R (cm);Z (cm);dR (cm)",selectPhi));
hists2D.h_Z_neg = MakeRZHist(hists3D.h_Z_neg, Form("hIntDistortionZ_negz_%s",tag.Data()), Form("dz negz @ #phi=%.3f;R (cm);Z (cm);dz (cm)",selectPhi));

hists2D.h_N_pos->SetDirectory(0);
hists2D.h_P_pos->SetDirectory(0);
hists2D.h_R_pos->SetDirectory(0);
hists2D.h_Z_pos->SetDirectory(0);
hists2D.h_N_neg->SetDirectory(0);
hists2D.h_P_neg->SetDirectory(0);
hists2D.h_R_neg->SetDirectory(0);
hists2D.h_Z_neg->SetDirectory(0);

plot2D_Xbin(hists3D.h_N_pos, hists2D.h_N_pos, selectPhi);
plot2D_Xbin(hists3D.h_P_pos, hists2D.h_P_pos, selectPhi);
plot2D_Xbin(hists3D.h_R_pos, hists2D.h_R_pos, selectPhi);
plot2D_Xbin(hists3D.h_Z_pos, hists2D.h_Z_pos, selectPhi);
plot2D_Xbin(hists3D.h_N_neg, hists2D.h_N_neg, selectPhi);
plot2D_Xbin(hists3D.h_P_neg, hists2D.h_P_neg, selectPhi);
plot2D_Xbin(hists3D.h_R_neg, hists2D.h_R_neg, selectPhi);
plot2D_Xbin(hists3D.h_Z_neg, hists2D.h_Z_neg, selectPhi);

int nMaskedPos = Mask2DWithN(hists2D.h_N_pos, hists2D.h_P_pos);
Mask2DWithN(hists2D.h_N_pos, hists2D.h_R_pos);
Mask2DWithN(hists2D.h_N_pos, hists2D.h_Z_pos);
int nMaskedNeg = Mask2DWithN(hists2D.h_N_neg, hists2D.h_P_neg);
Mask2DWithN(hists2D.h_N_neg, hists2D.h_R_neg);
Mask2DWithN(hists2D.h_N_neg, hists2D.h_Z_neg);
std::cout << "Mask N=0 bins in P/R/Z: posz = " << nMaskedPos
          << ", negz = " << nMaskedNeg
          << std::endl;

SetCommonZRange({hists2D.h_N_pos, hists2D.h_N_neg}, true);
SetCommonZRange({hists2D.h_P_pos, hists2D.h_P_neg});
SetCommonZRange({hists2D.h_R_pos, hists2D.h_R_neg});
SetCommonZRange({hists2D.h_Z_pos, hists2D.h_Z_neg});

TCanvas* can = new TCanvas("can","",6400,2600);
can->Divide(4,2);
can->cd(1); Draw2D(hists2D.h_N_pos);
can->cd(2); Draw2D(hists2D.h_P_pos, hists2D.h_N_pos);
can->cd(3); Draw2D(hists2D.h_R_pos, hists2D.h_N_pos);
can->cd(4); Draw2D(hists2D.h_Z_pos, hists2D.h_N_pos);
can->cd(5); Draw2D(hists2D.h_N_neg);
can->cd(6); Draw2D(hists2D.h_P_neg, hists2D.h_N_neg);
can->cd(7); Draw2D(hists2D.h_R_neg, hists2D.h_N_neg);
can->cd(8); Draw2D(hists2D.h_Z_neg, hists2D.h_N_neg);
gPad->RedrawAxis();
can->Update();
can->SaveAs(Form("figure/rz_NPRZ_from3D_%d_%s_weighted_%s.pdf", run, tag.Data(), phiType.Data()));
delete can;

TFile* outfile = new TFile(Form("Rootfiles/hist_NPRZ_2Drz_%d_%s_weighted_%s.root", run, tag.Data(), phiType.Data()), "recreate");
outfile->cd();
hists2D.h_N_pos->Write();
hists2D.h_P_pos->Write();
hists2D.h_R_pos->Write();
hists2D.h_Z_pos->Write();
hists2D.h_N_neg->Write();
hists2D.h_P_neg->Write();
hists2D.h_R_neg->Write();
hists2D.h_Z_neg->Write();
outfile->Write();
outfile->Close();

delete file_3D_map;
}
