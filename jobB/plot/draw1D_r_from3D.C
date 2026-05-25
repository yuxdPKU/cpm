#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);
TH1* SubtractHistograms(const TH1* h1, const TH1* h2, const char* name = "diff");
void draw1Dmap(TString filename, TString tag, double selectZ, TH1*& h_R_pos, TH1*& h_R_neg, TH1*& h_P_pos, TH1*& h_P_neg, TH1*& h_Z_pos, TH1*& h_Z_neg, TH1*& h_N_pos, TH1*& h_N_neg, int color=1, bool convert_RP_2_P=false);
void draw1Dmap_2D(TString filename, TString tag, TH1*& h_R_pos, TH1*& h_R_neg, TH1*& h_P_pos, TH1*& h_P_neg, TH1*& h_Z_pos, TH1*& h_Z_neg, int color);

void plot1D_Zbin(TH3* h3, TH1* h1, float z, bool convert_RP_2_P=false)
{
  int xbin = h3->GetXaxis()->FindBin((4.55503+4.85652)/2.);
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
  int xbin = h2->GetXaxis()->FindBin(4.55073);
  for (int i = 1; i <= h2->GetNbinsY(); i++)
  {
      int xminbin = h2->GetXaxis()->FindBin(4.55073);
      int xmaxbin = h2->GetXaxis()->FindBin(4.84711);
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

void draw1D_r_from3D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

const int nrun = 1;
//int runs[nrun] = {79507,79509,79510,79511,79512,79513,79514,79515,79516};
int runs[nrun] = {79516};

//get map from cdb lamination study
std::string cdbfilename_lamination[nrun];
TH1 *hcdb_lamination_N_pos[nrun], *hcdb_lamination_R_pos[nrun], *hcdb_lamination_P_pos[nrun], *hcdb_lamination_Z_pos[nrun];
TH1 *hcdb_lamination_N_neg[nrun], *hcdb_lamination_R_neg[nrun], *hcdb_lamination_P_neg[nrun], *hcdb_lamination_Z_neg[nrun];

//get map from cdb moduleedge study
std::string cdbfilename_moduleedge[nrun];
TH1 *hcdb_moduleedge_N_pos[nrun], *hcdb_moduleedge_R_pos[nrun], *hcdb_moduleedge_P_pos[nrun], *hcdb_moduleedge_Z_pos[nrun];
TH1 *hcdb_moduleedge_N_neg[nrun], *hcdb_moduleedge_R_neg[nrun], *hcdb_moduleedge_P_neg[nrun], *hcdb_moduleedge_Z_neg[nrun];

for (int i = 0; i < nrun; i++)
{
  auto rc = recoConsts::instance();
  rc->set_uint64Flag("TIMESTAMP", runs[i]);
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");

  cdbfilename_lamination[i] = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
  cout<<"Lamination CDB file path: "<<cdbfilename_lamination[i].c_str()<<endl;

  cdbfilename_moduleedge[i] = CDBInterface::instance()->getUrl("TPC_Module_Edge");
  cout<<"Module-edge CDB file path: "<<cdbfilename_moduleedge[i].c_str()<<endl;

  draw1Dmap_2D(cdbfilename_lamination[i].c_str(), Form("%d_lamination_cdb",runs[i]), hcdb_lamination_R_pos[i], hcdb_lamination_R_neg[i], hcdb_lamination_P_pos[i], hcdb_lamination_P_neg[i], hcdb_lamination_Z_pos[i], hcdb_lamination_Z_neg[i], 1);
}

std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
for (const auto& selectZ : selectZs)
{

TH1 *h_N_pos_method1[nrun], *h_R_pos_method1[nrun], *h_P_pos_method1[nrun], *h_Z_pos_method1[nrun];
TH1 *h_N_neg_method1[nrun], *h_R_neg_method1[nrun], *h_P_neg_method1[nrun], *h_Z_neg_method1[nrun];

TH1 *h_N_pos_method2[nrun], *h_R_pos_method2[nrun], *h_P_pos_method2[nrun], *h_Z_pos_method2[nrun];
TH1 *h_N_neg_method2[nrun], *h_R_neg_method2[nrun], *h_P_neg_method2[nrun], *h_Z_neg_method2[nrun];

for (int i = 0; i < nrun; i++)
{
  draw1Dmap(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Si_TPOT_fit/run3pp_newSiFieldonAlignment_newTPOTzfAlignment/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]), Form("%d_MI",runs[i]), selectZ, h_R_pos_method1[i], h_R_neg_method1[i], h_P_pos_method1[i], h_P_neg_method1[i], h_Z_pos_method1[i], h_Z_neg_method1[i], h_N_pos_method1[i], h_N_neg_method1[i], 2);

  draw1Dmap(Form("/sphenix/u/xyu3/workarea/cpm/jobB/output/run%d/run%d_B3_average_correction_histograms.root",runs[i],runs[i]), Form("%d_cpm",runs[i]), selectZ, h_R_pos_method2[i], h_R_neg_method2[i], h_P_pos_method2[i], h_P_neg_method2[i], h_Z_pos_method2[i], h_Z_neg_method2[i], h_N_pos_method2[i], h_N_neg_method2[i], 4);
}

TH1 *h_N_pos_method3, *h_R_pos_method3, *h_P_pos_method3, *h_Z_pos_method3;
TH1 *h_N_neg_method3, *h_R_neg_method3, *h_P_neg_method3, *h_Z_neg_method3;

draw1Dmap(Form("/sphenix/u/xyu3/workarea/cpm/jobB/output/sim/sim_B3_average_correction_histograms.root"), Form("sim_cpm"), selectZ, h_R_pos_method3, h_R_neg_method3, h_P_pos_method3, h_P_neg_method3, h_Z_pos_method3, h_Z_neg_method3, h_N_pos_method3, h_N_neg_method3, 1);

std::vector<TH1*> hists_P; hists_P.clear();
std::vector<TH1*> hists_R; hists_R.clear();
std::vector<TH1*> hists_Z; hists_Z.clear();
for (int i=0; i<nrun; i++)
{
	hists_P.push_back(h_P_neg_method1[i]);
	hists_P.push_back(h_P_pos_method1[i]);
	hists_P.push_back(h_P_neg_method2[i]);
	hists_P.push_back(h_P_pos_method2[i]);
	hists_P.push_back(hcdb_lamination_P_pos[i]);
	hists_P.push_back(hcdb_lamination_P_neg[i]);
}
for (int i=0; i<nrun; i++)
{
	hists_R.push_back(h_R_neg_method1[i]);
	hists_R.push_back(h_R_pos_method1[i]);
	hists_R.push_back(h_R_neg_method2[i]);
	hists_R.push_back(h_R_pos_method2[i]);
	hists_R.push_back(hcdb_lamination_R_pos[i]);
	hists_R.push_back(hcdb_lamination_R_neg[i]);
}
for (int i=0; i<nrun; i++)
{
	hists_Z.push_back(h_Z_neg_method1[i]);
	hists_Z.push_back(h_Z_pos_method1[i]);
	hists_Z.push_back(h_Z_neg_method2[i]);
	hists_Z.push_back(h_Z_pos_method2[i]);
	hists_Z.push_back(hcdb_lamination_Z_pos[i]);
	hists_Z.push_back(hcdb_lamination_Z_neg[i]);
}
std::pair<double,double> yrange_P = SetCommonYRange(hists_P);
std::pair<double,double> yrange_R = SetCommonYRange(hists_R);
std::pair<double,double> yrange_Z = SetCommonYRange(hists_Z);

TFile* ofile = new TFile(Form("Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos_method1[i]->Write();
  h_R_pos_method2[i]->Write();
  h_R_neg_method1[i]->Write();
  h_R_neg_method2[i]->Write();
  hcdb_lamination_R_neg[i]->Write();
}
ofile->Write();

TFile* ofile2 = new TFile(Form("Rootfiles/hist_dphi_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile2->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos_method1[i]->Write();
  h_P_pos_method2[i]->Write();
  h_P_neg_method1[i]->Write();
  h_P_neg_method2[i]->Write();
}
ofile2->Write();

TFile* ofile3 = new TFile(Form("Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile3->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos_method1[i]->Write();
  h_Z_pos_method2[i]->Write();
  h_Z_neg_method1[i]->Write();
  h_Z_neg_method2[i]->Write();
}
ofile3->Write();

for (int i=0; i<nrun; i++)
{
TCanvas* can_MI_CPM_LAM = new TCanvas("can_MI_CPM_LAM","",2400,1200);
can_MI_CPM_LAM->Divide(3,2);
can_MI_CPM_LAM->cd(1);
gPad->SetLogy(0);
h_P_pos_method1[i]->Draw("hist,e,same");
h_P_pos_method2[i]->Draw("hist,same");
hcdb_lamination_P_pos[i]->Draw("hist,same");
can_MI_CPM_LAM->cd(2);
gPad->SetLogy(0);
h_R_pos_method1[i]->Draw("hist,e,same");
h_R_pos_method2[i]->Draw("hist,same");
hcdb_lamination_R_pos[i]->Draw("hist,same");
can_MI_CPM_LAM->cd(3);
gPad->SetLogy(0);
h_Z_pos_method1[i]->Draw("hist,e,same");
h_Z_pos_method2[i]->Draw("hist,same");
hcdb_lamination_Z_pos[i]->Draw("hist,same");
can_MI_CPM_LAM->cd(4);
gPad->SetLogy(0);
h_P_neg_method1[i]->Draw("hist,e,same");
h_P_neg_method2[i]->Draw("hist,same");
hcdb_lamination_P_neg[i]->Draw("hist,same");
can_MI_CPM_LAM->cd(5);
gPad->SetLogy(0);
h_R_neg_method1[i]->Draw("hist,e,same");
h_R_neg_method2[i]->Draw("hist,same");
hcdb_lamination_R_neg[i]->Draw("hist,same");
can_MI_CPM_LAM->cd(6);
gPad->SetLogy(0);
h_Z_neg_method1[i]->Draw("hist,e,same");
h_Z_neg_method2[i]->Draw("hist,same");
hcdb_lamination_Z_neg[i]->Draw("hist,same");

gPad->RedrawAxis();

can_MI_CPM_LAM->Update();
can_MI_CPM_LAM->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_%d_MI_CPM_LAM.pdf",(int)selectZ,runs[i]));

TLegend *legend_MI_CPM_LAM = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_MI_CPM_LAM->SetHeader(Form("Run %d",runs[i]));
legend_MI_CPM_LAM->AddEntry(h_P_pos_method1[i], Form("Matrix Inversion"), "l");
legend_MI_CPM_LAM->AddEntry(h_P_pos_method2[i], Form("CPM"), "l");
legend_MI_CPM_LAM->AddEntry(hcdb_lamination_P_pos[i], Form("Lamination"), "l");
legend_MI_CPM_LAM->Draw();
TCanvas* can_MI_CPM_LAM_leg = new TCanvas("can_MI_CPM_LAM_leg","",2500,1000);
legend_MI_CPM_LAM->Draw();
can_MI_CPM_LAM_leg->SaveAs(Form("figure/resid_vsR_from3D_leg_%d_MI_CPM_LAM.pdf",runs[i]));

delete can_MI_CPM_LAM;
delete can_MI_CPM_LAM_leg;
}

TCanvas* can_sim = new TCanvas("can_sim","",2400,1200);
can_sim->Divide(4,2);
can_sim->cd(1);
gPad->SetLogy(0);
h_N_pos_method3->Draw("hist,same");
can_sim->cd(2);
gPad->SetLogy(0);
h_P_pos_method3->SetMaximum(0.03);
h_P_pos_method3->SetMinimum(-0.01);
h_P_pos_method3->Draw("hist,same");
can_sim->cd(3);
gPad->SetLogy(0);
h_R_pos_method3->SetMaximum(0.5);
h_R_pos_method3->SetMinimum(-0.5);
h_R_pos_method3->Draw("hist,same");
can_sim->cd(4);
gPad->SetLogy(0);
h_Z_pos_method3->SetMaximum(0.1);
h_Z_pos_method3->SetMinimum(-0.1);
h_Z_pos_method3->Draw("hist,same");
can_sim->cd(5);
gPad->SetLogy(0);
h_N_neg_method3->Draw("hist,same");
can_sim->cd(6);
gPad->SetLogy(0);
h_P_neg_method3->SetMaximum(0.03);
h_P_neg_method3->SetMinimum(-0.01);
h_P_neg_method3->Draw("hist,same");
can_sim->cd(7);
gPad->SetLogy(0);
h_R_neg_method3->SetMaximum(0.5);
h_R_neg_method3->SetMinimum(-0.5);
h_R_neg_method3->Draw("hist,same");
can_sim->cd(8);
gPad->SetLogy(0);
h_Z_neg_method3->SetMaximum(0.1);
h_Z_neg_method3->SetMinimum(-0.1);
h_Z_neg_method3->Draw("hist,same");

gPad->RedrawAxis();

can_sim->Update();
can_sim->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_sim.pdf",(int)selectZ));

delete can_sim;

}

}

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms)
{
    Double_t yMin = TMath::Infinity();
    Double_t yMax = -TMath::Infinity();

    for (TH1* h : histograms) {
        if (!h) continue;
	//recover default max and min
	h->SetMinimum();
	h->SetMaximum();
        yMin = TMath::Min(yMin, h->GetMinimum());
        yMax = TMath::Max(yMax, h->GetMaximum());
    }

    if (yMin == TMath::Infinity()) yMin = 0;
    if (yMax == -TMath::Infinity()) yMax = 1;

    if (yMin > 0) { yMin *= 0.8; }
    else if (yMin < 0) { yMin *= 1.2; }

    if (yMax > 0) { yMax *= 1.2; }
    else if (yMax < 0) { yMax *= 0.8; }

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

void draw1Dmap(TString filename, TString tag, double selectZ,
  TH1*& h_R_pos, TH1*& h_R_neg,
  TH1*& h_P_pos, TH1*& h_P_neg,
  TH1*& h_Z_pos, TH1*& h_Z_neg,
  TH1*& h_N_pos, TH1*& h_N_neg,
  int color=1,
  bool convert_RP_2_P=false)
{
  TFile* file_3D_map = new TFile(filename,"");
  if (!file_3D_map || file_3D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  TH3* h_R_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  TH3* h_P_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  TH3* h_Z_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  TH3* h_N_prz_pos = (TH3*) file_3D_map->Get("hentries_posz");
  TH3* h_R_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  TH3* h_P_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  TH3* h_Z_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");
  TH3* h_N_prz_neg = (TH3*) file_3D_map->Get("hentries_negz");

  if (!h_R_prz_pos || !h_P_prz_pos || !h_Z_prz_pos ||
      !h_R_prz_neg || !h_P_prz_neg || !h_Z_prz_neg) {
      std::cerr << "Error: missing 3D histograms in file" << std::endl;
      delete file_3D_map;
      return;
  }

  h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos->GetYaxis()->GetNbins(),h_N_prz_pos->GetYaxis()->GetXmin(),h_N_prz_pos->GetYaxis()->GetXmax());
  h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos->GetYaxis()->GetNbins(),h_R_prz_pos->GetYaxis()->GetXmin(),h_R_prz_pos->GetYaxis()->GetXmax());
  h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_pos->GetYaxis()->GetNbins(),h_P_prz_pos->GetYaxis()->GetXmin(),h_P_prz_pos->GetYaxis()->GetXmax());
  h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos->GetYaxis()->GetNbins(),h_Z_prz_pos->GetYaxis()->GetXmin(),h_Z_prz_pos->GetYaxis()->GetXmax());
  h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),h_N_prz_neg->GetYaxis()->GetNbins(),h_N_prz_neg->GetYaxis()->GetXmin(),h_N_prz_neg->GetYaxis()->GetXmax());
  h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg->GetYaxis()->GetNbins(),h_R_prz_neg->GetYaxis()->GetXmin(),h_R_prz_neg->GetYaxis()->GetXmax());
  h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ Z=-%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_neg->GetYaxis()->GetNbins(),h_P_prz_neg->GetYaxis()->GetXmin(),h_P_prz_neg->GetYaxis()->GetXmax());
  h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg->GetYaxis()->GetNbins(),h_Z_prz_neg->GetYaxis()->GetXmin(),h_Z_prz_neg->GetYaxis()->GetXmax());

  // do not save in the file_3D_map
  // directly saved in the stack memory
  h_N_pos->SetDirectory(0);
  h_R_pos->SetDirectory(0);
  h_P_pos->SetDirectory(0);
  h_Z_pos->SetDirectory(0);
  h_N_neg->SetDirectory(0);
  h_R_neg->SetDirectory(0);
  h_P_neg->SetDirectory(0);
  h_Z_neg->SetDirectory(0);

  plot1D_Zbin(h_N_prz_pos,h_N_pos,selectZ, false);
  plot1D_Zbin(h_R_prz_pos,h_R_pos,selectZ, false);
  plot1D_Zbin(h_P_prz_pos,h_P_pos,selectZ, convert_RP_2_P);
  plot1D_Zbin(h_Z_prz_pos,h_Z_pos,selectZ, false);
  plot1D_Zbin(h_N_prz_neg,h_N_neg,-selectZ, false);
  plot1D_Zbin(h_R_prz_neg,h_R_neg,-selectZ, false);
  plot1D_Zbin(h_P_prz_neg,h_P_neg,-selectZ, convert_RP_2_P);
  plot1D_Zbin(h_Z_prz_neg,h_Z_neg,-selectZ, false);

  h_N_pos->SetLineColor(color); h_N_pos->SetLineWidth(1); h_N_pos->SetFillColor(0); h_N_pos->SetMarkerColor(color);
  h_R_pos->SetLineColor(color); h_R_pos->SetLineWidth(1); h_R_pos->SetFillColor(0); h_R_pos->SetMarkerColor(color);
  h_P_pos->SetLineColor(color); h_P_pos->SetLineWidth(1); h_P_pos->SetFillColor(0); h_P_pos->SetMarkerColor(color);
  h_Z_pos->SetLineColor(color); h_Z_pos->SetLineWidth(1); h_Z_pos->SetFillColor(0); h_Z_pos->SetMarkerColor(color);
  h_N_neg->SetLineColor(color); h_N_neg->SetLineWidth(1); h_N_neg->SetFillColor(0); h_R_neg->SetMarkerColor(color);
  h_R_neg->SetLineColor(color); h_R_neg->SetLineWidth(1); h_R_neg->SetFillColor(0); h_R_neg->SetMarkerColor(color);
  h_P_neg->SetLineColor(color); h_P_neg->SetLineWidth(1); h_P_neg->SetFillColor(0); h_P_neg->SetMarkerColor(color);
  h_Z_neg->SetLineColor(color); h_Z_neg->SetLineWidth(1); h_Z_neg->SetFillColor(0); h_Z_neg->SetMarkerColor(color);

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
