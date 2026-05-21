#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);
TH1* SubtractHistograms(const TH1* h1, const TH1* h2, const char* name = "diff");
void draw1Dmap(TString filename, TString tag, double selectR, TH1*& h_R_pos, TH1*& h_R_neg, TH1*& h_P_pos, TH1*& h_P_neg, TH1*& h_Z_pos, TH1*& h_Z_neg, int color=1, bool convert_RP_2_P=false);
void draw1Dmap_2D(TString filename, TString tag, double selectR, TLine*& l_R_pos, TLine*& l_R_neg, TLine*& l_P_pos, TLine*& l_P_neg, TLine*& l_Z_pos, TLine*& l_Z_neg, int color=1);

void plot1D_Zbin(TH3* h3, TH1* h1, float y, bool convert_RP_2_P=false)
{
  int xbin = h3->GetXaxis()->FindBin((4.55503+4.85652)/2.);
  int ybin = h3->GetYaxis()->FindBin(y);
  for (int i = 1; i <= h3->GetNbinsZ(); i++)
  {
      double value = h3->GetBinContent(xbin, ybin, i);
      double error = h3->GetBinError(xbin, ybin, i);
      if (convert_RP_2_P)
      {
        value /= y;
        error /= y;
      }
      h1->SetBinContent(i, value);
      h1->SetBinError(i, error);
  }
}

void plot1D_line(TH2* h2, TLine*& line, float y, bool is_north_or_south)
{
  int ybin = h2->GetYaxis()->FindBin(y);
  int xminbin = h2->GetXaxis()->FindBin(4.55073);
  int xmaxbin = h2->GetXaxis()->FindBin(4.84711);
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

void draw1D_z_from3D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::pair<double, double> phirange_NCO = {0.609928,0.917104};
std::pair<double, double> phirange_NCI = {0.0276558,0.566484};
std::pair<double, double> phirange_SCI = {-0.568656,-0.0325273};
std::pair<double, double> phirange_SCO = {-0.788212,-0.613404};

const int nrun = 1;
//int runs[nrun] = {79507,79509,79510,79511,79512,79513,79514,79515,79516};
int runs[nrun] = {79516};

std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
for (const auto& selectR : selectRs)
{

//get map from cdb lamination study
std::string cdbfilename_lamination[nrun];
TLine *lcdb_lamination_N_pos[nrun], *lcdb_lamination_R_pos[nrun], *lcdb_lamination_P_pos[nrun], *lcdb_lamination_Z_pos[nrun];
TLine *lcdb_lamination_N_neg[nrun], *lcdb_lamination_R_neg[nrun], *lcdb_lamination_P_neg[nrun], *lcdb_lamination_Z_neg[nrun];

TH1 *h_N_pos_method1[nrun], *h_R_pos_method1[nrun], *h_P_pos_method1[nrun], *h_Z_pos_method1[nrun];
TH1 *h_N_neg_method1[nrun], *h_R_neg_method1[nrun], *h_P_neg_method1[nrun], *h_Z_neg_method1[nrun];

TH1 *h_N_pos_method2[nrun], *h_R_pos_method2[nrun], *h_P_pos_method2[nrun], *h_Z_pos_method2[nrun];
TH1 *h_N_neg_method2[nrun], *h_R_neg_method2[nrun], *h_P_neg_method2[nrun], *h_Z_neg_method2[nrun];

TLine *l_nco_o_R[nrun], *l_nco_i_R[nrun], *l_nci_o_R[nrun], *l_nci_i_R[nrun];
TLine *l_nco_o_P[nrun], *l_nco_i_P[nrun], *l_nci_o_P[nrun], *l_nci_i_P[nrun];
TLine *l_nco_o_Z[nrun], *l_nco_i_Z[nrun], *l_nci_o_Z[nrun], *l_nci_i_Z[nrun];
TLine *l_nco_o_RP[nrun], *l_nco_i_RP[nrun], *l_nci_o_RP[nrun], *l_nci_i_RP[nrun];
TLine *l_sci_i_R[nrun], *l_sci_o_R[nrun], *l_sco_i_R[nrun], *l_sco_o_R[nrun];
TLine *l_sci_i_P[nrun], *l_sci_o_P[nrun], *l_sco_i_P[nrun], *l_sco_o_P[nrun];
TLine *l_sci_i_Z[nrun], *l_sci_o_Z[nrun], *l_sco_i_Z[nrun], *l_sco_o_Z[nrun];
TLine *l_sci_i_RP[nrun], *l_sci_o_RP[nrun], *l_sco_i_RP[nrun], *l_sco_o_RP[nrun];

for (int i = 0; i < nrun; i++)
{

  auto rc = recoConsts::instance();
  rc->set_uint64Flag("TIMESTAMP", runs[i]);
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");

  cdbfilename_lamination[i] = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
  cout<<"Lamination CDB file path: "<<cdbfilename_lamination[i].c_str()<<endl;

  draw1Dmap_2D(cdbfilename_lamination[i].c_str(), Form("%d_lamination_cdb",runs[i]), selectR, lcdb_lamination_R_pos[i], lcdb_lamination_R_neg[i], lcdb_lamination_P_pos[i], lcdb_lamination_P_neg[i], lcdb_lamination_Z_pos[i], lcdb_lamination_Z_neg[i], 1);

  draw1Dmap(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Si_TPOT_fit/run3pp_newSiFieldonAlignment_newTPOTzfAlignment/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]), Form("%d_MI",runs[i]), selectR, h_R_pos_method1[i], h_R_neg_method1[i], h_P_pos_method1[i], h_P_neg_method1[i], h_Z_pos_method1[i], h_Z_neg_method1[i], 2);

  draw1Dmap(Form("/sphenix/u/xyu3/workarea/cpm/jobB/output/run%d/run%d_B3_average_correction_histograms.root",runs[i],runs[i]), Form("%d_cpm",runs[i]), selectR, h_R_pos_method2[i], h_R_neg_method2[i], h_P_pos_method2[i], h_P_neg_method2[i], h_Z_pos_method2[i], h_Z_neg_method2[i], 4);

}

TH1 *h_N_pos_method3, *h_R_pos_method3, *h_P_pos_method3, *h_Z_pos_method3;
TH1 *h_N_neg_method3, *h_R_neg_method3, *h_P_neg_method3, *h_Z_neg_method3;
draw1Dmap(Form("/sphenix/u/xyu3/workarea/cpm/jobB/output/sim/sim_B3_average_correction_histograms.root"), Form("sim_cpm"), selectR, h_R_pos_method3, h_R_neg_method3, h_P_pos_method3, h_P_neg_method3, h_Z_pos_method3, h_Z_neg_method3, 1);

std::vector<TH1*> hists_P; hists_P.clear();
std::vector<TH1*> hists_R; hists_R.clear();
std::vector<TH1*> hists_Z; hists_Z.clear();
for (int i=0; i<nrun; i++)
{
	hists_P.push_back(h_P_neg_method1[i]);
	hists_P.push_back(h_P_pos_method1[i]);
	hists_P.push_back(h_P_neg_method2[i]);
	hists_P.push_back(h_P_pos_method2[i]);
}
for (int i=0; i<nrun; i++)
{
	hists_R.push_back(h_R_neg_method1[i]);
	hists_R.push_back(h_R_pos_method1[i]);
	hists_R.push_back(h_R_neg_method2[i]);
	hists_R.push_back(h_R_pos_method2[i]);
}
for (int i=0; i<nrun; i++)
{
	hists_Z.push_back(h_Z_neg_method1[i]);
	hists_Z.push_back(h_Z_pos_method1[i]);
	hists_Z.push_back(h_Z_neg_method2[i]);
	hists_Z.push_back(h_Z_pos_method2[i]);
}
std::pair<double,double> yrange_P = SetCommonYRange(hists_P);
std::pair<double,double> yrange_R = SetCommonYRange(hists_R);
std::pair<double,double> yrange_Z = SetCommonYRange(hists_Z);

//special case

//for (int i=0; i<(hists_Z.size()); i++)
//{
//	hists_Z[i]->SetMaximum(0.5);
//	hists_Z[i]->SetMinimum(-0.5);
//}

//special case

for (int i=0; i<nrun; i++)
{
  l_nco_o_R[i] = new TLine(selectR*tan(phirange_NCO.second),yrange_R.first,selectR*tan(phirange_NCO.second),yrange_R.second);
  l_nco_o_P[i] = new TLine(selectR*tan(phirange_NCO.second),yrange_P.first,selectR*tan(phirange_NCO.second),yrange_P.second);
  l_nco_o_Z[i] = new TLine(selectR*tan(phirange_NCO.second),yrange_Z.first,selectR*tan(phirange_NCO.second),yrange_Z.second);

  l_nco_i_R[i] = new TLine(selectR*tan(phirange_NCO.first),yrange_R.first,selectR*tan(phirange_NCO.first),yrange_R.second);
  l_nco_i_P[i] = new TLine(selectR*tan(phirange_NCO.first),yrange_P.first,selectR*tan(phirange_NCO.first),yrange_P.second);
  l_nco_i_Z[i] = new TLine(selectR*tan(phirange_NCO.first),yrange_Z.first,selectR*tan(phirange_NCO.first),yrange_Z.second);

  l_nci_o_R[i] = new TLine(selectR*tan(phirange_NCI.second),yrange_R.first,selectR*tan(phirange_NCI.second),yrange_R.second);
  l_nci_o_P[i] = new TLine(selectR*tan(phirange_NCI.second),yrange_P.first,selectR*tan(phirange_NCI.second),yrange_P.second);
  l_nci_o_Z[i] = new TLine(selectR*tan(phirange_NCI.second),yrange_Z.first,selectR*tan(phirange_NCI.second),yrange_Z.second);

  l_nci_i_R[i] = new TLine(selectR*tan(phirange_NCI.first),yrange_R.first,selectR*tan(phirange_NCI.first),yrange_R.second);
  l_nci_i_P[i] = new TLine(selectR*tan(phirange_NCI.first),yrange_P.first,selectR*tan(phirange_NCI.first),yrange_P.second);
  l_nci_i_Z[i] = new TLine(selectR*tan(phirange_NCI.first),yrange_Z.first,selectR*tan(phirange_NCI.first),yrange_Z.second);

  l_sci_i_R[i] = new TLine(selectR*tan(phirange_SCI.second),yrange_R.first,selectR*tan(phirange_SCI.second),yrange_R.second);
  l_sci_i_P[i] = new TLine(selectR*tan(phirange_SCI.second),yrange_P.first,selectR*tan(phirange_SCI.second),yrange_P.second);
  l_sci_i_Z[i] = new TLine(selectR*tan(phirange_SCI.second),yrange_Z.first,selectR*tan(phirange_SCI.second),yrange_Z.second);

  l_sci_o_R[i] = new TLine(selectR*tan(phirange_SCI.first),yrange_R.first,selectR*tan(phirange_SCI.first),yrange_R.second);
  l_sci_o_P[i] = new TLine(selectR*tan(phirange_SCI.first),yrange_P.first,selectR*tan(phirange_SCI.first),yrange_P.second);
  l_sci_o_Z[i] = new TLine(selectR*tan(phirange_SCI.first),yrange_Z.first,selectR*tan(phirange_SCI.first),yrange_Z.second);

  l_sco_i_R[i] = new TLine(selectR*tan(phirange_SCO.second),yrange_R.first,selectR*tan(phirange_SCO.second),yrange_R.second);
  l_sco_i_P[i] = new TLine(selectR*tan(phirange_SCO.second),yrange_P.first,selectR*tan(phirange_SCO.second),yrange_P.second);
  l_sco_i_Z[i] = new TLine(selectR*tan(phirange_SCO.second),yrange_Z.first,selectR*tan(phirange_SCO.second),yrange_Z.second);

  l_sco_o_R[i] = new TLine(selectR*tan(phirange_SCO.first),yrange_R.first,selectR*tan(phirange_SCO.first),yrange_R.second);
  l_sco_o_P[i] = new TLine(selectR*tan(phirange_SCO.first),yrange_P.first,selectR*tan(phirange_SCO.first),yrange_P.second);
  l_sco_o_Z[i] = new TLine(selectR*tan(phirange_SCO.first),yrange_Z.first,selectR*tan(phirange_SCO.first),yrange_Z.second);
}

TFile* ofile = new TFile(Form("Rootfiles/hist_dr_1Dz_R%d.root",(int)selectR),"recreate");
ofile->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos_method1[i]->Write();
  h_R_pos_method2[i]->Write();
  h_R_neg_method1[i]->Write();
  h_R_neg_method2[i]->Write();
}
ofile->Write();

TFile* ofile2 = new TFile(Form("Rootfiles/hist_dphi_1Dz_R%d.root",(int)selectR),"recreate");
ofile2->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos_method1[i]->Write();
  h_P_pos_method2[i]->Write();
  h_P_neg_method1[i]->Write();
  h_P_neg_method2[i]->Write();
}
ofile2->Write();

TFile* ofile3 = new TFile(Form("Rootfiles/hist_dz_1Dz_R%d.root",(int)selectR),"recreate");
ofile3->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos_method1[i]->Write();
  h_Z_pos_method2[i]->Write();
  h_Z_neg_method1[i]->Write();
  h_Z_neg_method2[i]->Write();
}
ofile3->Write();

for (int i = 0; i < nrun; i++)
{
TCanvas* can_MI_CPM_LAM = new TCanvas("can_MI_CPM_LAM","",2400,1200);
can_MI_CPM_LAM->Divide(3,2);
can_MI_CPM_LAM->cd(1);
gPad->SetLogy(0);
h_P_pos_method1[i]->Draw("hist,e,same");
h_P_pos_method2[i]->Draw("hist,same");
l_nco_i_P[i]->Draw();
l_nco_o_P[i]->Draw();
l_nci_i_P[i]->Draw();
l_nci_o_P[i]->Draw();
can_MI_CPM_LAM->cd(2);
gPad->SetLogy(0);
h_R_pos_method1[i]->Draw("hist,e,same");
h_R_pos_method2[i]->Draw("hist,same");
l_nco_i_R[i]->Draw();
l_nco_o_R[i]->Draw();
l_nci_i_R[i]->Draw();
l_nci_o_R[i]->Draw();
can_MI_CPM_LAM->cd(3);
gPad->SetLogy(0);
h_Z_pos_method1[i]->Draw("hist,e,same");
h_Z_pos_method2[i]->Draw("hist,same");
l_nco_i_Z[i]->Draw();
l_nco_o_Z[i]->Draw();
l_nci_i_Z[i]->Draw();
l_nci_o_Z[i]->Draw();
can_MI_CPM_LAM->cd(4);
gPad->SetLogy(0);
h_P_neg_method1[i]->Draw("hist,e,same");
h_P_neg_method2[i]->Draw("hist,same");
l_sco_i_P[i]->Draw();
l_sco_o_P[i]->Draw();
l_sci_i_P[i]->Draw();
l_sci_o_P[i]->Draw();
can_MI_CPM_LAM->cd(5);
gPad->SetLogy(0);
h_R_neg_method1[i]->Draw("hist,e,same");
h_R_neg_method2[i]->Draw("hist,same");
l_sco_i_R[i]->Draw();
l_sco_o_R[i]->Draw();
l_sci_i_R[i]->Draw();
l_sci_o_R[i]->Draw();
can_MI_CPM_LAM->cd(6);
gPad->SetLogy(0);
h_Z_neg_method1[i]->Draw("hist,e,same");
h_Z_neg_method2[i]->Draw("hist,same");
l_sco_i_Z[i]->Draw();
l_sco_o_Z[i]->Draw();
l_sci_i_Z[i]->Draw();
l_sci_o_Z[i]->Draw();

gPad->RedrawAxis();

can_MI_CPM_LAM->Update();
can_MI_CPM_LAM->SaveAs(Form("figure/resid_vsZ_from3D_atR%d_%d.pdf",(int)selectR,runs[i]));

TLegend *legend_MI_CPM_LAM = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_MI_CPM_LAM->SetHeader(Form("Run %d",runs[i]));
legend_MI_CPM_LAM->AddEntry(h_P_pos_method1[i], Form("Matrix Inversion"), "l");
legend_MI_CPM_LAM->AddEntry(h_P_pos_method2[i], Form("CPM"), "l");
legend_MI_CPM_LAM->Draw();
TCanvas* can_MI_CPM_LAM_leg = new TCanvas("can_MI_CPM_LAM_leg","",2500,1000);
legend_MI_CPM_LAM->Draw();
can_MI_CPM_LAM_leg->SaveAs(Form("figure/resid_vsZ_from3D_leg_%d_MI_CPM_LAM.pdf",runs[i]));

delete can_MI_CPM_LAM;
delete can_MI_CPM_LAM_leg;

}

TCanvas* can_sim = new TCanvas("can_sim","",2400,1200);
can_sim->Divide(3,2);
can_sim->cd(1);
gPad->SetLogy(0);
h_P_pos_method3->Draw("hist,same");
can_sim->cd(2);
gPad->SetLogy(0);
h_R_pos_method3->Draw("hist,same");
can_sim->cd(3);
gPad->SetLogy(0);
h_Z_pos_method3->Draw("hist,same");
can_sim->cd(4);
gPad->SetLogy(0);
h_P_neg_method3->Draw("hist,same");
can_sim->cd(5);
gPad->SetLogy(0);
h_R_neg_method3->Draw("hist,same");
can_sim->cd(6);
gPad->SetLogy(0);
h_Z_neg_method3->Draw("hist,same");

gPad->RedrawAxis();

can_sim->Update();
can_sim->SaveAs(Form("figure/resid_vsZ_from3D_atR%d_sim.pdf",(int)selectR));

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

void draw1Dmap(TString filename, TString tag, double selectR,
  TH1*& h_R_pos, TH1*& h_R_neg,
  TH1*& h_P_pos, TH1*& h_P_neg,
  TH1*& h_Z_pos, TH1*& h_Z_neg,
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
  //TH3* h_N_prz_pos = (TH3*) file_3D_map->Get("hentries_posz");
  TH3* h_R_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  TH3* h_P_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  TH3* h_Z_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");
  //TH3* h_N_prz_neg = (TH3*) file_3D_map->Get("hentries_negz");

  if (!h_R_prz_pos || !h_P_prz_pos || !h_Z_prz_pos ||
      !h_R_prz_neg || !h_P_prz_neg || !h_Z_prz_neg) {
      std::cerr << "Error: missing 3D histograms in file" << std::endl;
      delete file_3D_map;
      return;
  }

  //h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_pos->GetZaxis()->GetNbins(),h_N_prz_pos->GetZaxis()->GetXmin(),h_N_prz_pos->GetZaxis()->GetXmax());
  h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_pos->GetZaxis()->GetNbins(),h_R_prz_pos->GetZaxis()->GetXmin(),h_R_prz_pos->GetZaxis()->GetXmax());
  h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_pos->GetZaxis()->GetNbins(),h_P_prz_pos->GetZaxis()->GetXmin(),h_P_prz_pos->GetZaxis()->GetXmax());
  h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_pos->GetZaxis()->GetNbins(),h_Z_prz_pos->GetZaxis()->GetXmin(),h_Z_prz_pos->GetZaxis()->GetXmax());
  //h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_neg->GetZaxis()->GetNbins(),h_N_prz_neg->GetZaxis()->GetXmin(),h_N_prz_neg->GetZaxis()->GetXmax());
  h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_neg->GetZaxis()->GetNbins(),h_R_prz_neg->GetZaxis()->GetXmin(),h_R_prz_neg->GetZaxis()->GetXmax());
  h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_neg->GetZaxis()->GetNbins(),h_P_prz_neg->GetZaxis()->GetXmin(),h_P_prz_neg->GetZaxis()->GetXmax());
  h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_neg->GetZaxis()->GetNbins(),h_Z_prz_neg->GetZaxis()->GetXmin(),h_Z_prz_neg->GetZaxis()->GetXmax());

  // do not save in the file_3D_map
  // directly saved in the stack memory
  //h_N_pos->SetDirectory(0);
  h_R_pos->SetDirectory(0);
  h_P_pos->SetDirectory(0);
  h_Z_pos->SetDirectory(0);
  //h_N_neg->SetDirectory(0);
  h_R_neg->SetDirectory(0);
  h_P_neg->SetDirectory(0);
  h_Z_neg->SetDirectory(0);

  //plot1D_Zbin(h_N_prz_pos,h_N_pos,selectR, false);
  plot1D_Zbin(h_R_prz_pos,h_R_pos,selectR, false);
  plot1D_Zbin(h_P_prz_pos,h_P_pos,selectR, convert_RP_2_P);
  plot1D_Zbin(h_Z_prz_pos,h_Z_pos,selectR, false);
  //plot1D_Zbin(h_N_prz_neg,h_N_neg,selectR, false);
  plot1D_Zbin(h_R_prz_neg,h_R_neg,selectR, false);
  plot1D_Zbin(h_P_prz_neg,h_P_neg,selectR, convert_RP_2_P);
  plot1D_Zbin(h_Z_prz_neg,h_Z_neg,selectR, false);

  //h_N_pos->SetLineColor(color); h_N_pos->SetLineWidth(1); h_N_pos->SetFillColor(0); h_N_pos->SetMarkerColor(color);
  h_R_pos->SetLineColor(color); h_R_pos->SetLineWidth(1); h_R_pos->SetFillColor(0); h_R_pos->SetMarkerColor(color);
  h_P_pos->SetLineColor(color); h_P_pos->SetLineWidth(1); h_P_pos->SetFillColor(0); h_P_pos->SetMarkerColor(color);
  h_Z_pos->SetLineColor(color); h_Z_pos->SetLineWidth(1); h_Z_pos->SetFillColor(0); h_Z_pos->SetMarkerColor(color);
  //h_N_neg->SetLineColor(color); h_N_neg->SetLineWidth(1); h_N_neg->SetFillColor(0); h_R_neg->SetMarkerColor(color);
  h_R_neg->SetLineColor(color); h_R_neg->SetLineWidth(1); h_R_neg->SetFillColor(0); h_R_neg->SetMarkerColor(color);
  h_P_neg->SetLineColor(color); h_P_neg->SetLineWidth(1); h_P_neg->SetFillColor(0); h_P_neg->SetMarkerColor(color);
  h_Z_neg->SetLineColor(color); h_Z_neg->SetLineWidth(1); h_Z_neg->SetFillColor(0); h_Z_neg->SetMarkerColor(color);

  delete file_3D_map;

  return;
}

void draw1Dmap_2D(TString filename, TString tag, double selectR,
  TLine*& l_R_pos, TLine*& l_R_neg,
  TLine*& l_P_pos, TLine*& l_P_neg,
  TLine*& l_Z_pos, TLine*& l_Z_neg,
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
