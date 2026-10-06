#include "utilities.h"

void prepare_1D_plot(int run = 79516, bool useWeighted = false)
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TString phiType[3] = { "central", "west", "east" };

  const int nfiles = 7;
  oneDHist hists[nfiles]{};
  const TString weightType = useWeighted ? "weighted" : "unweighted";
  const TString outputDir = Form("Rootfiles/%s", weightType.Data());
  if (gSystem->mkdir(outputDir, true) != 0 && gSystem->AccessPathName(outputDir))
  {
    std::cerr << "Failed to create output directory " << outputDir << std::endl;
    return;
  }
  std::cout << "Use CPM " << weightType << " correction maps" << std::endl;
  TString names[nfiles] = 
  {
    "data_reco_FITacts",
    "data_reco_FITgenfit",
    "Lamination_CDB",
    "CM_Devon_reference75073",
    "CM_Devon_reference83319",
    "PHGARFIELD_Hugo",
    "PHGARFIELD_Hugo_new"
  };

  for (int itype = 0; itype < 3; itype++)
  {
    double selectPhi = phi_type2val(phiType[itype]);
    if (std::isnan(selectPhi)) return;

    std::cout << "Use phiType = " << phiType[itype]
          << ", phi = " << selectPhi
          << std::endl;

    std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
    for (const auto& selectZ : selectZs)
    {
      for (int i=0; i<nfiles; i++)
      {
        if (i<2)
        {
          draw1Dmap_atZ(
            Form("../../output/%s/run%d_CPM_%s_maxdca0p2_B3_average_correction_histograms.root",
                 names[i].Data(), run, weightType.Data()),
            Form("CPM_%s", names[i].Data()),
            selectZ, selectPhi, hists[i], kBlack);
        }
        else if (i==2)
        {
          auto rc = recoConsts::instance();
          rc->set_uint64Flag("TIMESTAMP", run);
          rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
          std::string cdbfilename_lamination = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
          cout<<"CDB file path: "<<cdbfilename_lamination.c_str()<<endl;

          draw1Dmap_2D(cdbfilename_lamination.c_str(), names[i], selectZ, selectPhi, hists[i], kBlack);
        }
	else if (i==3)
	{
	  //Devon's map
          draw1Dmap_2D("/gpfs/mnt/gpfs02/sphenix/user/dloomis/CentralMembraneAnalysis/CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079516_reference75073.root", names[i], selectZ, selectPhi, hists[i], kBlack);
	}
	else if (i==4)
	{
	  //Devon's map
          draw1Dmap_2D("/gpfs/mnt/gpfs02/sphenix/user/dloomis/CentralMembraneAnalysis/CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079516_reference83319.root", names[i], selectZ, selectPhi, hists[i], kBlack);
	}
	else if (i==5)
	{
	  //Hugo's map
          draw1Dmap_atZ("/phenix/u/hpereira/sphenix/work/g4simulations/Rootfiles/Distortions-00079513_beam_induced_phgarfield.root", names[i], selectZ, selectPhi, hists[i], kBlack);
	}
	else if (i==6)
	{
	  //Hugo's map
          draw1Dmap_atZ("/phenix/u/hpereira/sphenix/work/g4simulations/Rootfiles/Distortions-00079513_beam_induced_phgarfield_new.root", names[i], selectZ, selectPhi, hists[i], kBlack);
	}
      }

      TFile* ofile = new TFile(Form("%s/hist_1Dr_Z%d_%s.root", outputDir.Data(), (int)selectZ, phiType[itype].Data()), "recreate");
      ofile->cd();
      for (int i=0; i<nfiles; i++)
      {
        if (i<2) {hists[i].h_N_pos->Write(); hists[i].h_N_neg->Write();}
        hists[i].h_R_pos->Write(); hists[i].h_R_neg->Write();
        hists[i].h_P_pos->Write(); hists[i].h_P_neg->Write();
        hists[i].h_Z_pos->Write(); hists[i].h_Z_neg->Write();
      }
      ofile->Write();
      ofile->Close();
      delete ofile;

    }
    
    std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
    for (const auto& selectR : selectRs)
    {
      for (int i=0; i<nfiles; i++)
      {
        if (i<2)
        {
          draw1Dmap_atR(
            Form("../../output/%s/run%d_CPM_%s_maxdca0p2_B3_average_correction_histograms.root",
                 names[i].Data(), run, weightType.Data()),
            Form("CPM_%s", names[i].Data()),
            selectR, selectPhi, hists[i], kBlack);
        }
	else if (i==5)
	{
	  //Hugo's map
          draw1Dmap_atR("/phenix/u/hpereira/sphenix/work/g4simulations/Rootfiles/Distortions-00079513_beam_induced_phgarfield.root", names[i], selectR, selectPhi, hists[i], kBlack);
	}
	else if (i==6)
	{
	  //Hugo's map
          draw1Dmap_atR("/phenix/u/hpereira/sphenix/work/g4simulations/Rootfiles/Distortions-00079513_beam_induced_phgarfield_new.root", names[i], selectR, selectPhi, hists[i], kBlack);
	}
      }

      TFile* ofile = new TFile(Form("%s/hist_1Dz_R%d_%s.root", outputDir.Data(), (int)selectR, phiType[itype].Data()), "recreate");
      ofile->cd();
      for (int i=0; i<nfiles; i++)
      {
        if (i>=2 && i<=4) continue;
        hists[i].h_N_pos->Write(); hists[i].h_N_neg->Write();
        hists[i].h_R_pos->Write(); hists[i].h_R_neg->Write();
        hists[i].h_P_pos->Write(); hists[i].h_P_neg->Write();
        hists[i].h_Z_pos->Write(); hists[i].h_Z_neg->Write();
      }
      ofile->Write();
      ofile->Close();
      delete ofile;
    
    }

  }

}
