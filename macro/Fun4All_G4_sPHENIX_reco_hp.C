#include <caloreco/RawClusterBuilderTopo.h>

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>

#include <g4main/PHG4SimpleEventGenerator.h>
#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <trackreco/PHTrackPruner.h>
#include <PHCPMTpcCalibration.h>

/*
// own modules
#include <g4eval_hp/EventCounter_hp.h>
#include <g4eval_hp/SimEvaluator_hp.h>
#include <g4eval_hp/MicromegasEvaluator_hp.h>
#include <g4eval_hp/MicromegasTrackEvaluator_hp.h>
#include <g4eval_hp/TrackingEvaluator_hp.h>
*/

// local macros
#include <Calo_Calib.C>
#include <G4Setup_sPHENIX.C>
#include <G4_Global.C>
#include <QA.C>
#include <Trkr_QA.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_Reco.C>

R__LOAD_LIBRARY(libfun4all.so)
//R__LOAD_LIBRARY(libg4eval_hp.so)
R__LOAD_LIBRARY(libcpm.so)

//____________________________________________________________________
int Fun4All_G4_sPHENIX_reco_hp(
    const int nEvents = 10,
    const std::string inputFile = "~/hftg01/DST_FOR_DISTORTION/SimulationDST/G4Hits-00000.root",
    const std::string outdir = "root/",
    const std::string outfilename = "dst_sim_acts",
    const int index = 0,
    const int stepsize = 10,
    const int segment = 0,
    const bool writeMiniDst = true,
    const bool writePrunedSeedsToMiniDst = false)
{
  // print inputs
  std::cout << "Fun4All_G4_sPHENIX_reco_hp - nEvents: " << nEvents << std::endl;
  std::cout << "Fun4All_G4_sPHENIX_reco_hp - inputFile: " << inputFile << std::endl;

  // options
  Enable::PIPE = true;
  Enable::MBD = true;
  Enable::MBDRECO = true;
  Enable::MBDFAKE = false;
  Enable::PLUGDOOR = false;

  // enable all absorbers
  // this is equivalent to the old "absorberactive" flag
  Enable::ABSORBER = false;

  // central tracking
  Enable::MVTX = true;
  Enable::INTT = true;
  Enable::TPC = true;
  Enable::MICROMEGAS = true;
  Enable::BLACKHOLE = true;

  Enable::CEMC = false;
  Enable::HCALIN = false;
  Enable::MAGNET = true;
  Enable::HCALOUT = false;

  // deactivate all absorbers
  Enable::CEMC_ABSORBER = false;
  Enable::HCALIN_ABSORBER = false;
  Enable::MAGNET_ABSORBER = false;
  Enable::HCALOUT_ABSORBER = false;

  G4TPC::ENABLE_STATIC_DISTORTIONS = false;

  //G4TPC::ENABLE_STATIC_DISTORTIONS = true;
  G4TPC::DISTORTIONS_USE_PHI_AS_RADIANS = false;
  G4TPC::ENABLE_REACHES_READOUT = false;
  G4TPC::static_distortion_filename = "/phenix/u/hpereira/sphenix/work/g4simulations/distortion_maps/average_minus_static_distortion_converted.root";

  G4TPC::ENABLE_STATIC_CORRECTIONS = false;
  // G4TPC::correction_filename = "/phenix/u/hpereira/sphenix/work/g4simulations/distortion_maps/average_minus_static_distortion_inverted_10.root";

  // distortion reconstruction
  G4TRACKING::SC_CALIBMODE = true;
  G4TRACKING::SC_USE_MICROMEGAS = true;

  std::cout<< "Fun4All_CombinedDataReconstruction - tpc_drift_velocity_sim: " << G4TPC::tpc_drift_velocity_sim << std::endl;
  std::cout<< "Fun4All_CombinedDataReconstruction - tpc_drift_velocity_reco: " << G4TPC::tpc_drift_velocity_reco << std::endl;

  int runnumber=0;
  const std::string outputBase = outfilename + "_" + std::to_string(runnumber) + "-" + std::to_string(segment) + ".root";
  const std::string outDir = outdir + "/inReconstruction/" + std::to_string(runnumber) + "/";
  const std::string outputDirMove = outdir + "/Reconstructed/" + std::to_string(runnumber) + "/";
  const std::string makeDirectoryMove = "mkdir -p " + outputDirMove;
  std::string makeDirectory = "mkdir -p " + outDir;
  system(makeDirectory.c_str());
  TString outfile = outDir + outputBase;
  std::cout << "outfile " << outfile << std::endl;
  std::string theOutfile = outfile.Data();

  // server
  auto se = Fun4AllServer::instance();
  // se->Verbosity(2);

  // input manager
  auto in = new Fun4AllDstInputManager("DSTin");
  in->fileopen(inputFile);
  se->registerInputManager(in);

  // reco const
  auto rc = recoConsts::instance();

  // make sure to printout random seeds for reproducibility
  PHRandomSeed::Verbosity(1);

  // rc->set_IntFlag("RANDOMSEED",1);

  // condition database
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG",CDB::global_tag);
  rc->set_uint64Flag("TIMESTAMP",CDB::timestamp);

  // BBC
  //Mbd_Reco();

  // cells
  Mvtx_Cells();
  Intt_Cells();
  TPC_Cells();
  Micromegas_Cells();

  // tracking initialization
  TrackingInit();

  // digitizer and clustering
  Mvtx_Clustering();
  Intt_Clustering();
  TPC_Clustering();
  Micromegas_Clustering();

  // silicon seeding
  auto silicon_Seeding = new PHActsSiliconSeeding;
  silicon_Seeding->Verbosity(0);
  silicon_Seeding->setStrobeRange(-5,5);
  silicon_Seeding->seedAnalysis(false);
  silicon_Seeding->setinttRPhiSearchWindow(0.2);
  silicon_Seeding->setinttZSearchWindow(1.0);
  se->registerSubsystem(silicon_Seeding);

  auto merger = new PHSiliconSeedMerger;
  merger->Verbosity(0);
  se->registerSubsystem(merger);

  // TPC seeding
  auto seeder = new PHCASeeding("PHCASeeding");
  seeder->Verbosity(0);
  seeder->SetLayerRange(7, 55);
  seeder->SetSearchWindow(2.,0.05); // z-width and phi-width, default in macro at 1.5 and 0.05
  seeder->SetClusAdd_delta_window(3.0,0.06); //  (0.5, 0.005) are default; sdzdr_cutoff, d2/dr2(phi)_cutoff
  seeder->SetMinHitsPerCluster(0);
  seeder->SetMinClustersPerTrack(3);
  seeder->useFixedClusterError(true);
  se->registerSubsystem(seeder);

  // expand stubs in the TPC using simple kalman filter
  auto cprop = new PHSimpleKFProp("PHSimpleKFProp");
  cprop->useFixedClusterError(true);
  cprop->set_max_window(5.);
  cprop->Verbosity(0);
  se->registerSubsystem(cprop);

  // matching to silicons
  auto silicon_match = new PHSiliconTpcTrackMatching;
  silicon_match->Verbosity(0);

  // narrow matching windows
  silicon_match->set_x_search_window(0.36);
  silicon_match->set_y_search_window(0.36);
  silicon_match->set_z_search_window(2.5);
  silicon_match->set_phi_search_window(0.014);
  silicon_match->set_eta_search_window(0.0091);
  silicon_match->set_test_windows_printout(false);
  se->registerSubsystem(silicon_match);

  // matching with micromegas
  auto mm_match = new PHMicromegasTpcTrackMatching;
  mm_match->Verbosity(0);
  mm_match->set_rphi_search_window_lyr1(3.0);
  mm_match->set_rphi_search_window_lyr2(15.0);

  mm_match->set_z_search_window_lyr1(30.0);
  mm_match->set_z_search_window_lyr2(3.0);

  mm_match->set_min_tpc_layer(38);             // layer in TPC to start projection fit
  mm_match->set_test_windows_printout(false);  // used for tuning search windows only
  se->registerSubsystem(mm_match);

  // track fit
  se->registerSubsystem(new PHTpcDeltaZCorrection);

  // perform final track fit with ACTS
  auto actsFit = new PHActsTrkFitter;
  actsFit->Verbosity(0);
  actsFit->commissioning(G4TRACKING::use_alignment);

  // fit with Micromegas and Silicon ONLY
  actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
  actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);

  actsFit->set_use_clustermover(true);
  actsFit->useActsEvaluator(false);
  actsFit->useOutlierFinder(false);
  actsFit->setFieldMap(G4MAGNET::magfield_tracking);
  se->registerSubsystem(actsFit);

  auto cleaner = new PHTrackCleaner();
  cleaner->Verbosity(0);
  se->registerSubsystem(cleaner);

  //prune acts full tracks, create new SvtxTrackMap
  auto trackpruner = new PHTrackPruner;
  trackpruner->Verbosity(0);
  trackpruner->set_pruned_svtx_seed_map_name("PrunedSvtxTrackSeedContainer");
  trackpruner->set_track_pt_low_cut(0.5);
  trackpruner->set_track_quality_high_cut(100);
  trackpruner->set_nmvtx_clus_low_cut(3);
  trackpruner->set_nintt_clus_low_cut(2);
  trackpruner->set_ntpc_clus_low_cut(35);
  trackpruner->set_ntpot_clus_low_cut(1);
  trackpruner->set_nmvtx_states_low_cut(3);
  trackpruner->set_nintt_states_low_cut(2);
  trackpruner->set_ntpc_states_low_cut(35);
  trackpruner->set_ntpot_states_low_cut(1);
  se->registerSubsystem(trackpruner);

  // perform final track fit with ACTS
  // Si-TPOT fit
  auto actsFit_SiTpotFit = new PHActsTrkFitter;
  actsFit_SiTpotFit->Verbosity(0);
  actsFit_SiTpotFit->commissioning(G4TRACKING::use_alignment);
  // in calibration mode, fit only Silicons and Micromegas hits
  actsFit_SiTpotFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
  actsFit_SiTpotFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
  actsFit_SiTpotFit->set_svtx_seed_map_name("PrunedSvtxTrackSeedContainer");
  actsFit_SiTpotFit->set_pp_mode(TRACKING::pp_mode);
  actsFit_SiTpotFit->set_use_clustermover(true);  // default is true for now
  actsFit_SiTpotFit->useActsEvaluator(false);
  actsFit_SiTpotFit->useOutlierFinder(false);
  actsFit_SiTpotFit->setFieldMap(G4MAGNET::magfield_tracking);
  se->registerSubsystem(actsFit_SiTpotFit);

  std::string cpmstring;
  std::string cpmmindststring;
  std::string cpmmindstfinalstring;
  if (G4TRACKING::SC_CALIBMODE)
  {
    auto cpmreco = new PHCPMTpcCalibration;
    const TString cpmoutfile = theOutfile + "_CPMVoxelContainer.root";
    cpmstring = cpmoutfile.Data();
    cpmmindststring = theOutfile + "_cpm_mini_dst.root";
    cpmmindstfinalstring = outputDirMove + outputBase + "_cpm_mini_dst.root";

    cpmreco->setOutputfile(cpmstring);
    //cpmreco->setClusterSource(inputclusterFile);//cluster is generated on the fly
    cpmreco->setTrackSource(writeMiniDst ? cpmmindstfinalstring : "");
    cpmreco->setRunSegment(runnumber, segment);
    cpmreco->setTrackMapName("SvtxSiliconMMTrackMap");
    cpmreco->setWriteRecords(true);
    cpmreco->setWriteQARecords(true);
    cpmreco->setMinPt(0.5);
    cpmreco->requireCrossing(false);
    cpmreco->requireTPOT(true);
    // CPM does not apply the legacy PHTpcResiduals central-membrane requirement.
    cpmreco->disableAverageCorr();
    cpmreco->setGridDimensions(36, 16, 80);
    se->registerSubsystem(cpmreco);

    if (!writeMiniDst)
    {
      std::cout << "Fun4All_G4_sPHENIX_reco_hp - writeMiniDst is false. "
                << "CPM snapshots remain usable, but SvtxTrack object rehydration is disabled."
                << std::endl;
    }
    if (writeMiniDst)
    {
      auto out = new Fun4AllDstOutputManager("CPMMiniDstOutput", cpmmindststring);
      out->AddNode("Sync");
      out->AddNode("EventHeader");
      out->AddNode("SvtxSiliconMMTrackMap");
      if (writePrunedSeedsToMiniDst)
      {
        out->AddNode("PrunedSvtxTrackSeedContainer");
      }
      se->registerOutputManager(out);
    }
  }

  Enable::QA = true;
  if (Enable::QA)
  {
    Distortions_QA();
  }

  // process events
  se->skip(stepsize * index);
  se->run(nEvents);

  // terminate
  se->End();
  se->PrintTimer();
  CDBInterface::instance()->Print();

  std::string qaOutputFileName;
  if (Enable::QA)
  {
    TString qaname = theOutfile + "_qa.root";
    qaOutputFileName = qaname.Data();
    QAHistManagerDef::saveQARootFile(qaOutputFileName);
  }

  std::ifstream file_cpm(cpmstring.c_str(), std::ios::binary | std::ios::ate);
  if (file_cpm.good() && (file_cpm.tellg() > 100))
  {
    system(makeDirectoryMove.c_str());
    std::string moveOutput = "mv " + cpmstring + " " + outputDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  std::ifstream file_cpmmindst(cpmmindststring.c_str(), std::ios::binary | std::ios::ate);
  if (file_cpmmindst.good() && (file_cpmmindst.tellg() > 100))
  {
    system(makeDirectoryMove.c_str());
    std::string moveOutput = "mv " + cpmmindststring + " " + outputDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  std::ifstream file_qa(qaOutputFileName.c_str(), std::ios::binary | std::ios::ate);
  if (file_qa.good() && (file_qa.tellg() > 100))
  {
    system(makeDirectoryMove.c_str());
    std::string moveOutput = "mv " + qaOutputFileName + " " + outputDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
