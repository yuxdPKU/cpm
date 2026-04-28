/*
 * CPM Job A Fun4All macro.
 *
 * This macro follows currentworkflow/Fun4All_TrackAnalysis.C and replaces the
 * PHTpcResiduals average-correction extraction block with PHCPMTpcCalibration.
 */

#include <fun4all/Fun4AllUtils.h>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <GlobalVariables.C>
#include <QA.C>
#include <Trkr_QA.C>
#include <Trkr_Clustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <cdbobjects/CDBTTree.h>

#include <PHCPMTpcCalibration.h>

#include <trackingqa/TpcSeedsQA.h>

#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>

// #include <distortionanalysis/DistortionAnalysis.h>
#include <trackreco/PHTrackPruner.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <utility>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libDistortionAnalysis.so)
R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libtpcqa.so)
R__LOAD_LIBRARY(libcpm.so)

void Fun4All_CPMTrackAnalysis(
    const int nEvents = 10,
    const std::string clusterfilename = "DST_TRKR_CLUSTER_run3pp_ana532_2025p009_v001-00079516-00000.root",
    const std::string outdir = "root/",
    const std::string outfilename = "clusters_seeds",
    const int index = 0,
    const int stepsize = 10,
    const bool convertSeeds = false,
    const bool writeMiniDst = true)
{
  std::string inputclusterFile = clusterfilename;

  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
  std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;

  std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(clusterfilename);
  int runnumber = runseg.first;
  int segment = runseg.second;

  auto rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  TpcReadoutInit(runnumber);
  // these lines show how to override the drift velocity and time offset values set in TpcReadoutInit
  // G4TPC::tpc_drift_velocity_reco = 0.0073844; // cm/ns
  // TpcClusterZCrossingCorrection::_vdrift = G4TPC::tpc_drift_velocity_reco;
  // G4TPC::tpc_tzero_reco = -5*50;  // ns
  std::cout << " run: " << runnumber
            << " samples: " << TRACKING::reco_tpc_maxtime_sample
            << " pre: " << TRACKING::reco_tpc_time_presample
            << " vdrift: " << G4TPC::tpc_drift_velocity_reco
            << std::endl;

  // distortion calibration mode
  /*
   * set to true to enable residuals in the TPC with
   * TPC clusters not participating to the ACTS track fit
   */
  G4TRACKING::SC_CALIBMODE = true;
  G4TRACKING::SC_USE_MICROMEGAS = true;
  TRACKING::pp_mode = true;

  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  const std::string outputBase = outfilename + "_" + std::to_string(runnumber) + "-" + std::to_string(segment) + ".root";
  const std::string outDir = outdir + "/inReconstruction/" + std::to_string(runnumber) + "/";
  const std::string outputDirMove = outdir + "/Reconstructed/" + std::to_string(runnumber) + "/";
  const std::string makeDirectoryMove = "mkdir -p " + outputDirMove;
  std::string makeDirectory = "mkdir -p " + outDir;
  system(makeDirectory.c_str());
  TString outfile = outDir + outputBase;
  std::cout << "outfile " << outfile << std::endl;
  std::string theOutfile = outfile.Data();

  auto se = Fun4AllServer::instance();
  se->Verbosity(1);

  auto ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;

  //to turn on the default static corrections, enable the two lines below
  G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;

  //to turn on the average corrections, enable the three lines below
  //note: these are designed to be used only if static corrections are also applied
  G4TPC::ENABLE_AVERAGE_CORRECTIONS = true;
   // to use a custom file instead of the database file:
  G4TPC::average_correction_filename = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
  std::cout << "Average distortion map used: " << G4TPC::average_correction_filename << std::endl;

  G4MAGNET::magfield_rescale = 1;
  TrackingInit();

  auto hitsinclus = new Fun4AllDstInputManager("ClusterInputManager");
  hitsinclus->fileopen(inputclusterFile);
  se->registerInputManager(hitsinclus);

  /*
   * Begin Track Seeding
   */

  /*
   * Silicon Seeding
   */

  /*
  auto silicon_Seeding = new PHActsSiliconSeeding;
  silicon_Seeding->Verbosity(0);
  silicon_Seeding->searchInIntt();
  silicon_Seeding->setinttRPhiSearchWindow(0.4);
  silicon_Seeding->setinttZSearchWindow(1.6);
  silicon_Seeding->seedAnalysis(false);
  se->registerSubsystem(silicon_Seeding);
  */

  auto silicon_Seeding = new PHActsSiliconSeeding;
  silicon_Seeding->Verbosity(0);
  silicon_Seeding->setStrobeRange(-5,5);
  // these get us to about 83% INTT > 1
  silicon_Seeding->setinttRPhiSearchWindow(0.4);
  silicon_Seeding->setinttZSearchWindow(2.0);
  silicon_Seeding->seedAnalysis(false);
  se->registerSubsystem(silicon_Seeding);

  auto merger = new PHSiliconSeedMerger;
  merger->Verbosity(0);
  se->registerSubsystem(merger);

  /*
   * Tpc Seeding
   */
  auto seeder = new PHCASeeding("PHCASeeding");
  double fieldstrength = std::numeric_limits<double>::quiet_NaN();  // set by isConstantField if constant
  bool ConstField = isConstantField(G4MAGNET::magfield_tracking, fieldstrength);
  if (ConstField)
  {
    seeder->useConstBField(true);
    seeder->constBField(fieldstrength);
  }
  else
  {
    seeder->set_field_dir(-1 * G4MAGNET::magfield_rescale);
    seeder->useConstBField(false);
    seeder->magFieldFile(G4MAGNET::magfield_tracking);  // to get charge sign right
  }
  seeder->Verbosity(0);
  seeder->SetLayerRange(7, 55);
  seeder->SetSearchWindow(2.,0.05); // z-width and phi-width, default in macro at 1.5 and 0.05
  seeder->SetClusAdd_delta_window(3.0,0.06); //  (0.5, 0.005) are default; sdzdr_cutoff, d2/dr2(phi)_cutoff
  //seeder->SetNClustersPerSeedRange(4,60); // default is 6, 6
  seeder->SetMinHitsPerCluster(0);
  seeder->SetMinClustersPerTrack(3);
  seeder->useFixedClusterError(true);
  seeder->set_pp_mode(true);
  se->registerSubsystem(seeder);

  // expand stubs in the TPC using simple kalman filter
  auto cprop = new PHSimpleKFProp("PHSimpleKFProp");
  cprop->set_field_dir(G4MAGNET::magfield_rescale);
  if (ConstField)
  {
    cprop->useConstBField(true);
    cprop->setConstBField(fieldstrength);
  }
  else
  {
    cprop->magFieldFile(G4MAGNET::magfield_tracking);
    cprop->set_field_dir(-1 * G4MAGNET::magfield_rescale);
  }
  cprop->useFixedClusterError(true);
  cprop->set_max_window(5.);
  cprop->Verbosity(0);
  cprop->set_pp_mode(true);
  se->registerSubsystem(cprop);

  // Always apply preliminary distortion corrections to TPC clusters before silicon matching
  // and refit the trackseeds. Replace KFProp fits with the new fit parameters in the TPC seeds.
  auto prelim_distcorr = new PrelimDistortionCorrection;
  prelim_distcorr->set_pp_mode(true);
  prelim_distcorr->Verbosity(0);
  se->registerSubsystem(prelim_distcorr);

  /*
   * Track Matching between silicon and TPC
   */
  // The normal silicon association methods
  // Match the TPC track stubs from the CA seeder to silicon track stubs from PHSiliconTruthTrackSeeding
  auto silicon_match = new PHSiliconTpcTrackMatching;
  silicon_match->Verbosity(0);
  silicon_match->set_pp_mode(TRACKING::pp_mode);
  if (G4TPC::ENABLE_AVERAGE_CORRECTIONS)
  {
    // for general tracking
    // Eta/Phi window is determined by 3 sigma window
    // X/Y/Z window is determined by 4 sigma window
    silicon_match->window_deta.set_posQoverpT_maxabs({-0.014,0.0331,0.48});
    silicon_match->window_deta.set_negQoverpT_maxabs({-0.006,0.0235,0.52});
    silicon_match->set_deltaeta_min(0.03);
    silicon_match->window_dphi.set_QoverpT_range({-0.15,0,0}, {0.15,0,0});
    silicon_match->window_dx.set_QoverpT_maxabs({3.0,0,0});
    silicon_match->window_dy.set_QoverpT_maxabs({3.0,0,0});
    silicon_match->window_dz.set_posQoverpT_maxabs({1.138,0.3919,0.84});
    silicon_match->window_dz.set_negQoverpT_maxabs({0.719,0.6485,0.65});
    silicon_match->set_crossing_deltaz_max(30);
    silicon_match->set_crossing_deltaz_min(2);

    // for distortion correction using SI-TPOT fit and track pT>0.5
    if (G4TRACKING::SC_CALIBMODE)
    {
      silicon_match->window_deta.set_posQoverpT_maxabs({0.016,0.0060,1.13});
      silicon_match->window_deta.set_negQoverpT_maxabs({0.022,0.0022,1.44});
      silicon_match->set_deltaeta_min(0.03);
      silicon_match->window_dphi.set_QoverpT_range({-0.15,0,0}, {0.09,0,0});
      silicon_match->window_dx.set_QoverpT_maxabs({2.0,0,0});
      silicon_match->window_dy.set_QoverpT_maxabs({1.5,0,0});
      silicon_match->window_dz.set_posQoverpT_maxabs({1.213,0.0211,2.09});
      silicon_match->window_dz.set_negQoverpT_maxabs({1.307,0.0001,4.52});
      silicon_match->set_crossing_deltaz_min(1.2);
    }
  }
  silicon_match->print_windows(true);
  se->registerSubsystem(silicon_match);

  // Match TPC track stubs from CA seeder to clusters in the micromegas layers
  auto mm_match = new PHMicromegasTpcTrackMatching;
  mm_match->Verbosity(0);
  mm_match->set_pp_mode(TRACKING::pp_mode);
  //mm_match->set_rphi_search_window_lyr1(3.);
  mm_match->set_rphi_search_window_lyr1(1.5);//test value
  mm_match->set_rphi_search_window_lyr2(15.0);
  mm_match->set_z_search_window_lyr1(30.0);
  mm_match->set_z_search_window_lyr2(3.);

  mm_match->set_min_tpc_layer(38);             // layer in TPC to start projection fit
  mm_match->set_test_windows_printout(false);  // used for tuning search windows only
  se->registerSubsystem(mm_match);

  /*
   * Either converts seeds to tracks with a straight line/helix fit
   * or run the full Acts track kalman filter fit
   */
  if (G4TRACKING::convert_seeds_to_svtxtracks)
  {
    auto converter = new TrackSeedTrackMapConverter;
    // Default set to full SvtxTrackSeeds. Can be set to
    // SiliconTrackSeedContainer or TpcTrackSeedContainer
    converter->setTrackSeedName("SvtxTrackSeedContainer");
    converter->setFieldMap(G4MAGNET::magfield_tracking);
    converter->Verbosity(0);
    se->registerSubsystem(converter);
  }
  else
  {
    auto deltazcorr = new PHTpcDeltaZCorrection;
    deltazcorr->Verbosity(0);
    se->registerSubsystem(deltazcorr);

    // perform final track fit with ACTS
    auto actsFit = new PHActsTrkFitter;
    actsFit->Verbosity(0);
    actsFit->commissioning(G4TRACKING::use_alignment);
    // in calibration mode, fit only Silicons and Micromegas hits
//    actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
//    actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    actsFit->fitSiliconMMs(false);//by default
    actsFit->setUseMicromegas(true);//by default
    actsFit->set_pp_mode(TRACKING::pp_mode);
    actsFit->set_use_clustermover(true);  // default is true for now
    actsFit->useActsEvaluator(false);
    actsFit->useOutlierFinder(false);
    actsFit->setFieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(actsFit);

    auto cleaner = new PHTrackCleaner();
    cleaner->Verbosity(0);
    cleaner->set_pp_mode(TRACKING::pp_mode);
    se->registerSubsystem(cleaner);

  }

  auto finder = new PHSimpleVertexFinder;
  finder->Verbosity(0);
  
  //new cuts
  finder->setDcaCut(0.05);
  finder->setTrackPtCut(0.1);
  finder->setBeamLineCut(1);
  finder->setTrackQualityCut(300);
  finder->setNmvtxRequired(3);
  finder->setOutlierPairCut(0.10);
  
  se->registerSubsystem(finder);

  // Propagate track positions to the vertex position
  auto vtxProp = new PHActsVertexPropagator;
  vtxProp->Verbosity(0);
  vtxProp->fieldMap(G4MAGNET::magfield_tracking);
  se->registerSubsystem(vtxProp);
 
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
    auto cpmreco = new cpm::PHCPMTpcCalibration;
    const TString cpmoutfile = theOutfile + "_CPMVoxelContainer.root";
    cpmstring = cpmoutfile.Data();
    cpmmindststring = theOutfile + "_cpm_mini_dst.root";
    cpmmindstfinalstring = outputDirMove + outputBase + "_cpm_mini_dst.root";

    cpmreco->setOutputfile(cpmstring);
    cpmreco->setClusterSource(inputclusterFile);
    cpmreco->setTrackSource(writeMiniDst ? cpmmindstfinalstring : inputclusterFile);
    cpmreco->setRunSegment(runnumber, segment);
    cpmreco->setTrackMapName("SvtxSiliconMMTrackMap");
    cpmreco->setMinPt(0.5);
    cpmreco->requireCrossing(false);
    cpmreco->requireTPOT(true);
    cpmreco->requireCM(true);
    cpmreco->disableAverageCorr();
    cpmreco->setGridDimensions(36, 16, 80);
    se->registerSubsystem(cpmreco);

    if (writeMiniDst)
    {
      auto out = new Fun4AllDstOutputManager("CPMMiniDstOutput", cpmmindststring);
      out->AddNode("Sync");
      out->AddNode("EventHeader");
      //out->AddNode("TRKR_CLUSTER");
      out->AddNode("SvtxSiliconMMTrackMap");
      out->AddNode("PrunedSvtxTrackSeedContainer");
      se->registerOutputManager(out);
    }
  }

  Enable::QA = true;
  if (Enable::QA)
  {
    Distortions_QA();
  }

  se->skip(stepsize * index);
  se->run(nEvents);
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

  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}
