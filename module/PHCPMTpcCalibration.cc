#include "PHCPMTpcCalibration.h"

#include "CPMPairUtils.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <ffaobjects/EventHeader.h>
#include <ffaobjects/SyncDefs.h>
#include <ffaobjects/SyncObject.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>

#include <TFile.h>
#include <TTree.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <array>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>

namespace
{
  struct CPMRecordTreeFields
  {
    std::string cluster_source;
    std::string track_source;
    std::string track_map_name;

    int run = -1;
    int segment = -1;
    int sync_event = -1;
    int event_sequence = -1;
    unsigned long long stream_event_ordinal = 0;

    unsigned int track_id = InvalidTrackId;
    int charge = 0;
    float pt = std::numeric_limits<float>::quiet_NaN();
    float quality = std::numeric_limits<float>::quiet_NaN();
    unsigned short n_mvtx = 0;
    unsigned short n_intt = 0;
    unsigned short n_tpc = 0;
    unsigned short n_tpot = 0;
    unsigned short n_mvtx_states = 0;
    unsigned short n_intt_states = 0;
    unsigned short n_tpc_states = 0;
    unsigned short n_tpot_states = 0;

    unsigned long long cluskey = InvalidClusterKey;
    unsigned long long hitsetkey = InvalidHitSetKey;
    unsigned int subsurfkey = InvalidSubSurfKey;
    unsigned short layer = 0;
    unsigned short side = 0;

    int iphi = -1;
    int ir = -1;
    int iz = -1;

    double cluster_x = std::numeric_limits<double>::quiet_NaN();
    double cluster_y = std::numeric_limits<double>::quiet_NaN();
    double cluster_z = std::numeric_limits<double>::quiet_NaN();
    double voxel_x = std::numeric_limits<double>::quiet_NaN();
    double voxel_y = std::numeric_limits<double>::quiet_NaN();
    double voxel_z = std::numeric_limits<double>::quiet_NaN();
    double offset_x = std::numeric_limits<double>::quiet_NaN();
    double offset_y = std::numeric_limits<double>::quiet_NaN();
    double offset_z = std::numeric_limits<double>::quiet_NaN();

    double state_pathlength = std::numeric_limits<double>::quiet_NaN();
    double state_local_x = std::numeric_limits<double>::quiet_NaN();
    double state_local_y = std::numeric_limits<double>::quiet_NaN();
    double state_x = std::numeric_limits<double>::quiet_NaN();
    double state_y = std::numeric_limits<double>::quiet_NaN();
    double state_z = std::numeric_limits<double>::quiet_NaN();
    double state_px = std::numeric_limits<double>::quiet_NaN();
    double state_py = std::numeric_limits<double>::quiet_NaN();
    double state_pz = std::numeric_limits<double>::quiet_NaN();
    double state_covariance[36] = {};

    bool has_crossing = false;
    bool passes_tpot = false;
    bool passes_track_quality = false;
    bool passes_geometry = false;

    void copy_from(const TrackStateRecord& record)
    {
      cluster_source = record.event_ref.cluster_source;
      track_source = record.event_ref.track_source;
      track_map_name = record.track_ref.track_map_name;

      run = record.event_ref.run;
      segment = record.event_ref.segment;
      sync_event = record.event_ref.sync_event;
      event_sequence = record.event_ref.event_sequence;
      stream_event_ordinal = record.event_ref.stream_event_ordinal;

      track_id = record.track_ref.track_id;
      charge = record.track.charge;
      pt = record.track.pt;
      quality = record.track.quality;
      n_mvtx = record.track.n_mvtx;
      n_intt = record.track.n_intt;
      n_tpc = record.track.n_tpc;
      n_tpot = record.track.n_tpot;
      n_mvtx_states = record.track.n_mvtx_states;
      n_intt_states = record.track.n_intt_states;
      n_tpc_states = record.track.n_tpc_states;
      n_tpot_states = record.track.n_tpot_states;

      cluskey = record.cluster_ref.cluskey;
      hitsetkey = record.cluster_ref.hitsetkey;
      subsurfkey = record.cluster_ref.subsurfkey;
      layer = record.cluster_ref.layer;
      side = record.cluster_ref.side;

      iphi = record.voxel.iphi;
      ir = record.voxel.ir;
      iz = record.voxel.iz;

      cluster_x = record.cluster.corrected_position.x;
      cluster_y = record.cluster.corrected_position.y;
      cluster_z = record.cluster.corrected_position.z;
      voxel_x = record.cluster.voxel_center.x;
      voxel_y = record.cluster.voxel_center.y;
      voxel_z = record.cluster.voxel_center.z;
      offset_x = record.cluster.cluster_minus_voxel_center.x;
      offset_y = record.cluster.cluster_minus_voxel_center.y;
      offset_z = record.cluster.cluster_minus_voxel_center.z;

      state_pathlength = record.state.pathlength;
      state_local_x = record.state.local_x;
      state_local_y = record.state.local_y;
      state_x = record.state.position.x;
      state_y = record.state.position.y;
      state_z = record.state.position.z;
      state_px = record.state.momentum.x;
      state_py = record.state.momentum.y;
      state_pz = record.state.momentum.z;
      for (std::size_t i = 0; i < 36; ++i)
      {
        state_covariance[i] = record.state.covariance[i];
      }

      has_crossing = record.selection.has_crossing;
      passes_tpot = record.selection.passes_tpot;
      passes_track_quality = record.selection.passes_track_quality;
      passes_geometry = record.selection.passes_geometry;
    }
  };

  void book_cpm_record_tree(TTree& tree, CPMRecordTreeFields& fields)
  {
    tree.Branch("cluster_source", &fields.cluster_source);
    tree.Branch("track_source", &fields.track_source);
    tree.Branch("track_map_name", &fields.track_map_name);

    tree.Branch("run", &fields.run);
    tree.Branch("segment", &fields.segment);
    tree.Branch("sync_event", &fields.sync_event);
    tree.Branch("event_sequence", &fields.event_sequence);
    tree.Branch("stream_event_ordinal", &fields.stream_event_ordinal);

    tree.Branch("track_id", &fields.track_id);
    tree.Branch("charge", &fields.charge);
    tree.Branch("pt", &fields.pt);
    tree.Branch("quality", &fields.quality);
    tree.Branch("n_mvtx", &fields.n_mvtx);
    tree.Branch("n_intt", &fields.n_intt);
    tree.Branch("n_tpc", &fields.n_tpc);
    tree.Branch("n_tpot", &fields.n_tpot);
    tree.Branch("n_mvtx_states", &fields.n_mvtx_states);
    tree.Branch("n_intt_states", &fields.n_intt_states);
    tree.Branch("n_tpc_states", &fields.n_tpc_states);
    tree.Branch("n_tpot_states", &fields.n_tpot_states);

    tree.Branch("cluskey", &fields.cluskey);
    tree.Branch("hitsetkey", &fields.hitsetkey);
    tree.Branch("subsurfkey", &fields.subsurfkey);
    tree.Branch("layer", &fields.layer);
    tree.Branch("side", &fields.side);

    tree.Branch("iphi", &fields.iphi);
    tree.Branch("ir", &fields.ir);
    tree.Branch("iz", &fields.iz);

    tree.Branch("cluster_x", &fields.cluster_x);
    tree.Branch("cluster_y", &fields.cluster_y);
    tree.Branch("cluster_z", &fields.cluster_z);
    tree.Branch("voxel_x", &fields.voxel_x);
    tree.Branch("voxel_y", &fields.voxel_y);
    tree.Branch("voxel_z", &fields.voxel_z);
    tree.Branch("offset_x", &fields.offset_x);
    tree.Branch("offset_y", &fields.offset_y);
    tree.Branch("offset_z", &fields.offset_z);

    tree.Branch("state_pathlength", &fields.state_pathlength);
    tree.Branch("state_local_x", &fields.state_local_x);
    tree.Branch("state_local_y", &fields.state_local_y);
    tree.Branch("state_x", &fields.state_x);
    tree.Branch("state_y", &fields.state_y);
    tree.Branch("state_z", &fields.state_z);
    tree.Branch("state_px", &fields.state_px);
    tree.Branch("state_py", &fields.state_py);
    tree.Branch("state_pz", &fields.state_pz);
    tree.Branch("state_covariance", fields.state_covariance, "state_covariance[36]/D");

    tree.Branch("has_crossing", &fields.has_crossing);
    tree.Branch("passes_tpot", &fields.passes_tpot);
    tree.Branch("passes_track_quality", &fields.passes_track_quality);
    tree.Branch("passes_geometry", &fields.passes_geometry);
  }
}

PHCPMTpcCalibration::PHCPMTpcCalibration(const std::string& name)
  : SubsysReco(name)
{
}

int PHCPMTpcCalibration::Init(PHCompositeNode* /*topNode*/)
{
  std::cout << "PHCPMTpcCalibration::Init"
            << " outputfile: " << m_outputfile
            << " output object: " << m_outputObjectName
            << " trackmap: " << m_trackmapname
            << " grid: (" << m_phiBins << ", " << m_rBins << ", " << m_zBins << ")"
            << " running batch size: " << m_runningBatchSize
            << " write records: " << m_writeRecords
            << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCPMTpcCalibration::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_zMax = m_tGeometry->get_max_driftlength() + m_tGeometry->get_CM_halfwidth();
  m_zMin = -m_zMax;
  m_correctionContainer.set_grid_dimensions(m_phiBins, m_rBins, m_zBins);
  m_correctionContainer.set_grid_range(m_rMin, m_rMax, m_zMin, m_zMax);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCPMTpcCalibration::process_event(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  const int returnValue = processTracks();
  ++m_event;
  return returnValue;
}

int PHCPMTpcCalibration::End(PHCompositeNode* /*topNode*/)
{
  flushPendingBatches();

  std::cout << "PHCPMTpcCalibration::End"
            << " records: " << m_voxelContainer.record_count()
            << " voxels: " << m_voxelContainer.voxel_count()
            << " write_records: " << m_writeRecords
            << " batches: " << m_batches
            << " candidate pairs: " << m_candidate_pairs
            << " accepted pairs: " << m_accepted_pairs
            << " final flush batches: " << m_final_flush_batches
            << " outputfile: " << m_outputfile
            << std::endl;

  std::cout << "PHCPMTpcCalibration::End"
            << " track statistics total: " << m_total_tracks
            << " accepted: " << m_accepted_tracks
            << std::endl;

  std::cout << "PHCPMTpcCalibration::End"
            << " state statistics total: " << m_total_states
            << " accepted: " << m_accepted_states
            << std::endl;

  std::cout << "PHCPMTpcCalibration::End"
            << " final flush records positive: " << m_final_flush_positive_records
            << " negative: " << m_final_flush_negative_records
            << " unflushed single-charge positive: " << m_unflushed_positive_records
            << " negative: " << m_unflushed_negative_records
            << std::endl;

  return writeOutput();
}

void PHCPMTpcCalibration::setGridDimensions(const int phiBins, const int rBins, const int zBins)
{
  m_phiBins = phiBins;
  m_rBins = rBins;
  m_zBins = zBins;
}

int PHCPMTpcCalibration::getNodes(PHCompositeNode* topNode)
{
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE << "No TRKR_CLUSTER node on node tree. Exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "ActsGeometry not on node tree. Exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);
  if (!m_trackMap)
  {
    std::cout << PHWHERE << " " << m_trackmapname << " not on node tree. Exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_syncObject = findNode::getClass<SyncObject>(topNode, syncdefs::SYNCNODENAME);
  m_eventHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");

  m_globalPositionWrapper.loadNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCPMTpcCalibration::processTracks()
{
  if (Verbosity())
  {
    std::cout << "PHCPMTpcCalibration::processTracks - track map size "
              << m_trackMap->size() << std::endl;
  }

  for (const auto& [trackKey, track] : *m_trackMap)
  {
    ++m_total_tracks;
    if (!checkTrack(track))
    {
      continue;
    }

    ++m_accepted_tracks;

    std::map<VoxelId, TrackStateRecord> bestRecordsByVoxel;
    for (auto stateIter = track->begin_states(); stateIter != track->end_states(); ++stateIter)
    {
      const auto* state = stateIter->second;
      ++m_total_states;

      if (!checkState(state))
      {
        continue;
      }

      const auto cluskey = state->get_cluskey();
      auto* cluster = m_clusterContainer->findCluster(cluskey);
      if (!cluster)
      {
        continue;
      }

      const auto crossing = track->get_crossing();
      const auto actsPosition = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(
          cluskey, cluster, crossing);

      const Vector3 clusterPosition{
          actsPosition.x(),
          actsPosition.y(),
          actsPosition.z()};

      VoxelId voxel;
      if (!getVoxelId(clusterPosition, voxel))
      {
        continue;
      }

      auto record = makeRecord(trackKey, track, state, cluster, clusterPosition, voxel);
      auto [iter, inserted] = bestRecordsByVoxel.emplace(voxel, record);
      if (!inserted && isCloserToVoxelCenter(record, iter->second))
      {
        iter->second = record;
      }
    }

    for (auto& [voxel, record] : bestRecordsByVoxel)
    {
      (void) voxel;
      addRecord(std::move(record));
      ++m_accepted_states;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHCPMTpcCalibration::addRecord(TrackStateRecord record)
{
  const auto voxel = record.voxel;
  if (m_writeRecords)
  {
    m_voxelContainer.add(record);
  }

  if (record.track.charge > 0)
  {
    m_pendingRecords[voxel].positive.push_back(std::move(record));
  }
  else if (record.track.charge < 0)
  {
    m_pendingRecords[voxel].negative.push_back(std::move(record));
  }
  else
  {
    return;
  }

  processPendingBatches(voxel);
}

void PHCPMTpcCalibration::processPendingBatches(const VoxelId& voxel)
{
  if (m_runningBatchSize == 0)
  {
    return;
  }

  auto iter = m_pendingRecords.find(voxel);
  if (iter == m_pendingRecords.end())
  {
    return;
  }

  auto& pending = iter->second;
  if (pending.positive.size() < m_runningBatchSize ||
      pending.negative.size() < m_runningBatchSize)
  {
    return;
  }

  processBatch(voxel);
  pending.positive.clear();
  pending.negative.clear();
}

void PHCPMTpcCalibration::flushPendingBatches()
{
  m_final_flush_batches = 0;
  m_final_flush_positive_records = 0;
  m_final_flush_negative_records = 0;
  m_unflushed_positive_records = 0;
  m_unflushed_negative_records = 0;

  for (auto& [voxel, pending] : m_pendingRecords)
  {
    if (!pending.positive.empty() && !pending.negative.empty())
    {
      ++m_final_flush_batches;
      m_final_flush_positive_records += pending.positive.size();
      m_final_flush_negative_records += pending.negative.size();
      processBatch(voxel);
      pending.positive.clear();
      pending.negative.clear();
      continue;
    }

    m_unflushed_positive_records += pending.positive.size();
    m_unflushed_negative_records += pending.negative.size();
  }
}

void PHCPMTpcCalibration::processBatch(const VoxelId& voxel)
{
  auto iter = m_pendingRecords.find(voxel);
  if (iter == m_pendingRecords.end())
  {
    return;
  }

  auto& pending = iter->second;
  const int cellIndex = m_correctionContainer.get_cell_index(voxel.iphi, voxel.ir, voxel.iz);
  if (cellIndex < 0)
  {
    return;
  }

  CPMPairOptions options;
  options.solver = CPMPairSolver::Helix;
  options.min_pt = m_minPt;
  options.max_pair_dca = m_maxPairDca;
  options.magnetic_field_z = m_magneticFieldZ;

  const auto positiveRecords = pending.positive.size();
  const auto negativeRecords = pending.negative.size();
  const auto voxelCenter = getVoxelCenter(voxel);

  ++m_batches;
  for (std::size_t ipos = 0; ipos < positiveRecords; ++ipos)
  {
    const auto& positive = pending.positive[ipos];

    for (std::size_t ineg = 0; ineg < negativeRecords; ++ineg)
    {
      const auto& negative = pending.negative[ineg];
      if (sameTrack(positive, negative))
      {
        continue;
      }

      const auto result = computeCPMPair(
          makeCPMPairInput(positive),
          makeCPMPairInput(negative),
          voxelCenter,
          options);
      if (result.status == CPMPairStatus::PtRejected ||
          result.status == CPMPairStatus::InvalidWeight)
      {
        continue;
      }

      ++m_candidate_pairs;
      if (!result.accepted())
      {
        continue;
      }

      m_correctionContainer.add_sample(
          cellIndex,
          result.delta_r,
          result.delta_phi,
          result.delta_rphi,
          result.delta_z,
          result.dca,
          result.pair_weight);
      ++m_accepted_pairs;
    }
  }
}

bool PHCPMTpcCalibration::checkTrack(const SvtxTrack* track) const
{
  if (!track)
  {
    return false;
  }

  if (m_requireCrossing && track->get_crossing() != 0)
  {
    return false;
  }

  if (track->get_pt() < m_minPt)
  {
    return false;
  }

  if (m_requireTPOT &&
      countTrackClusters(track, TrkrDefs::micromegasId) == 0 &&
      countTrackStates(track, TrkrDefs::micromegasId) == 0)
  {
    return false;
  }

  return true;
}

bool PHCPMTpcCalibration::checkState(const SvtxTrackState* state) const
{
  if (!state)
  {
    return false;
  }

  const auto cluskey = state->get_cluskey();
  if (cluskey == InvalidClusterKey)
  {
    return false;
  }

  return TrkrDefs::getTrkrId(cluskey) == TrkrDefs::tpcId;
}

bool PHCPMTpcCalibration::getVoxelId(const Vector3& position, VoxelId& voxel) const
{
  double phi = std::atan2(position.y, position.x);
  if (phi < 0.0)
  {
    phi += m_phiMax;
  }

  if (phi < m_phiMin || phi >= m_phiMax)
  {
    return false;
  }

  const double radius = std::sqrt(position.x * position.x + position.y * position.y);
  if (radius < m_rMin || radius >= m_rMax)
  {
    return false;
  }

  if (position.z < m_zMin || position.z >= m_zMax)
  {
    return false;
  }

  voxel.iphi = static_cast<int>(m_phiBins * (phi - m_phiMin) / (m_phiMax - m_phiMin));
  voxel.ir = static_cast<int>(m_rBins * (radius - m_rMin) / (m_rMax - m_rMin));
  voxel.iz = static_cast<int>(m_zBins * (position.z - m_zMin) / (m_zMax - m_zMin));

  return voxel.valid();
}

int PHCPMTpcCalibration::writeOutput() const
{
  auto output = std::unique_ptr<TFile>(TFile::Open(m_outputfile.c_str(), "RECREATE"));
  if (!output || output->IsZombie())
  {
    std::cout << "PHCPMTpcCalibration::writeOutput - failed to open "
              << m_outputfile << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  std::unique_ptr<TTree> records;
  std::unique_ptr<CPMRecordTreeFields> fields;
  if (m_writeRecords)
  {
    records = std::make_unique<TTree>("cpm_records", "CPM ACTS-ready voxel records");
    fields = std::make_unique<CPMRecordTreeFields>();
    book_cpm_record_tree(*records, *fields);

    for (const auto& [voxel, voxelRecords] : m_voxelContainer)
    {
      (void) voxel;
      for (const auto& record : voxelRecords)
      {
        fields->copy_from(record);
        records->Fill();
      }
    }
  }

  TTree metadata("cpm_metadata", "CPM Job A metadata");
  int phiBins = m_phiBins;
  int rBins = m_rBins;
  int zBins = m_zBins;
  double rMin = m_rMin;
  double rMax = m_rMax;
  double zMin = m_zMin;
  double zMax = m_zMax;
  double minPt = m_minPt;
  unsigned long long totalTracks = m_total_tracks;
  unsigned long long acceptedTracks = m_accepted_tracks;
  unsigned long long totalStates = m_total_states;
  unsigned long long acceptedStates = m_accepted_states;
  unsigned long long candidatePairs = m_candidate_pairs;
  unsigned long long acceptedPairs = m_accepted_pairs;
  unsigned long long batches = m_batches;
  unsigned long long finalFlushBatches = m_final_flush_batches;
  unsigned long long finalFlushPositiveRecords = m_final_flush_positive_records;
  unsigned long long finalFlushNegativeRecords = m_final_flush_negative_records;
  unsigned long long unflushedPositiveRecords = m_unflushed_positive_records;
  unsigned long long unflushedNegativeRecords = m_unflushed_negative_records;
  unsigned int runningBatchSize = m_runningBatchSize;
  double maxPairDca = m_maxPairDca;
  double magneticFieldZ = m_magneticFieldZ;
  bool writeRecords = m_writeRecords;

  metadata.Branch("phi_bins", &phiBins);
  metadata.Branch("r_bins", &rBins);
  metadata.Branch("z_bins", &zBins);
  metadata.Branch("r_min", &rMin);
  metadata.Branch("r_max", &rMax);
  metadata.Branch("z_min", &zMin);
  metadata.Branch("z_max", &zMax);
  metadata.Branch("min_pt", &minPt);
  metadata.Branch("total_tracks", &totalTracks);
  metadata.Branch("accepted_tracks", &acceptedTracks);
  metadata.Branch("total_states", &totalStates);
  metadata.Branch("accepted_states", &acceptedStates);
  metadata.Branch("candidate_pairs", &candidatePairs);
  metadata.Branch("accepted_pairs", &acceptedPairs);
  metadata.Branch("batches", &batches);
  metadata.Branch("final_flush_batches", &finalFlushBatches);
  metadata.Branch("final_flush_positive_records", &finalFlushPositiveRecords);
  metadata.Branch("final_flush_negative_records", &finalFlushNegativeRecords);
  metadata.Branch("unflushed_positive_records", &unflushedPositiveRecords);
  metadata.Branch("unflushed_negative_records", &unflushedNegativeRecords);
  metadata.Branch("running_batch_size", &runningBatchSize);
  metadata.Branch("max_pair_dca", &maxPairDca);
  metadata.Branch("magnetic_field_z", &magneticFieldZ);
  metadata.Branch("write_records", &writeRecords);
  metadata.Fill();

  output->cd();
  if (records)
  {
    records->Write();
  }
  m_correctionContainer.Write(m_outputObjectName.c_str());
  metadata.Write();
  output->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

TrackStateRecord PHCPMTpcCalibration::makeRecord(
    const unsigned int trackKey,
    const SvtxTrack* track,
    const SvtxTrackState* state,
    const TrkrCluster* cluster,
    const Vector3& clusterPosition,
    const VoxelId& voxel) const
{
  TrackStateRecord record;

  const auto cluskey = state->get_cluskey();
  const auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);

  record.event_ref = makeEventReference();
  record.track_ref.track_id = trackKey;
  record.track_ref.state_cluskey = cluskey;
  record.track_ref.track_map_name = m_trackmapname;
  record.track = makeTrackSummary(track);

  record.cluster_ref.cluskey = cluskey;
  record.cluster_ref.hitsetkey = hitsetkey;
  record.cluster_ref.subsurfkey = cluster->getSubSurfKey();
  record.cluster_ref.layer = TrkrDefs::getLayer(cluskey);
  record.cluster_ref.side = TpcDefs::getSide(cluskey);

  record.voxel = voxel;
  record.cluster.corrected_position = clusterPosition;
  record.cluster.voxel_center = getVoxelCenter(voxel);
  record.cluster.cluster_minus_voxel_center = clusterPosition - record.cluster.voxel_center;

  record.state.pathlength = state->get_pathlength();
  record.state.local_x = state->get_localX();
  record.state.local_y = state->get_localY();
  record.state.position = {state->get_x(), state->get_y(), state->get_z()};
  record.state.momentum = {state->get_px(), state->get_py(), state->get_pz()};
  record.state.covariance = copyCovariance(state);

  record.selection.has_crossing = track->get_crossing() == 0;
  record.selection.passes_tpot =
      !m_requireTPOT ||
      countTrackClusters(track, TrkrDefs::micromegasId) > 0 ||
      countTrackStates(track, TrkrDefs::micromegasId) > 0;
  record.selection.passes_track_quality = true;
  record.selection.passes_geometry = record.cluster_ref.subsurfkey != InvalidSubSurfKey;

  return record;
}

EventReference PHCPMTpcCalibration::makeEventReference() const
{
  EventReference out;
  out.cluster_source = m_cluster_source;
  out.track_source = m_track_source;
  out.run = m_run;
  out.segment = m_segment;
  out.stream_event_ordinal = m_event;

  if (m_syncObject)
  {
    out.sync_event = m_syncObject->EventNumber();
  }

  if (m_eventHeader)
  {
    if (out.run < 0)
    {
      out.run = m_eventHeader->get_RunNumber();
    }
    out.event_sequence = m_eventHeader->get_EvtSequence();
  }

  return out;
}

TrackSummary PHCPMTpcCalibration::makeTrackSummary(const SvtxTrack* track) const
{
  TrackSummary out;
  out.charge = track->get_charge();
  out.pt = track->get_pt();
  out.quality = track->get_quality();
  out.n_mvtx = countTrackClusters(track, TrkrDefs::mvtxId);
  out.n_intt = countTrackClusters(track, TrkrDefs::inttId);
  out.n_tpc = countTrackClusters(track, TrkrDefs::tpcId);
  out.n_tpot = countTrackClusters(track, TrkrDefs::micromegasId);
  out.n_mvtx_states = countTrackStates(track, TrkrDefs::mvtxId);
  out.n_intt_states = countTrackStates(track, TrkrDefs::inttId);
  out.n_tpc_states = countTrackStates(track, TrkrDefs::tpcId);
  out.n_tpot_states = countTrackStates(track, TrkrDefs::micromegasId);
  return out;
}

Vector3 PHCPMTpcCalibration::getVoxelCenter(const VoxelId& voxel) const
{
  const double phi = m_phiMin + (voxel.iphi + 0.5) * (m_phiMax - m_phiMin) / m_phiBins;
  const double radius = m_rMin + (voxel.ir + 0.5) * (m_rMax - m_rMin) / m_rBins;
  const double z = m_zMin + (voxel.iz + 0.5) * (m_zMax - m_zMin) / m_zBins;

  return {radius * std::cos(phi), radius * std::sin(phi), z};
}

Matrix6 PHCPMTpcCalibration::copyCovariance(const SvtxTrackState* state)
{
  Matrix6 out{};
  for (int row = 0; row < 6; ++row)
  {
    for (int col = 0; col < 6; ++col)
    {
      out[row * 6 + col] = state->get_error(row, col);
    }
  }
  return out;
}

unsigned int PHCPMTpcCalibration::countTrackStates(const SvtxTrack* track, const unsigned int trkrId)
{
  unsigned int out = 0;
  for (auto iter = track->begin_states(); iter != track->end_states(); ++iter)
  {
    const auto* state = iter->second;
    if (state && TrkrDefs::getTrkrId(state->get_cluskey()) == trkrId)
    {
      ++out;
    }
  }
  return out;
}

unsigned int PHCPMTpcCalibration::countTrackClusters(const SvtxTrack* track, const unsigned int trkrId)
{
  unsigned int out = 0;
  for (auto iter = track->begin_cluster_keys(); iter != track->end_cluster_keys(); ++iter)
  {
    if (TrkrDefs::getTrkrId(*iter) == trkrId)
    {
      ++out;
    }
  }
  return out;
}

bool PHCPMTpcCalibration::sameTrack(
    const TrackStateRecord& lhs,
    const TrackStateRecord& rhs)
{
  return lhs.event_ref.cluster_source == rhs.event_ref.cluster_source &&
         lhs.event_ref.track_source == rhs.event_ref.track_source &&
         lhs.event_ref.run == rhs.event_ref.run &&
         lhs.event_ref.segment == rhs.event_ref.segment &&
         lhs.event_ref.sync_event == rhs.event_ref.sync_event &&
         lhs.event_ref.event_sequence == rhs.event_ref.event_sequence &&
         lhs.event_ref.stream_event_ordinal == rhs.event_ref.stream_event_ordinal &&
         lhs.track_ref.track_id == rhs.track_ref.track_id;
}

bool PHCPMTpcCalibration::isCloserToVoxelCenter(
    const TrackStateRecord& candidate,
    const TrackStateRecord& current)
{
  const double candidateDistance2 = offsetMagnitude2(candidate);
  const double currentDistance2 = offsetMagnitude2(current);
  if (std::isfinite(candidateDistance2) &&
      std::isfinite(currentDistance2) &&
      candidateDistance2 != currentDistance2)
  {
    return candidateDistance2 < currentDistance2;
  }
  if (std::isfinite(candidateDistance2) != std::isfinite(currentDistance2))
  {
    return std::isfinite(candidateDistance2);
  }
  return candidate.cluster_ref.cluskey < current.cluster_ref.cluskey;
}

double PHCPMTpcCalibration::offsetMagnitude2(const TrackStateRecord& record)
{
  const auto& offset = record.cluster.cluster_minus_voxel_center;
  return offset.x * offset.x + offset.y * offset.y + offset.z * offset.z;
}
