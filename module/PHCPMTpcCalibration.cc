#include "PHCPMTpcCalibration.h"

#include "PHCPMTpcCalibrationQA.h"

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

    unsigned long long cluskey = InvalidClusterKey;
    unsigned long long hitsetkey = InvalidHitSetKey;
    unsigned int subsurfkey = InvalidSubSurfKey;

    int iphi = -1;
    int ir = -1;
    int iz = -1;

    double voxel_x = std::numeric_limits<double>::quiet_NaN();
    double voxel_y = std::numeric_limits<double>::quiet_NaN();
    double voxel_z = std::numeric_limits<double>::quiet_NaN();
    double offset_x = std::numeric_limits<double>::quiet_NaN();
    double offset_y = std::numeric_limits<double>::quiet_NaN();
    double offset_z = std::numeric_limits<double>::quiet_NaN();

    double state_x = std::numeric_limits<double>::quiet_NaN();
    double state_y = std::numeric_limits<double>::quiet_NaN();
    double state_z = std::numeric_limits<double>::quiet_NaN();
    double state_px = std::numeric_limits<double>::quiet_NaN();
    double state_py = std::numeric_limits<double>::quiet_NaN();
    double state_pz = std::numeric_limits<double>::quiet_NaN();

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

      cluskey = record.cluster_ref.cluskey;
      hitsetkey = record.cluster_ref.hitsetkey;
      subsurfkey = record.cluster_ref.subsurfkey;

      iphi = record.voxel.iphi;
      ir = record.voxel.ir;
      iz = record.voxel.iz;

      voxel_x = record.cluster.voxel_center.x;
      voxel_y = record.cluster.voxel_center.y;
      voxel_z = record.cluster.voxel_center.z;
      offset_x = record.cluster.cluster_minus_voxel_center.x;
      offset_y = record.cluster.cluster_minus_voxel_center.y;
      offset_z = record.cluster.cluster_minus_voxel_center.z;

      state_x = record.state.position.x;
      state_y = record.state.position.y;
      state_z = record.state.position.z;
      state_px = record.state.momentum.x;
      state_py = record.state.momentum.y;
      state_pz = record.state.momentum.z;
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

    tree.Branch("cluskey", &fields.cluskey);
    tree.Branch("hitsetkey", &fields.hitsetkey);
    tree.Branch("subsurfkey", &fields.subsurfkey);

    tree.Branch("iphi", &fields.iphi);
    tree.Branch("ir", &fields.ir);
    tree.Branch("iz", &fields.iz);

    tree.Branch("voxel_x", &fields.voxel_x);
    tree.Branch("voxel_y", &fields.voxel_y);
    tree.Branch("voxel_z", &fields.voxel_z);
    tree.Branch("offset_x", &fields.offset_x);
    tree.Branch("offset_y", &fields.offset_y);
    tree.Branch("offset_z", &fields.offset_z);

    tree.Branch("state_x", &fields.state_x);
    tree.Branch("state_y", &fields.state_y);
    tree.Branch("state_z", &fields.state_z);
    tree.Branch("state_px", &fields.state_px);
    tree.Branch("state_py", &fields.state_py);
    tree.Branch("state_pz", &fields.state_pz);
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
            << " trackmap: " << m_trackmapname
            << " grid: (" << m_phiBins << ", " << m_rBins << ", " << m_zBins << ")"
            << " write records: " << m_writeRecords
            << " write QA records: " << m_writeQARecords
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
  std::cout << "PHCPMTpcCalibration::End"
            << " records: " << m_voxelContainer.record_count()
            << " voxels: " << m_voxelContainer.voxel_count()
            << " write_records: " << m_writeRecords
            << " write_qa_records: " << m_writeQARecords
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
  if (m_writeRecords || m_writeQARecords)
  {
    m_voxelContainer.add(std::move(record));
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

  std::cout << "PHCPMTpcCalibration::writeOutput - writing "
            << m_outputfile
            << " records: " << m_voxelContainer.record_count()
            << " qa_records: " << m_writeQARecords
            << std::endl;

  output->cd();
  if (m_writeRecords)
  {
    CPMRecordTreeFields fields;
    TTree records("cpm_records", "CPM ACTS-ready voxel records");
    records.SetDirectory(nullptr);
    book_cpm_record_tree(records, fields);

    for (const auto& [voxel, voxelRecords] : m_voxelContainer)
    {
      (void) voxel;
      for (const auto& record : voxelRecords)
      {
        fields.copy_from(record);
        records.Fill();
      }
    }

    output->cd();
    records.Write();
  }

  TTree metadata("cpm_metadata", "CPM Job A metadata");
  metadata.SetDirectory(nullptr);
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
  bool writeRecords = m_writeRecords;
  bool writeQARecords = m_writeQARecords;

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
  metadata.Branch("write_records", &writeRecords);
  metadata.Branch("write_qa_records", &writeQARecords);
  metadata.Fill();

  output->cd();
  if (m_writeQARecords)
  {
    PHCPMTpcCalibrationQA::write_records(*output, m_voxelContainer);
  }
  metadata.Write();
  output->Close();

  std::cout << "PHCPMTpcCalibration::writeOutput - done "
            << m_outputfile << std::endl;

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
