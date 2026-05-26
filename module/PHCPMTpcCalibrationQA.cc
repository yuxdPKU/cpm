#include "PHCPMTpcCalibrationQA.h"

#include <TFile.h>
#include <TTree.h>

#include <cstddef>
#include <limits>

namespace
{
  struct CPMQARecordTreeFields
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

  void book_tree(TTree& tree, CPMQARecordTreeFields& fields)
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

void PHCPMTpcCalibrationQA::write_records(
    TFile& output,
    const VoxelContainer& records,
    const std::string& tree_name)
{
  output.cd();
  CPMQARecordTreeFields fields;
  TTree tree(tree_name.c_str(), "CPM detailed track-state QA records");
  tree.SetDirectory(nullptr);
  book_tree(tree, fields);

  for (const auto& [voxel, voxel_records] : records)
  {
    (void) voxel;
    for (const auto& record : voxel_records)
    {
      fields.copy_from(record);
      tree.Fill();
    }
  }

  tree.Write();
}
