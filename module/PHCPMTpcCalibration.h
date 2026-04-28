#ifndef CPM_PHCPMTPCCALIBRATION_H
#define CPM_PHCPMTPCCALIBRATION_H

#include "CPMVoxelContainer.h"

#include <fun4all/SubsysReco.h>
#include <tpc/TpcGlobalPositionWrapper.h>

#include <cstdint>
#include <string>

class ActsGeometry;
class EventHeader;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxTrackState;
class SyncObject;
class TrkrCluster;
class TrkrClusterContainer;

namespace cpm
{
  class PHCPMTpcCalibration : public SubsysReco
  {
   public:
    explicit PHCPMTpcCalibration(const std::string& name = "PHCPMTpcCalibration");
    ~PHCPMTpcCalibration() override = default;

    int Init(PHCompositeNode* topNode) override;
    int InitRun(PHCompositeNode* topNode) override;
    int process_event(PHCompositeNode* topNode) override;
    int End(PHCompositeNode* topNode) override;

    void setOutputfile(const std::string& outputfile) { m_outputfile = outputfile; }
    void setClusterSource(const std::string& value) { m_cluster_source = value; }
    void setTrackSource(const std::string& value) { m_track_source = value; }
    void setTrackMapName(const std::string& value) { m_trackmapname = value; }
    void setRunSegment(const int run, const int segment)
    {
      m_run = run;
      m_segment = segment;
    }
    void setMinPt(const double value) { m_minPt = value; }
    void requireCrossing(const bool value = true) { m_requireCrossing = value; }
    void requireTPOT(const bool value = true) { m_requireTPOT = value; }
    void requireCM(const bool value = true) { m_requireCM = value; }

    void setGridDimensions(int phiBins, int rBins, int zBins);

    void disableModuleEdgeCorr()
    {
      m_globalPositionWrapper.set_enable_module_edge_corr(false);
    }

    void disableStaticCorr()
    {
      m_globalPositionWrapper.set_enable_static_corr(false);
    }

    void disableAverageCorr()
    {
      m_globalPositionWrapper.set_enable_average_corr(false);
    }

    void disableFluctuationCorr()
    {
      m_globalPositionWrapper.set_enable_fluctuation_corr(false);
    }

   private:
    int getNodes(PHCompositeNode* topNode);
    int processTracks();
    int writeOutput() const;

    bool checkTrack(const SvtxTrack* track) const;
    bool checkState(const SvtxTrackState* state) const;
    bool getVoxelId(const Vector3& position, VoxelId& voxel) const;

    TrackStateRecord makeRecord(
        unsigned int trackKey,
        const SvtxTrack* track,
        const SvtxTrackState* state,
        const TrkrCluster* cluster,
        const Vector3& clusterPosition,
        const VoxelId& voxel) const;

    EventReference makeEventReference() const;
    TrackSummary makeTrackSummary(const SvtxTrack* track) const;
    SurfaceSnapshot makeSurfaceSnapshot(const TrkrCluster* cluster, ClusterKey cluskey) const;

    Vector3 getVoxelCenter(const VoxelId& voxel) const;
    static Matrix6 copyCovariance(const SvtxTrackState* state);
    static unsigned int countTrackStates(const SvtxTrack* track, unsigned int trkrId);
    static unsigned int countTrackClusters(const SvtxTrack* track, unsigned int trkrId);

    std::string m_trackmapname = "SvtxSiliconMMTrackMap";
    std::string m_outputfile = "CPMVoxelContainer.root";
    std::string m_cluster_source;
    std::string m_track_source;
    int m_run = -1;
    int m_segment = -1;

    SvtxTrackMap* m_trackMap = nullptr;
    ActsGeometry* m_tGeometry = nullptr;
    TrkrClusterContainer* m_clusterContainer = nullptr;
    SyncObject* m_syncObject = nullptr;
    EventHeader* m_eventHeader = nullptr;

    TpcGlobalPositionWrapper m_globalPositionWrapper;
    VoxelContainer m_voxelContainer;

    int m_phiBins = 36;
    int m_rBins = 16;
    int m_zBins = 80;

    static constexpr double m_phiMin = 0.0;
    static constexpr double m_phiMax = 6.28318530717958647692;
    static constexpr double m_rMin = 20.0;
    static constexpr double m_rMax = 78.0;

    double m_zMin = 0.0;
    double m_zMax = 0.0;
    double m_minPt = 0.5;

    bool m_requireCrossing = false;
    bool m_requireTPOT = true;
    bool m_requireCM = true;

    std::uint64_t m_event = 0;

    std::uint64_t m_total_tracks = 0;
    std::uint64_t m_accepted_tracks = 0;
    std::uint64_t m_total_states = 0;
    std::uint64_t m_accepted_states = 0;
  };
}

#endif
