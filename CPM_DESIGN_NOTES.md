# Crossing Point Method Design Notes

## Project Goal

Build a Crossing Point Method (CPM) calibration prototype for sPHENIX TPC space-charge distortion correction.

The intended workflow is to read sPHENIX DST files, collect distorted TPC cluster positions into a voxel map, use matched silicon-TPOT tracks as approximately undistorted reference trajectories, estimate the undistorted voxel position from track-pair crossing points, and write a per-voxel correction map.

## Method Summary

- TPC clusters are binned in distorted TPC space.
- Clusters in the same voxel are assumed to share one local correction vector, but they are not treated as the same physical point.
- For each cluster-track entry, the intra-voxel offset `cluster - voxel center` must be applied before comparing track-pair crossings.
- Tracks should be grouped by charge sign before pairwise crossing calculations.
- The voxel correction is `r_undistorted_voxel - r_distorted_voxel`.
- A weighted average over accepted pairwise crossing estimates should define the voxel correction.

## Inspected Files

- `/Users/yuxd/Desktop/CPM/CPM_METHOD.md`
- `/Users/yuxd/Desktop/CPM/TASK_CPM.md`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpccalib/PHTpcResiduals.h`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpccalib/PHTpcResiduals.cc`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpc/TpcGlobalPositionWrapper.h`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpc/TpcGlobalPositionWrapper.cc`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/trackreco/MakeSourceLinks.cc`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/trackbase/TrackFitUtils.h`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/trackbase/TrackFitUtils.cc`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/trackreco/PHSimpleVertexFinder.cc`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/TrackingDiagnostics/KshortReconstruction.cc`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpccalib/TpcSpaceChargeMatrixContainer.h`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpccalib/TpcSpaceChargeMatrixContainerv2.h`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpccalib/TpcSpaceChargeMatrixContainerv2.cc`
- `/Users/yuxd/Desktop/CPM/macros/TrackingProduction/Fun4All_FullReconstruction.C`
- `/Users/yuxd/Desktop/CPM/macros/TrackingProduction/run3pp/Fun4All_raw_hit_KFP.C`
- `/Users/yuxd/Desktop/CPM/macros/TrackingProduction/run3auau/Fun4All_PRDFReconstruction.C`
- `/Users/yuxd/Desktop/CPM/ProdFlow/short/run3pp/tracking_code/Fun4All_JobA.C`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/trackreco/PHActsTrkFitter.cc`
- `/Users/yuxd/Desktop/CPM/macros/common/Trkr_RecoInit.C`
- `/Users/yuxd/Desktop/CPM/macros/common/Trkr_Reco.C`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpc/TpcLoadDistortionCorrection.h`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpc/TpcLoadDistortionCorrection.cc`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpccalib/TpcSpaceChargeMatrixInversion.h`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpccalib/TpcSpaceChargeMatrixInversion.cc`
- `/Users/yuxd/Desktop/CPM/macros/calibrations/tpc/jobB/DistortionCorrectionMatrixInversion.C`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpccalib/TpcLaminationFitting.h`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/tpccalib/TpcLaminationFitting.cc`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/trackreco/ActsPropagator.h`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/trackreco/ActsPropagator.cc`
- `/Users/yuxd/Desktop/CPM/coresoftware/offline/packages/trackbase_historic/SvtxTrackState.h`
- `/Users/yuxd/Desktop/CPM/acts/Core/include/Acts/Surfaces/LineSurface.hpp`
- `/Users/yuxd/Desktop/CPM/currentworkflow/Fun4All_TrackAnalysis.C`
- `/Users/yuxd/Desktop/CPM/currentworkflow/DistortionCorrectionMatrixInversion.C`
- External reference macro:
  `https://github.com/yuxdPKU/TPCdistortion/blob/main/Si_TPOT_fit/staticCorrOn_scale1_CDBavgCorrOn_oldTPCalignment_newTPOTalignment/macro/Fun4All_TrackAnalysis.C`
- External reference matrix-inversion macro:
  `https://github.com/yuxdPKU/TPCdistortion/blob/main/Si_TPOT_fit/staticCorrOn_scale1_CDBavgCorrOn_oldTPCalignment_newTPOTalignment/jobB/DistortionCorrectionMatrixInversion.C`

The external GitHub files were later copied into `/Users/yuxd/Desktop/CPM/currentworkflow`, so the CPM plan should now be based on these local copies rather than the earlier partial GitHub fetch.

## Relevant Existing Components

### `PHTpcResiduals`

`PHTpcResiduals` is a `SubsysReco` module used as the closest existing model for the CPM prototype.

Key lifecycle behavior:

- `Init()` prints configuration and resets counters.
- `InitRun()` loads nodes, creates output nodes, and initializes TPC z limits from `ActsGeometry`.
- `process_event()` calls `processTracks()` once per event.
- `End()` writes a `TpcSpaceChargeMatrixContainer` object to a ROOT file.

Important inputs:

- Track map node: default `SvtxSiliconMMTrackMap`, configurable through `setTrackMapName`.
- Cluster node: `TRKR_CLUSTER`.
- Geometry node: `ActsGeometry`.
- TPC global-position helper: `TpcGlobalPositionWrapper::loadNodes(topNode)`.

Important behavior:

- Track cluster keys are collected from `track->get_silicon_seed()` and `track->get_tpc_seed()`.
- Track states are accessed through `track->begin_states()` and matched to cluster keys with `SvtxTrackState::get_cluskey()`.
- TPC cluster positions are obtained with `m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluskey, cluster, crossing)`.
- The existing residual method extrapolates the fitted state linearly to the TPC cluster radius and accumulates matrix terms, rather than storing track-pair information.
- Voxel indexing uses `(phi, r, z)` bins with hard-coded TPC radial limits of 20-78 cm and z limits from geometry.

### `TpcGlobalPositionWrapper`

`TpcGlobalPositionWrapper` converts cluster coordinates to global coordinates and, for TPC clusters, applies:

- crossing z correction through `TpcClusterZCrossingCorrection::correctZ`;
- module-edge distortion correction if the node is loaded and enabled;
- static distortion correction if the node is loaded and enabled;
- average distortion correction if the node is loaded and enabled;
- fluctuation distortion correction if the node is loaded and enabled.

The wrapper loads these correction nodes when present:

- `TpcDistortionCorrectionContainerModuleEdge`
- `TpcDistortionCorrectionContainerStatic`
- `TpcDistortionCorrectionContainerAverage`
- `TpcDistortionCorrectionContainerFluctuation`

The correction actually applied by `getGlobalPositionDistortionCorrected` depends on which distortion nodes are loaded and which wrapper-level enable flags are left on. It should not be described as unconditionally "fully corrected".

For the current CPM task, the target is the average distortion correction. `PHTpcResiduals` is the existing average-correction extraction scheme, and CPM is a parallel replacement candidate for that scheme. In average-correction studies, module-edge and static corrections are expected to be applied as the baseline/pre-corrections. The new average correction should be understood relative to those effects, not as a replacement for them.

### `TrackFitUtils` and Existing Track-Pair DCA Code

`TrackFitUtils` provides circle fits, helix-to-point PCA utilities, helix tangent approximations, and surface/line intersection helpers. The helix PCA code explicitly uses an approximate PCA in x-y to the fitted circle and z from a z-vs-r line, followed by a local straight-line approximation.

`PHSimpleVertexFinder` and `KshortReconstruction` contain two-track DCA/PCA calculations, but these are local straight-line calculations using track positions and momentum directions. They may be useful for first diagnostics, but they do not satisfy the desired final CPM crossing solver by themselves.

### `TpcSpaceChargeMatrixContainer`

`TpcSpaceChargeMatrixContainer` stores accumulated matrix terms for existing space-charge residual inversion. Version 2 stores per-cell entries plus full and reduced left-hand/right-hand matrix arrays for `(drphi, dz, dr)` style fits.

This matrix container is the direct output of `PHTpcResiduals`, not the final average correction object.

### `TpcSpaceChargeMatrixInversion`

`TpcSpaceChargeMatrixInversion` is the post-processing step for the current `PHTpcResiduals` average-correction workflow.

Current chain:

- Job A produces many `TpcSpaceChargeMatrices*.root` files containing `TpcSpaceChargeMatrixContainer`.
- `DistortionCorrectionMatrixInversion.C` loads those files into `TpcSpaceChargeMatrixInversion`.
- `calculate_distortion_corrections()` solves the per-voxel matrices and creates 3D histograms named like `hDistortionP_rec`, `hDistortionR_rec`, and `hDistortionZ_rec`.
- These histograms are split into `negz` and `posz`, guard bins are added, and the result is stored in a `TpcDistortionCorrectionContainer`-style object using `hIntDistortionP`, `hIntDistortionR`, `hIntDistortionZ`, and `hentries`.
- `extrapolate_distortion_corrections()` can use central membrane distortion corrections to extrapolate from TPOT acceptance to the rest of the TPC acceptance.
- `save_distortion_corrections()` writes histograms in the same format consumed by `TpcLoadDistortionCorrection` as an average correction file.

For CPM, the natural replacement point is likely after the track/cluster accumulation: CPM can produce average-correction histograms directly, bypassing the matrix-container inversion step, but it should still write the same histogram names and coordinate conventions expected by `TpcDistortionCorrectionContainer`.

### `TpcLaminationFitting`

`TpcLaminationFitting` is another producer of `TpcDistortionCorrectionContainerAverage`, based on laser/central membrane lamination information. It reads module-edge and static correction containers as inputs and creates average-style 2D correction histograms. This reinforces the convention that average correction is extracted on top of module-edge and static corrections.

### Production and Fitting Macros

The available tracking production macros live under `/Users/yuxd/Desktop/CPM/macros/TrackingProduction`, not under `coresoftware/macro/TrackingProduction`.

Relevant macro behavior:

- `Fun4All_FullReconstruction.C` and run-3 AuAu PRDF reconstruction configure `PHActsTrkFitter`.
- In space-charge calibration mode, the macros call `actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE)`.
- The same calibration-mode block registers `PHTpcResiduals`, sets `setMinPt(0.2)`, and sets grid dimensions to `(36, 48, 80)`.
- `run3pp/Fun4All_raw_hit_KFP.C` enables module-edge and static corrections by default. The average correction lines are present but commented out, with a note that average corrections are designed to be used only if static corrections are also applied.
- `TrackingInit()` registers `TpcLoadDistortionCorrection` when any of module-edge, static, or average correction flags are enabled.
- `TpcLoadDistortionCorrection` maps correction types to node names:
  `Static -> TpcDistortionCorrectionContainerStatic`,
  `Average -> TpcDistortionCorrectionContainerAverage`,
  `Fluctuation -> TpcDistortionCorrectionContainerFluctuation`,
  `ModuleEdge -> TpcDistortionCorrectionContainerModuleEdge`.
- The short run-3 p+p `ProdFlow` JobA writes seed containers and builds separate `SiliconSvtxTrackMap` and `TpcSvtxTrackMap`, but that short workflow is not the direct `SvtxSiliconMMTrackMap` calibration path.
- Current inspected output-node lists do not persist `SvtxSiliconMMTrackMap`; they generally write `SvtxTrackMap`, seed containers, and `TRKR_CLUSTER`. CPM should therefore run in the same Fun4All job after `PHActsTrkFitter` creates `SvtxSiliconMMTrackMap`, unless the macro is explicitly modified to write that node.

### External PHTpcResiduals Workflow References

The user's current `TPCdistortion` reference workflow should be treated as the practical production template for CPM, not just as an example macro.

Observed Job A structure in `currentworkflow/Fun4All_TrackAnalysis.C`:

- `G4TRACKING::SC_CALIBMODE = true`, `G4TRACKING::SC_USE_MICROMEGAS = true`, and `TRACKING::pp_mode = true`.
- Module-edge, static, and average correction loading are all enabled globally. The average correction file is taken from CDB through `TPC_LAMINATION_FIT_CORRECTION`.
- The first `PHActsTrkFitter` is a normal fit with `fitSiliconMMs(false)` and `setUseMicromegas(true)`.
- `PHSimpleVertexFinder`, `PHActsVertexPropagator`, and `PHTrackPruner` then select/prune tracks. The pruned seed map is named `PrunedSvtxTrackSeedContainer`.
- A second `PHActsTrkFitter` performs the Si-TPOT calibration fit with `fitSiliconMMs(G4TRACKING::SC_CALIBMODE)`, `setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS)`, and `set_svtx_seed_map_name("PrunedSvtxTrackSeedContainer")`.
- `PHTpcResiduals` is registered after this second fitter, so it uses the Si-TPOT reference fit and the TPC states propagated from that fit.
- `PHTpcResiduals` explicitly calls `disableAverageCorr()`. This means that even though the average correction is enabled globally for earlier reconstruction/matching, the residual extraction itself fills TPC positions with the average component disabled in the wrapper.
- The residual cuts are tuned for this workflow: `setMinPt(0.5)`, `requireCrossing(false)`, `requireCM(true)`, `setPCAzcut(10)`, `setEtacut(0.25)`, `setMaxTrackAlpha(0.6)`, `setMaxTrackBeta(1.5)`, `setMaxTrackResidualDrphi(2)`, `setMaxTrackResidualDz(5)`, `setMinRPhiErr(0.005)`, and `setMinZErr(0.01)`.
- The active reconstructed distortion grid is `(phi, r, z) = (36, 16, 80)`. The earlier `setGridDimensions(48)` call is superseded by `setGridDimensions(36, 16, 80)`.

Observed Job B structure in `currentworkflow/DistortionCorrectionMatrixInversion.C`:

- TPOT acceptance is configured with explicit hard-coded phi and theta ranges for central, east, and west TPOT regions, rather than using `load_tpot_geometry`.
- Input files are read from `../root/Reconstructed/<run>/clusters_seeds_<run>-*.root_PhTpcResiduals.root`.
- The default run argument is `53285`.
- A central-membrane input file path is declared, but the CM loading calls are currently commented for the 3D inversions.
- The 1D inversion block is present but commented out.
- The 2D radius-z chain is active in the matrix-inversion workflow as a test branch. It reads `TpcSpaceChargeMatrixContainer_2D_radius_z`, sets `min_cluster_count = 50`, calls `calculate_distortion_corrections()`, calls `extrapolate_distortion_corrections()`, and writes `Rootfiles/Distortions_2D_mm_<run>_rz.root`. This is not part of the CPM target design.
- The full 3D chain is active. It reads `TpcSpaceChargeMatrixContainer`, sets `min_cluster_count = 1000`, calls `calculate_distortion_corrections()`, and writes `Rootfiles/Distortions_full_mm_<run>.root`.
- Two reduced 3D chains are active: `ReducedInversion_phi` writes `Rootfiles/Distortions_full_mm_<run>_phi.root`, and `ReducedInversion_z` writes `Rootfiles/Distortions_full_mm_<run>_z.root`. The variable name for the z reduced inversion is `spaceChargeMatrixInversion_r`, but the mode and output name indicate the z-reduced path.

Design consequences for CPM Job A:

- CPM Job A should mirror the two-fit production structure: first build/prune good full tracks, then run the second Si-TPOT-only `PHActsTrkFitter`, then register CPM after the second fitter when `SvtxSiliconMMTrackMap` is available.
- The baseline distortion configuration should follow the average-correction production convention. Module-edge and static corrections are loaded as pre-corrections. The average correction may be enabled globally for matching and preliminary reconstruction, but CPM's cluster-position wrapper should disable the average component when producing a new average correction, analogous to `PHTpcResiduals::disableAverageCorr()`.
- The CPM module should expose macro-level setters analogous to `PHTpcResiduals`: output file, track map name, minimum track pT, grid dimensions, and optional enable/disable controls for the wrapper-level correction components.
- For comparison with the user's current production chain, the first CPM grid should use `(phi, r, z) = (36, 16, 80)`. The broader `(36, 48, 80)` grid remains useful as the default from other macros, so the grid must be a macro-level setting.
- CPM should reuse the same first-pass selection logic where possible: pruned seed map `PrunedSvtxTrackSeedContainer`, `pT > 0.5`, at least one TPOT cluster/state, and the same crossing/CM requirements if the required node information is available.
- Job A output should be segment-friendly. It can either write a compact CPM intermediate container or write per-segment correction/QA histograms directly, but it should not require all input DSTs to be processed in one Fun4All job.

Design consequences for CPM Job B:

- CPM should keep a second-stage macro analogous to `DistortionCorrectionMatrixInversion.C`, even if CPM does not solve the same residual matrix equations.
- Job B should aggregate outputs from many Job A segments, apply statistics/quality cuts, create final per-voxel correction histograms, split negative and positive z, add guard bins if needed, and write the average-correction histogram names consumed by `TpcDistortionCorrectionContainer`.
- Job B should own TPOT-acceptance and central-membrane extrapolation policy. The current workflow uses explicit hard-coded TPOT phi/theta ranges; CPM should support explicit override ranges as the production-compatible mode, with geometry-derived TPOT ranges as an optional convenience mode.
- CPM Job B should produce only the 3D correction-map version. The matrix-inversion workflow's 2D r-z branch is a test of that method and should not be carried into CPM. Reduced phi/z-like diagnostic variants can be considered later only if they help QA the 3D CPM result, but they are not part of the v1 target.
- The final file should be loadable through `G4TPC::average_correction_filename` and `TpcLoadDistortionCorrection` as an average correction file, rather than requiring a CPM-specific consumer in reconstruction.

This leads to a two-stage CPM design:

1. CPM Job A, `SubsysReco` stage: collect corrected TPC cluster positions, silicon-TPOT reference states, local voxel offsets, charge sign, and line-PoCA QA quantities from each event.
2. CPM Job B, post-processing stage: merge Job A outputs, compute robust 3D voxel corrections, perform optional TPOT/CM extrapolation in 3D, and write `TpcDistortionCorrectionContainerAverage`-compatible output.

### CPM Intermediate Data Model

Job A should not write only the minimum fields needed for local line-line PoCA. The intermediate file should be ACTS-ready, so that a later ACTS trajectory crossing solver can be implemented in Job B without rerunning the Fun4All reconstruction step.

For each accepted TPC cluster-track entry, Job A should store:

- event identity: source filename, run, segment, Sync event number if available, `EventHeader::get_EvtSequence()` if available, and the local processed-event ordinal;
- track identity: track id, charge, pT, quality, n-MVTX/n-INTT/n-TPC/n-TPOT cluster and state counters if available;
- cluster identity: `cluskey`, `hitsetkey`, layer, side, sector-like indices if convenient, and `subsurfkey`;
- voxel identity: global voxel id plus `(iphi, ir, iz)`;
- cluster geometry: average-disabled corrected cluster global position, voxel center, and `cluster - voxel_center` offset;
- ACTS state copy from the Si-TPOT reference fit at the associated TPC cluster surface: pathlength, localX, localY, global x/y/z, momentum px/py/pz, and the full 6x6 state covariance from `SvtxTrackState::get_error(i,j)`;
- ACTS surface identity/diagnostics: surface geometry id fields if accessible (`volume`, `layer`, `sensitive`, `approach`, `boundary`), plus optional surface center and normal/local axes in global coordinates for debugging;
- selection flags: crossing availability, CM requirement result, TPOT requirement result, and any residual/angle cuts applied.

The `cluskey` alone is not enough for later ACTS surface reconstruction in Job B, because TPC surface lookup uses `hitsetkey + TrkrCluster::getSubSurfKey()`. Therefore Job A must persist `subsurfkey`. With `cluskey`, `hitsetkey`, `subsurfkey`, and the geometry loaded in Job B, CPM can recover the reference surface through `ActsGeometry::maps().getTpcSurface(hitsetkey, subsurfkey)` and rebuild `Acts::BoundTrackParameters` from the stored state position, momentum, charge, and covariance.

This schema deliberately supports two solver levels:

- v1 solver: use the stored global position/momentum as a local line and compute line-line PoCA.
- future ACTS solver: rebuild bound parameters on the stored/recovered surface and run an ACTS-based crossing/minimization algorithm.

If a future task requires a full ACTS refit in Job B, not just ACTS propagation/minimization from the stored Si-TPOT fit states, then the intermediate format will also need original silicon/TPOT measurement/source-link information. That is out of scope for CPM v1, but the v1 schema should not preclude adding a per-track measurement block later.

### Event Reference and Rehydration Strategy

The intermediate format should use a dual-track data model:

- snapshot fields are the primary data used by normal CPM production;
- reference fields make it possible to go back to the original or intermediate DST and recover the corresponding `SvtxTrack`, `SvtxTrackState`, and `TrkrCluster` for validation, debugging, or future full-fit studies.

The reference fields should include at least:

- source DST filename or file id for the cluster DST;
- source DST filename or file id for the seed/track DST or CPM mini-DST;
- run and segment;
- Sync event number, when the `Sync` node is present;
- `EventHeader::get_EvtSequence()`, when the `EventHeader` node is present;
- local processed-event ordinal within the input stream;
- track map name, expected to be `SvtxSiliconMMTrackMap` for the calibration fit;
- track id inside that map;
- state cluster key, expected to match `SvtxTrackState::get_cluskey()`;
- cluster key in `TRKR_CLUSTER`.

`SvtxSiliconMMTrackMap` can be persisted in a CPM mini-DST. The current workflow already contains a commented output block that writes `SvtxSiliconMMTrackMap`; CPM should make this an explicit optional output. In the cluster-DST production mode, the mini-DST should stay small and store `Sync`, `EventHeader`, and `SvtxSiliconMMTrackMap`. `TRKR_CLUSTER` is intentionally not stored in the mini-DST because it is large and can be recovered from the original input cluster DST recorded as `cluster_source`. Seed containers are not required for the local line-line PoCA solver and should remain optional; they are only needed if a later workflow wants to redo a full ACTS refit from source links rather than using the stored track states/snapshots.

For performance, Job B should not randomly read the DST for every voxel entry or every track pair. If DST rehydration is needed, use a separate Job B0 stage:

1. Scan all `CPMVoxelContainer` files and collect the referenced event/object requests.
2. Group and deduplicate requests by `(filename, run, segment, event id)`.
3. Sort requests in file/event order.
4. Sequentially read each source mini-DST once.
5. For each requested event, pull requested tracks/states from `track_source`/`SvtxSiliconMMTrackMap` and requested clusters from `cluster_source`/`TRKR_CLUSTER`.
6. Write an enriched or validation container, then run the actual crossing-point correction in Job B1.

This keeps the normal CPM path fast and reproducible from the snapshot fields, while still preserving a route back to full Fun4All objects. The important rule is that DST rehydration should be event-ordered and deduplicated; the unit of DST readback is an event, not a pair or a single track.

The `Sync` node helps identify and synchronize events. It stores event/run/segment information internally, and `SyncObject::EventNumber()` is publicly accessible. `EventHeader` also provides public run and event-sequence accessors. However, the standard `Fun4AllDstInputManager` is primarily a sequential event-loop and synchronization tool; it does not expose a simple public random-access API by `(run, segment, event)`. Therefore CPM should store enough reference metadata to build its own event index if random lookup is needed, or prefer sorted sequential rehydration from a mini-DST.

### `PHActsTrkFitter` Silicon-MM Fit

When `m_fitSiliconMMs` is enabled:

- `PHActsTrkFitter` creates or reuses the `SvtxSiliconMMTrackMap` node.
- TPC surfaces are skipped in `filterSourceLinks`, so the ACTS fit uses silicon and optional Micromegas information, not TPC measurements.
- If Micromegas use is required, tracks without Micromegas surfaces are skipped.
- After the silicon-MM fit succeeds, the resulting track keeps the TPC seed and silicon seed.
- The fitted track parameters are propagated to each TPC cluster surface from the associated TPC seed, and a `SvtxTrackState` is added for each TPC cluster key.

This means `SvtxSiliconMMTrackMap` appears to provide the exact data relationship CPM needs: an unbiased silicon-TPOT reference fit plus TPC cluster keys/states for distorted-cluster association.

### ACTS Propagation and Crossing Solver Notes

`ActsPropagator` can build `Acts::BoundTrackParameters` from either an `SvtxTrack` or an `SvtxTrackState`. The state-based helper takes:

- `SvtxTrackState* state`;
- track charge;
- a target/reference `Acts::Surface` pointer.

It can then propagate bound parameters to a sPHENIX layer or to a given ACTS surface.

This is useful for CPM because `PHActsTrkFitter`, in silicon-MM calibration mode, already propagates the silicon-MM fit to the TPC cluster surfaces and stores `SvtxTrackState` objects keyed by TPC cluster key.

ACTS `LineSurface` mathematically supports point-of-closest-approach to a line, but it is an abstract/protected-constructor base class in ACTS. The current sPHENIX `ActsPropagator` does not expose a ready-made "track-to-track PoCA" interface. A production ACTS-based solver will likely need either:

- a small custom target surface/helper representing the local line of the second trajectory;
- or an iterative minimization using repeated ACTS propagation to locally defined surfaces.

For an initial CPM prototype, a local line-line PoCA using the two `SvtxTrackState` positions and momenta is a reasonable validation implementation, as long as the crossing solver remains isolated behind a replaceable interface.

Decision: CPM v1 will use local line-line PoCA from `SvtxTrackState` as the first crossing-point solver. This is accepted as a framework-building approximation. The solver must be implemented behind a replaceable interface so that a later ACTS-based trajectory crossing/minimization method can replace it without rewriting voxel filling, accumulation, or output.

## Initial Design Implications

- CPM can initially follow the `SubsysReco` lifecycle pattern from `PHTpcResiduals`.
- `process_event()` should likely collect compact per-cluster/per-track entries into an internal voxel structure.
- Pairwise crossing calculations should run in post-processing rather than inside the Fun4All event loop. Job A should stream events and write the voxel container; Job B should read voxel records, optionally run Job B0 event-ordered rehydration, and then compute crossing-point corrections.
- The first prototype should reuse the existing node access pattern: `SvtxSiliconMMTrackMap`, `TRKR_CLUSTER`, `ActsGeometry`, and `TpcGlobalPositionWrapper`.
- The CPM voxel binning can start from the same `(phi, r, z)` convention as `PHTpcResiduals`, unless later studies show that another convention better matches correction-map output.
- CPM should treat module-edge and static corrections as the default baseline corrections for average-distortion extraction. The open coordinate question is not "raw versus fully corrected" in the abstract, but exactly which corrections are applied before filling the CPM voxel and which residual correction component CPM is meant to output.
- CPM should not need to write `TpcSpaceChargeMatrixContainer` unless we intentionally emulate the old inversion pipeline. Its direct output should be correction-component histograms compatible with `TpcDistortionCorrectionContainer`, plus QA histograms.
- Histogram naming and coordinate conventions should match the current average correction files: `hIntDistortionP_{negz,posz}`, `hIntDistortionR_{negz,posz}`, `hIntDistortionZ_{negz,posz}`, and `hentries_{negz,posz}` after splitting/guard-bin handling.
- Existing line-DCA code can be reused conceptually for validation plots, but the production crossing solver should remain modular and replaceable by an ACTS propagation/minimization implementation.
- The v1 crossing solver is local line-line PoCA from the two reference `SvtxTrackState` positions and momentum directions. Use the midpoint of the two closest points as the pair crossing estimate, and store the pair DCA as a QA/rejection quantity.
- The first production-comparison grid should match the current workflow `PHTpcResiduals` setting: `(phi, r, z) = (36, 16, 80)`. Keep `(36, 48, 80)` available as a macro option because it appears in other calibration macros.
- CPM track entries should store enough snapshot state to allow either local-line PoCA or later ACTS propagation/minimization: cluster key, hitset key, TPC subsurfkey, distorted cluster coordinate, voxel offset, charge, reference state local/global position, momentum, full state covariance, pathlength, track quality, and optional ACTS surface geometry id diagnostics.
- CPM track entries should also store reference metadata sufficient to recover the original Fun4All objects from a persisted mini-DST: source file id, run/segment/event identifiers, track map name, track id, state cluskey, and cluster cluskey.
- CPM v1 should be registered immediately after `PHActsTrkFitter` in SC calibration mode, analogous to `PHTpcResiduals`, rather than relying on a later DST pass.

## Open Questions

- Confirm the exact cluster coordinate convention for CPM voxel filling by tracing `PHTpcResiduals::disableAverageCorr()` through the wrapper in the current workflow. The working assumption is: module-edge + static corrected TPC position, average disabled for the extraction module.
- Should the final correction map be written as `(delta_r, delta_phi, delta_z)`, `(delta_x, delta_y, delta_z)`, or both for diagnostics?
- Which exact event identifier should be treated as canonical for CPM rehydration: Sync event number, `EventHeader::get_EvtSequence()`, local tree entry, or a composite of all three?

## Next Tasks

- Implemented first framework-independent CPM core skeleton in `module`: ACTS-ready records, voxel container, and local line-line PoCA validation test.
- Repository layout now follows the intended CPM package split: `module/` contains source/header/build files, including an autotools skeleton modeled after coresoftware packages, while `macro/` is reserved for Fun4All running macros.
- Drafted the first CPM Job A `SubsysReco` interface around the new record/container model.
- Promoted the CPM Job A macro from a helper skeleton to a complete workflow macro based on `currentworkflow/Fun4All_TrackAnalysis.C`, with `PHTpcResiduals` replaced by `PHCPMTpcCalibration`.
- Added optional CPM mini-DST output to the Job A macro. The current default mini-DST stores `Sync`, `EventHeader`, and `SvtxSiliconMMTrackMap`; `TRKR_CLUSTER` is recovered from the input cluster DST, and pruned seeds are optional because the v1 PoCA solver does not need seed objects.
- Added first Job A flat ROOT persistence: `cpm_records` stores one ACTS-ready state/voxel record per row, and `cpm_metadata` stores grid/counter metadata.
- Drafted the first CPM Job B0 event-index skeleton. It scans `cpm_records`, sorts/deduplicates event references, and writes `cpm_event_requests` plus `cpm_object_requests` for a later sequential mini-DST rehydration pass.
- Added file-list input support for B0 so large Condor productions can pass hundreds or thousands of Job A ROOT outputs without constructing a huge ROOT command line.
- Added a first B0 event-index QA macro that checks object/event consistency before the mini-DST rehydration pass.
- Added the first B1 local line-line PoCA macro. It reads Job A snapshots directly, groups records by 3D voxel, forms different-track pairs, and writes pair-level crossing QA without requiring seeds or `TRKR_CLUSTER` in the mini-DST.
- Draft the second CPM Job B0 rehydration/validation skeleton that reads `cpm_event_requests`, sequentially scans the mini-DST once, and writes an enriched validation output. This should be finalized against the server-side Fun4All input setup so the source-file identity and event synchronization policy are not guessed locally.
- Draft a CPM Job B1 macro skeleton that reads CPM segment or enriched outputs and writes 3D `TpcDistortionCorrectionContainerAverage`-compatible histograms, with explicit TPOT phi/theta range overrides.
- Define whether CPM should include an extrapolation step analogous to `TpcSpaceChargeMatrixInversion::extrapolate_distortion_corrections`, or only produce corrections in the acceptance where crossing-point statistics exist.
- Defer ACTS trajectory minimization until the CPM framework, voxel accumulation, output, and line-PoCA QA are working.
