# CPM Procedure For Review

This is a living technical note for the sPHENIX TPC Crossing Point Method
(CPM) average-distortion calibration prototype. It is intended to be updated as
the implementation and physics review evolve.

## Goal

CPM estimates the average TPC space-charge distortion map by comparing distorted
TPC cluster voxel positions with crossing points inferred from silicon-TPOT
reference tracks. The output target is the same average-correction histogram
format consumed by `TpcLoadDistortionCorrection`.

## Coordinate And Correction Convention

- The target correction is the average correction on top of the standard
  module-edge and static pre-corrections.
- During extraction, the TPC global-position wrapper should use module-edge and
  static corrections, but should disable the average correction being measured.
- The stored pair delta convention is:
  `distortion = voxel center - crossing point`.
- This matches the sign convention used by `TpcDistortionCorrection`, where
  `phi_new = phi - dphi`, `r_new = r - dr`, and `z_new = z - dz`.
- The final average-correction histograms use:
  `hIntDistortionR_{negz,posz}`, `hIntDistortionP_{negz,posz}`,
  `hIntDistortionZ_{negz,posz}`, and `hentries_{negz,posz}`.
- For average corrections, `hIntDistortionP` is currently filled as a phi-angle
  correction in radians, consistent with
  `G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = true`.

## Job A: Record Production

Job A runs inside the Fun4All reconstruction/calibration chain, immediately
after the silicon-TPOT calibration fit. For each accepted TPC cluster associated
with a silicon-TPOT reference track, it stores an ACTS-ready record containing:

- event identity: source file, run, segment, Sync event if available,
  EventHeader sequence if available, and stream ordinal;
- track identity: track id, charge, pT, quality, and related counters when
  available;
- cluster identity: cluster key, hitset key, layer-like information, and
  subsurface key when available;
- voxel identity: `(iphi, ir, iz)` in the 3D TPC correction grid;
- cluster geometry: corrected cluster position, voxel center, and
  `cluster - voxel center` offset;
- reference state snapshot: state position, momentum, local parameters, path
  length, and covariance.

The Job A output is segment-friendly: each Condor job writes one compact CPM
ROOT file, and Job B can read either a single file or a file list.

Job A also fills a `CPMCorrectionContainer` during reconstruction. For each
voxel, accepted records are kept in positive-charge and negative-charge FIFO
queues. Whenever both queues reach the configured running threshold, CPM forms
all opposite-charge pairs from the currently queued records in that voxel,
accumulates the resulting distortion samples into the container, and clears the
voxel queues before collecting the next group. This preserves event-order
locality and avoids the large offline track reordering needed by the diagnostic
B1 path.

The mergeable `CPMCorrectionContainer` and `cpm_metadata` are the default Job A
outputs. The flat `cpm_records` tree is a QA/diagnostic product, disabled by
default, and should be enabled with `writeCpmRecords=true` only for B0/B1
rehydration or offline-PoCA studies.

## Job B: Macro Roles

The Job B macros support two routes. The production route follows the
matrix-inversion style: Job A writes additive `CPMCorrectionContainer` objects,
and `CPM_ReconstructAverageCorrection.C` calls
`CPMAverageCorrectionReconstruction` to merge them, calculate average voxel
corrections, and write the final average-correction histograms plus QA objects
in one ROOT file. The diagnostic route starts from the optional `cpm_records`
tree produced with `writeCpmRecords=true` and uses
`CPM_QA_RunOfflineDiagnostics.C` as the single user-facing macro; it still
writes separate B0/B1/B2/B3 stage ROOT files for offline debugging.

Macro responsibilities:

- `CPM_QA_B0_BuildEventIndex.C`: scan one or more Job A `cpm_records` files and
  write `cpm_event_requests` plus `cpm_object_requests` for optional event-wise
  mini-DST rehydration studies.
- `CPM_QA_B0_CheckEventIndex.C`: validate the B0 event/object request index before
  any sequential readback is attempted.
- `CPM_QA_RunOfflineDiagnostics.C`: record-based diagnostic driver. It can run
  optional B0 QA and then writes the B1 PoCA, B2 voxel-correction, and B3
  histogram stage outputs from one macro call.
- `CPM_QA_B1_ComputePoCA.C`: diagnostic offline crossing-point step. It groups
  Job A records by voxel, forms opposite-charge pairs, computes PoCA with the
  helix or line solver, and writes pair-level QA plus batch-level correction
  sums.
- `CPM_QA_B2_AccumulateVoxelCorrections.C`: diagnostic B2 for B1 outputs. It
  accumulates B1 pair rows or B1 batch sums into one `cpm_voxel_corrections`
  row per voxel.
- `CPM_ReconstructAverageCorrection.C`: recommended production macro. It uses
  `CPMAverageCorrectionReconstruction` to merge Job A
  `CPMCorrectionContainer` objects with grid/range guards, calculate weighted
  or plain averages, and write the merged container, `cpm_voxel_corrections`,
  `cpm_metadata`, reconstructed QA histograms, final guarded histograms, and
  summary.
- `CPM_B2_MergeCorrectionContainers.C`: lower-level split-stage container debug
  macro. Its main product is the `cpm_voxel_corrections` TTree consumed by B3.
- `CPM_B3_WriteAverageCorrectionHistograms.C`: convert B2 voxel rows to the
  guarded half-TPC average-correction histograms expected by the correction
  loader.
- `CPM_B3_CheckAverageCorrectionHistograms.C`: verify that the required B3
  histograms and summary outputs exist and have valid dimensions.

The recommended production sequence is:

```text
Job A output with CPMCorrectionContainer
  -> CPM_ReconstructAverageCorrection.C
  -> CPM_B3_CheckAverageCorrectionHistograms.C
```

The diagnostic/offline-PoCA sequence is:

```text
Job A cpm_records, with writeCpmRecords=true
  -> CPM_QA_RunOfflineDiagnostics.C
```

## Job B1: Pair Construction And Crossing Estimates

B1 reads one or more Job A outputs and groups records by voxel.

Current implemented selection:

- require a configurable minimum pT, default `0.5 GeV`;
- require a nonzero charge;
- for each unique track in a voxel, keep exactly one record, chosen as the
  cluster/state closest to the voxel center;
- split selected records by charge sign;
- form only opposite-charge track pairs;
- apply the intra-voxel shift before PoCA:
  `line point = state position - (cluster - voxel center)`.

High-occupancy handling:

- The old first-N cap is replaced by deterministic hash-ordered batching.
- The raw per-voxel safety skip is disabled by default (`--b1-max-records 0`);
  high occupancy is controlled by batch size rather than by dropping the voxel.
- `--b1-max-pair-records N` means at most `N` selected positive records and
  `N` selected negative records per batch.
- The current default is `N = 10`; `N = 0` means one unlimited full-voxel batch.
- B1 writes per-voxel QA including total selected records, selected positive and
  negative records, batch count, batched positive/negative records, theoretical
  opposite-charge pairs, batched opposite-charge pairs, candidate pairs, and
  accepted pairs.
- B1 also writes `cpm_b1_batch_corrections`, one row per accepted batch. Each
  row stores the accepted pair count, total pair weights, effective entries,
  unweighted sums and sum-of-squares, weighted sums and sum-of-squares, and DCA
  QA. This is the persistent running-average input for B2.

Default crossing solver:

- The current default solver is the CPM-local numerical ideal-helix PoCA helper,
  `CPMHelixPoCA`.
- The first implementation uses a configurable uniform `Bz` estimate. The
  default is `1.4 T`. This is an explicit approximation for early algorithm
  development, not a replacement for later field-map or ACTS-based validation.
- The pair crossing estimate is the midpoint of the two closest points.
- The pair DCA is stored as QA and can be cut in B1/B2.
- The previous framework-building local line-line PoCA solver remains available
  through `--b1-crossing-solver line` for comparison.

Target crossing solver:

- The target solver is an ideal two-helix PoCA/crossing method.
- A local coresoftware survey found no ready-made two-helix PoCA utility.
  `TrackFitUtils` provides helix-to-point/surface helpers, while
  `KshortReconstruction` and `PHSimpleVertexFinder` provide line-line two-track
  PCA/DCA.
- CPM now has an initial CPM-local numerical helix PoCA helper,
  `CPMHelixPoCA`, with a narrow interface designed for later upstreaming. It is
  now the default B1 solver, while the local line-line solver remains available
  as a control option.

## Job B2: Voxel Accumulation

In the production route, voxel accumulation now lives inside
`CPMAverageCorrectionReconstruction`. It reads the `CPMCorrectionContainer`
objects written by Job A, checks grid compatibility through the container merge
guard, adds the stored per-voxel sums, and writes one `cpm_voxel_corrections`
row per accepted voxel in the same ROOT file as the final histograms.

The split-stage `CPM_B2_MergeCorrectionContainers.C` macro remains available
when the intermediate voxel TTree should be produced separately. The diagnostic
`CPM_QA_B2_AccumulateVoxelCorrections.C` can still read B1 batch-level
correction sums and accumulate one correction row per voxel. It can also read
pair-level rows for backwards compatibility and detailed QA studies.

Implemented averaging modes:

- weighted mode: use the curvature-proxy pair weight `1/(pt_a * pt_b)`;
- unweighted mode: use unit weights for a simple arithmetic mean.

The voxel TTree stores per-voxel QA including entries, sum of weights,
effective weighted entries, means, and RMS values for `delta_r`, `delta_phi`,
`delta_rphi`, `delta_z`, and DCA. Batch input is pair-equivalent: unweighted
mode combines batch sums with accepted-pair counts, while weighted mode combines
batch weighted sums with total pair weights. It therefore does not average all
batch means equally unless the batch sizes and weights happen to match.

## Job B3: Average-Correction Histograms

B3 converts voxel correction rows into average-correction histograms. In the
production route this happens inside `CPMAverageCorrectionReconstruction`; the
split-stage `CPM_B3_WriteAverageCorrectionHistograms.C` macro is retained for
QA and compatibility. Both paths fill reconstructed QA histograms, split
negative and positive z, apply the existing guard-bin style, and write the
histogram names expected by the average correction loader.

## Validation Plan

The final method should be validated on simulation, not only on real data:

- no input distortion: the recovered correction should be consistent with zero;
- controlled input distortion: the recovered correction should close on the
  known truth map;
- compare line-PoCA and helix-PoCA results in the same Job A/B framework;
- compare weighted and unweighted accumulation;
- quantify closure per component and per detector region, including TPOT
  acceptance and extrapolated regions.

## Questions For Expert Review

1. Should a same-charge control sample be kept for QA, given that the production
   estimator uses opposite-charge pairs?
2. Should the default pair weight remain curvature-proxy based, or should the
   helix version use a DCA/covariance-based weight?
3. Should high-occupancy batching merge all pair rows directly, or merge
   batch-level correction estimates with batch-level uncertainties?
4. What is the preferred robust estimator for outlier crossing points:
   DCA cut, trimmed mean, median-like estimator, or iterative rejection?
5. For the ideal-helix solver, is an analytic two-helix PoCA preferred, or is a
   numerical minimization of the two helical paths acceptable for production?
6. What simulation closure tests are necessary before comparing to real data?
