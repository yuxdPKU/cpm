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
  length, covariance, and optional ACTS surface diagnostics.

The Job A output is segment-friendly: each Condor job writes one compact CPM
ROOT file, and Job B can read either a single file or a file list.

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
- `--b1-max-pair-records N` means at most `N` selected positive records and
  `N` selected negative records per batch.
- `N = 0` means one unlimited full-voxel batch.
- B1 writes per-voxel QA including total selected records, selected positive and
  negative records, batch count, batched positive/negative records, theoretical
  opposite-charge pairs, batched opposite-charge pairs, candidate pairs, and
  accepted pairs.

Current crossing solver:

- The current framework-building solver is local line-line PoCA using the two
  reference-state positions and momentum directions.
- The pair crossing estimate is the midpoint of the two closest points.
- The pair DCA is stored as QA and can be cut in B1/B2.

Target crossing solver:

- The target solver is an ideal two-helix PoCA/crossing method.
- A local coresoftware survey found no ready-made two-helix PoCA utility.
  `TrackFitUtils` provides helix-to-point/surface helpers, while
  `KshortReconstruction` and `PHSimpleVertexFinder` provide line-line two-track
  PCA/DCA.
- CPM should therefore add a CPM-local helix PoCA module first, with a narrow
  interface designed for later upstreaming.

## Job B2: Voxel Accumulation

B2 reads accepted B1 pair rows and accumulates one correction row per voxel.

Implemented averaging modes:

- weighted mode: use the B1 curvature-proxy weight `1/(pt_a * pt_b)`;
- unweighted mode: use unit weights for a simple arithmetic mean.

B2 stores per-voxel QA including entries, sum of weights, effective weighted
entries, means, and RMS values for `delta_r`, `delta_phi`, `delta_rphi`,
`delta_z`, and DCA.

Open accumulation question:

- The current implementation accumulates all accepted pair rows directly.
- The expert-suggested running-average version may instead compute one
  correction and uncertainty per batch, then merge batch-level results into the
  voxel accumulator. This needs a physics/statistics decision because averaging
  all pairs and averaging batch means are not identical when batch sizes or
  weights differ.

## Job B3: Average-Correction Histograms

B3 converts B2 voxel rows into average-correction histograms. It uses the Job A
metadata for grid dimensions, fills reconstructed QA histograms, splits negative
and positive z, applies the existing guard-bin style, and writes the histogram
names expected by the average correction loader.

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

1. Is opposite-charge pairing always preferred, or should there be a diagnostic
   same-charge control sample?
2. Should the default pair weight remain curvature-proxy based, or should the
   helix version use a DCA/covariance-based weight?
3. Should high-occupancy batching merge all pair rows directly, or merge
   batch-level correction estimates with batch-level uncertainties?
4. What is the preferred robust estimator for outlier crossing points:
   DCA cut, trimmed mean, median-like estimator, or iterative rejection?
5. For the ideal-helix solver, is an analytic two-helix PoCA preferred, or is a
   numerical minimization of the two helical paths acceptable for production?
6. What simulation closure tests are necessary before comparing to real data?
