# cpm

Crossing Point Method realization for distortion correction in the sPHENIX experiment.

## Current Prototype

This repository currently contains the framework-independent CPM core:

- ACTS-ready CPM record/data structures under `module`;
- a voxel container keyed by `(iphi, ir, iz)`;
- a mergeable CPM correction container aligned with the matrix-inversion
  container workflow;
- a module-backed average-correction reconstruction step that merges Job A
  containers and writes the final histogram file;
- a solver-agnostic B1 PoCA driver with selectable `helix` and `line` modes;
- an initial CPM-local ideal-helix PoCA helper for review and validation;
- a small CMake test source for the PoCA and voxel-container basics.

Directory layout:

- `module/`: CPM module source files, headers, and build helper files.
- `macro/`: Fun4All running macros.
- `jobB/`: ROOT Job B macros and the `run_cpm_b_chain.sh` driver.
- `CPM_PROCEDURE_FOR_REVIEW.md`: living method note for external review.

The first complete Job A macro is `macro/Fun4All_CPMTrackAnalysis.C`. It follows
the current two-fit/prune tracking workflow and replaces the `PHTpcResiduals`
registration block after the second Si-TPOT fit with `PHCPMTpcCalibration`.

Job A writes a ROOT file with:

- `CPMCorrectionContainer`: additive per-voxel entries, pair-weight sums,
  unweighted delta sums, weighted delta sums, and QA sums filled during
  reconstruction in running opposite-charge batches;
- `cpm_metadata`: grid and counter metadata for the Job A segment.
- optional `cpm_records`: one ACTS-ready TPC state/voxel record per row, written
  only when `writeCpmRecords=true` is passed to the Job A macro.

The running container is the preferred production path. It follows the
`TpcSpaceChargeMatrixContainer` style: each segment output is independently
mergeable, and the offline step can choose either plain arithmetic means or
pair-weighted means from the same persisted sums. A voxel is flushed when both
charge queues reach the configured running threshold; the flush pairs all
currently queued positive records with all currently queued negative records,
then clears that voxel queue. The optional `cpm_records` tree is treated as a
QA/diagnostic product for B0/B1 studies and is disabled by default in the Job A
macro.

For the cluster-DST production mode, Job A stores `cluster_source` as the input
cluster DST and `track_source` as the optional CPM mini-DST containing
`SvtxSiliconMMTrackMap`. The mini-DST intentionally does not store
`TRKR_CLUSTER`; clusters are recovered from `cluster_source` if a later
rehydration pass needs them.

## Job B Macro Map

The `jobB/CPM_B*`, `jobB/CPM_QA_*`, and reconstruction macros are the offline CPM
post-processing tools. There are two supported routes:

```text
Recommended production route:
Job A output with CPMCorrectionContainer
  -> CPM_ReconstructAverageCorrection.C
  -> CPM_B3_CheckAverageCorrectionHistograms.C

Diagnostic/offline-PoCA route:
Job A cpm_records, with writeCpmRecords=true
  -> CPM_QA_RunOfflineDiagnostics.C
```

`CPM_QA_RunOfflineDiagnostics.C` is the single user-facing macro for the
record-based offline diagnostic path. It can run optional B0 event-index QA,
then writes separate stage files for B1 PoCA, B2 voxel accumulation, and B3
histogram output so each stage can be inspected offline.

`CPM_QA_B0_BuildEventIndex.C` and `CPM_QA_B0_CheckEventIndex.C` are optional
rehydration QA tools. They scan `cpm_records`, build an event-ordered list of
requested tracks/states/clusters, and check that the resulting event/object
index is self-consistent before any later mini-DST readback study. These macros
require Job A to have been run with `writeCpmRecords=true`.

`CPM_QA_B1_ComputePoCA.C` is the offline PoCA diagnostic step. It reads
`cpm_records`, groups records by voxel, applies the intra-voxel offset shift,
forms opposite-charge pairs, and computes crossing points with either the
default helix solver or the local line solver. It can write detailed pair rows
and batch-level correction sums for solver comparisons and QA.
It also requires `writeCpmRecords=true` in Job A.

`CPM_QA_B2_AccumulateVoxelCorrections.C` is the diagnostic B2 paired with B1. It
reads B1 pair rows or B1 batch sums and produces one `cpm_voxel_corrections`
row per voxel. It is useful for checking the offline PoCA path, but it is not
the preferred production merger now that Job A writes `CPMCorrectionContainer`.

`CPM_ReconstructAverageCorrection.C` is the recommended production entry point.
It is a thin macro around `CPMAverageCorrectionReconstruction`, which is the CPM
counterpart of the matrix-inversion offline reconstruction class. It reads the
`CPMCorrectionContainer` objects written by Job A, verifies compatible grids
through the container merge guard, adds the stored per-voxel sums, chooses either
pair-weighted or plain arithmetic voxel means, and writes one ROOT file with the
merged container, `cpm_voxel_corrections`, `cpm_metadata`, reconstructed QA
histograms, final guarded average-correction histograms, and a summary tree.

`CPM_B2_MergeCorrectionContainers.C` is now a lower-level compatibility/debug
macro for the split-stage workflow. Its most important output is the
`cpm_voxel_corrections` TTree consumed by B3; the merged container is retained
as an auxiliary object for inspection.

`CPM_B3_WriteAverageCorrectionHistograms.C` converts `cpm_voxel_corrections`
into the average-correction histograms consumed by the distortion correction
loader. It writes the guarded half-TPC histograms
`hIntDistortion{R,P,Z}_{negz,posz}` and `hentries_{negz,posz}`.

`CPM_B3_CheckAverageCorrectionHistograms.C` is the final lightweight QA check.
It verifies that the required B3 histograms exist and have valid dimensions, so
batch output failures can be caught before the correction map is used downstream.

`jobB/CPM_QA_B0_BuildEventIndex.C` reads one or more Job A output files and
builds:

- `cpm_event_requests`: unique events sorted by source/run/segment/event keys;
- `cpm_object_requests`: requested track/state/cluster references per event.

`jobB/CPM_QA_B0_CheckEventIndex.C` performs a light QA pass on that index before
mini-DST rehydration is attempted.

`jobB/CPM_QA_B1_ComputePoCA.C` reads one or more Job A outputs, groups records
by voxel, applies the intra-voxel offset shift, forms opposite-charge crossing
pairs, and writes batch-level correction sums plus optional pair QA rows. It
stores a pair weight `1/(pt_a*pt_b)`, which is proportional to the method weight
`(1/R_a)*(1/R_b)` for a fixed magnetic field. By default it prints one
diagnostic line per voxel with `(iphi, ir, iz)`, total `(phi, r, z)` bins,
record count, unique track count, unique track-pair count, charge-pair counts,
batch counts, batched pair counts, candidate pairs, accepted pairs, and the
processing status. It also writes `cpm_b1_voxel_summary`, a persistent
per-voxel QA tree, and `cpm_b1_batch_corrections`, a persistent per-batch tree
with unweighted and weighted sums, sum-of-squares, accepted pair counts, total
weights, and DCA QA. Optional pair-input controls can require
`pt >= --b1-min-pair-pt`, keep only the record closest to the voxel center for
each unique track, and then process the selected records in deterministic
hash-ordered opposite-charge batches. The raw per-voxel safety skip is disabled
by default with `--b1-max-records 0`; high occupancy is controlled by batch
size instead of skipping the voxel. `--b1-max-pair-records` sets the maximum
selected records per charge sign in each batch; the default is `10`, while `0`
keeps one unlimited full-voxel batch. The default crossing solver is the CPM
ideal-helix PoCA helper using a configurable first-pass uniform `Bz` estimate
that defaults to `1.4 T`. The previous v1 local line-line solver can still be
selected with `--b1-crossing-solver line` for comparison studies. Use
`--no-b1-write-pairs` to skip the large pair-level QA tree when only batch-level
accumulation is needed.

`jobB/CPM_QA_B2_AccumulateVoxelCorrections.C` reads one or more B1 outputs and
accumulates crossing deltas into voxel-level correction QA rows. In `auto` mode
it prefers `cpm_b1_batch_corrections` and falls back to pair rows for older B1
outputs. The batch mode preserves pair-equivalent weighted and unweighted sums
without requiring every pair row to be persisted. B2 can run either the default
curvature-proxy weighted average or a simple unweighted average through the Job
B driver.
The B1/B2 delta convention is `voxel center - crossing point`, matching the
distortion values subtracted by `TpcDistortionCorrection`.

`jobB/CPM_ReconstructAverageCorrection.C` is the production macro for the
container path. It calls `CPMAverageCorrectionReconstruction`, which owns the
formerly split B2/B3 production work: merging Job A containers, calculating
weighted or plain average voxel corrections, writing the voxel QA TTree, and
writing the guarded final average-correction histograms in the same output file.

`jobB/CPM_B2_MergeCorrectionContainers.C` is the lower-level split-stage
container debug macro. It reads one or more Job A files, merges their
`CPMCorrectionContainer` objects with grid/range guards, and writes the
`cpm_voxel_corrections` tree consumed by B3. The `use_pair_weights` argument
selects weighted or plain averaging at merge time without rerunning Job A.

`jobB/CPM_B3_WriteAverageCorrectionHistograms.C` converts the B2 voxel rows
into average-correction histograms named like `hIntDistortionR_{negz,posz}`,
`hIntDistortionP_{negz,posz}`, and `hIntDistortionZ_{negz,posz}`. It uses the
same `TpcSpaceChargeReconstructionHelper::split` and guard-bin handling as
`TpcSpaceChargeMatrixInversion`. For average corrections, `hIntDistortionP` is
filled with `mean_delta_phi` in radians, consistent with
`G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = true`. The phi axis follows the
existing `PHTpcResiduals` convention `[0, 2pi]`.

Example QA/offline-PoCA path:

These commands require Job A output produced with `writeCpmRecords=true`.

```sh
root -l -b -q 'jobB/CPM_QA_RunOfflineDiagnostics.C("jobA_CPMVoxelContainer.root","qa_output","seg0",false,true)'
root -l -b -q 'jobB/CPM_QA_RunOfflineDiagnostics.C("cpm_filelist.txt","qa_output","merged",true,true)'

# Lower-level stage macros are still available when a single stage needs to be rerun.
root -l -b -q 'jobB/CPM_QA_B0_BuildEventIndex.C("jobA_CPMVoxelContainer.root","CPM_QA_B0_event_index.root")'
root -l -b -q 'jobB/CPM_QA_B0_BuildEventIndex.C("cpm_filelist.txt","CPM_QA_B0_event_index.root",true)'
root -l -b -q 'jobB/CPM_QA_B0_CheckEventIndex.C("CPM_QA_B0_event_index.root")'
root -l -b -q 'jobB/CPM_QA_B1_ComputePoCA.C("jobA_CPMVoxelContainer.root","CPM_QA_B1_poca.root")'
root -l -b -q 'jobB/CPM_QA_B2_AccumulateVoxelCorrections.C("CPM_QA_B1_poca.root","CPM_QA_B2_voxel_corrections.root")'
root -l -b -q 'jobB/CPM_B3_WriteAverageCorrectionHistograms.C("CPM_QA_B2_voxel_corrections.root","CPM_B3_average_correction_histograms.root","jobA_CPMVoxelContainer.root")'
root -l -b -q 'jobB/CPM_B3_CheckAverageCorrectionHistograms.C("CPM_B3_average_correction_histograms.root")'
```

Example container-aligned production path:

```sh
root -l -b -q 'jobB/CPM_ReconstructAverageCorrection.C("jobA_CPMVoxelContainer.root","CPM_average_correction_weighted.root","CPMCorrectionContainer",true)'
root -l -b -q 'jobB/CPM_ReconstructAverageCorrectionFromList("cpm_filelist.txt","CPM_average_correction_plain.root","CPMCorrectionContainer",false)'
root -l -b -q 'jobB/CPM_B3_CheckAverageCorrectionHistograms.C("CPM_average_correction_weighted.root")'
```

The split-stage container debug path is still available if the intermediate
voxel TTree needs to be inspected separately:

```sh
root -l -b -q 'jobB/CPM_B2_MergeCorrectionContainers.C("jobA_CPMVoxelContainer.root","CPM_B2_voxel_corrections.root","CPMCorrectionContainer",true)'
root -l -b -q 'jobB/CPM_B3_WriteAverageCorrectionHistograms.C("CPM_B2_voxel_corrections.root","CPM_B3_average_correction_histograms.root","CPM_B2_voxel_corrections.root")'
```

For Condor production, run Job A once per DST/segment and write one
`*_CPMVoxelContainer.root` per job. Put those output filenames in
`cpm_filelist.txt`, one file per line, then run the container production macro
or the `--b2-containers` driver path on that list.
For ordinary container production, keep `writeCpmRecords=false`; pass `true` as
the seventh argument to `macro/run_clusterdst.sh` only when B0/B1 diagnostic
records are needed.

The full Job B chain can also be run with:

```sh
cd ~/workarea/cpm

jobB/run_cpm_b_chain.sh \
  --input macro/root/Reconstructed/79516/clusters_seeds_79516-0.root_CPMVoxelContainer.root \
  --out-dir output/jobB/run79516 \
  --prefix seg0 \
  --no-keep-intermediates

jobB/run_cpm_b_chain.sh \
  --input cpm_filelist.txt \
  --input-is-list \
  --out-dir output/jobB/run79516 \
  --prefix merged \
  --b2-containers \
  --b2-weighted \
  --run-b0-qa \
  --no-keep-intermediates

jobB/run_cpm_b_chain.sh \
  --input cpm_filelist.txt \
  --input-is-list \
  --out-dir output/jobB/run79516 \
  --prefix merged_plain \
  --b2-containers \
  --b2-unweighted \
  --no-keep-intermediates
```

By default the driver skips B0 QA, keeps the B1/B2 intermediate ROOT files,
and also writes a combined file named `OUT_DIR/PREFIX_B.root` containing the
B1, B2, and B3 outputs. With `--b2-containers`, the driver skips B1 and the
diagnostic B2 macro; the module-backed reconstruction writes the final B3-style
ROOT file directly, and the combined file contains that output plus optional B0
QA. Use
`--run-b0-qa` to include the B0 index/check step and include its output in the
combined file. Use `--no-combined-output` to disable the merged file, or
`--no-keep-intermediates` to remove B1/B2 after the combined file and B3
correction map are written. B1 per-voxel diagnostic printing is enabled by
default; use `--no-b1-print-voxel-summary` for quieter diagnostic B1 runs.

The recommended convention is to launch the script from the repository root
instead of from `macro/`, and to write ROOT outputs into a dedicated output
directory such as `output/jobB/run79516/`. This keeps `macro/` reserved for
source macros and avoids mixing `.C` files with generated `.root` files.

## Build and Test

Compilation should be checked inside an sPHENIX software environment. The module includes an autotools skeleton matching the usual coresoftware package style:

```sh
cd module
./autogen.sh --prefix="$MYINSTALL"
make
make install
```

There is also a lightweight CMake entry point for local development:

```sh
cmake -S module -B build
cmake --build build
ctest --test-dir build --output-on-failure
```
