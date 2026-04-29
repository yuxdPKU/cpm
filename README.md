# cpm

Crossing Point Method realization for distortion correction in the sPHENIX experiment.

## Current Prototype

This repository currently contains the framework-independent CPM core:

- ACTS-ready CPM record/data structures under `module`;
- a voxel container keyed by `(iphi, ir, iz)`;
- a first local line-line PoCA solver for CPM v1 validation;
- a small CMake test source for the PoCA and voxel-container basics.

Directory layout:

- `module/`: CPM module source files, headers, and build helper files.
- `macro/`: Fun4All running macros.
- `jobB/`: ROOT Job B macros and the `run_cpm_b_chain.sh` driver.

The first complete Job A macro is `macro/Fun4All_CPMTrackAnalysis.C`. It follows
the current two-fit/prune tracking workflow and replaces the `PHTpcResiduals`
registration block after the second Si-TPOT fit with `PHCPMTpcCalibration`.

Job A writes a flat ROOT file with:

- `cpm_records`: one ACTS-ready TPC state/voxel record per row;
- `cpm_metadata`: grid and counter metadata for the Job A segment.

For the cluster-DST production mode, Job A stores `cluster_source` as the input
cluster DST and `track_source` as the optional CPM mini-DST containing
`SvtxSiliconMMTrackMap`. The mini-DST intentionally does not store
`TRKR_CLUSTER`; clusters are recovered from `cluster_source` if a later
rehydration pass needs them.

`jobB/CPM_B0_BuildEventIndex.C` reads one or more Job A output files and
builds:

- `cpm_event_requests`: unique events sorted by source/run/segment/event keys;
- `cpm_object_requests`: requested track/state/cluster references per event.

`jobB/CPM_B0_CheckEventIndex.C` performs a light QA pass on that index before
mini-DST rehydration is attempted.

`jobB/CPM_B1_LocalLinePoCA.C` reads one or more Job A outputs, groups records
by voxel, applies the intra-voxel offset shift, forms same-charge local
line-line PoCA pairs, and writes the pair QA tree. It stores a pair weight
`1/(pt_a*pt_b)`, which is proportional to the method weight
`(1/R_a)*(1/R_b)` for a fixed magnetic field. By default it prints one
diagnostic line per voxel with `(iphi, ir, iz)`, total `(phi, r, z)` bins,
record count, unique track count, unique track-pair count, same-charge counts,
candidate pairs, accepted pairs, and the processing status. It also writes
`cpm_b1_voxel_summary`, a persistent per-voxel QA tree. Optional pair-input
controls can require `pt >= --b1-min-pair-pt`, keep only the record closest to
the voxel center for each unique track, and then keep a deterministic hash
sample of at most `--b1-max-pair-records` records per voxel before forming
pairs.

`jobB/CPM_B2_AccumulateVoxelCorrections.C` reads one or more B1 outputs and
accumulates pair-level PoCA deltas into voxel-level correction QA rows using
that curvature-proxy weighted average.
The B1/B2 delta convention is `voxel center - crossing point`, matching the
distortion values subtracted by `TpcDistortionCorrection`.

`jobB/CPM_B3_WriteAverageCorrectionHistograms.C` converts the B2 voxel rows
into average-correction histograms named like `hIntDistortionR_{negz,posz}`,
`hIntDistortionP_{negz,posz}`, and `hIntDistortionZ_{negz,posz}`. It uses the
same `TpcSpaceChargeReconstructionHelper::split` and guard-bin handling as
`TpcSpaceChargeMatrixInversion`. For average corrections, `hIntDistortionP` is
filled with `mean_delta_phi` in radians, consistent with
`G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = true`. The phi axis follows the
existing `PHTpcResiduals` convention `[0, 2pi]`.

Example B0/B1 preflight:

```sh
root -l -b -q 'jobB/CPM_B0_BuildEventIndex.C("jobA_CPMVoxelContainer.root","CPM_B0_event_index.root")'
root -l -b -q 'jobB/CPM_B0_BuildEventIndex.C("cpm_filelist.txt","CPM_B0_event_index.root",true)'
root -l -b -q 'jobB/CPM_B0_CheckEventIndex.C("CPM_B0_event_index.root")'
root -l -b -q 'jobB/CPM_B1_LocalLinePoCA.C("jobA_CPMVoxelContainer.root","CPM_B1_local_line_poca.root")'
root -l -b -q 'jobB/CPM_B2_AccumulateVoxelCorrections.C("CPM_B1_local_line_poca.root","CPM_B2_voxel_corrections.root")'
root -l -b -q 'jobB/CPM_B3_WriteAverageCorrectionHistograms.C("CPM_B2_voxel_corrections.root","CPM_B3_average_correction_histograms.root","jobA_CPMVoxelContainer.root")'
root -l -b -q 'jobB/CPM_B3_CheckAverageCorrectionHistograms.C("CPM_B3_average_correction_histograms.root")'
```

For Condor production, run Job A once per DST/segment and write one
`*_CPMVoxelContainer.root` per job. Put those output filenames in
`cpm_filelist.txt`, one file per line, then build the B0 index from the list.

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
  --b1-min-records-per-charge 10 \
  --b1-min-pair-pt 1.0 \
  --b1-max-pair-records 100 \
  --run-b0-qa \
  --no-keep-intermediates
```

By default the driver skips B0 QA, keeps the B1/B2 intermediate ROOT files,
and also writes a combined file named `OUT_DIR/PREFIX_B.root` containing the
B1, B2, and B3 outputs. Use `--run-b0-qa` to include the B0 index/check step
and include its output in the combined file. Use `--no-combined-output` to
disable the merged file, or `--no-keep-intermediates` to remove B1/B2 after the
combined file and B3 correction map are written. B1 per-voxel diagnostic
printing is enabled by default; use `--no-b1-print-voxel-summary` for quieter
large production runs.

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
