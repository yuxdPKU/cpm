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

`macro/CPM_B0_BuildEventIndex.C` reads one or more Job A output files and
builds:

- `cpm_event_requests`: unique events sorted by source/run/segment/event keys;
- `cpm_object_requests`: requested track/state/cluster references per event.

`macro/CPM_B0_CheckEventIndex.C` performs a light QA pass on that index before
mini-DST rehydration is attempted.

`macro/CPM_B1_LocalLinePoCA.C` reads one or more Job A outputs, groups records
by voxel, and writes the first local line-line PoCA pair QA tree.

`macro/CPM_B2_AccumulateVoxelCorrections.C` reads one or more B1 outputs and
accumulates pair-level PoCA deltas into voxel-level correction QA rows.

Example B0/B1 preflight:

```sh
root -l -b -q 'macro/CPM_B0_BuildEventIndex.C("jobA_CPMVoxelContainer.root","CPM_B0_event_index.root")'
root -l -b -q 'macro/CPM_B0_BuildEventIndex.C("cpm_filelist.txt","CPM_B0_event_index.root",true)'
root -l -b -q 'macro/CPM_B0_CheckEventIndex.C("CPM_B0_event_index.root")'
root -l -b -q 'macro/CPM_B1_LocalLinePoCA.C("jobA_CPMVoxelContainer.root","CPM_B1_local_line_poca.root")'
root -l -b -q 'macro/CPM_B2_AccumulateVoxelCorrections.C("CPM_B1_local_line_poca.root","CPM_B2_voxel_corrections.root")'
```

For Condor production, run Job A once per DST/segment and write one
`*_CPMVoxelContainer.root` per job. Put those output filenames in
`cpm_filelist.txt`, one file per line, then build the B0 index from the list.

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
