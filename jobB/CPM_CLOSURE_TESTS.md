# CPM Closure Tests

This note records the no-input-distortion CPM closure checks for the simulated
ACTS-fit and GenFit samples. The goal is to understand why the final B3
`hIntDistortionR` and `hIntDistortionZ` slices can be visibly offset from zero
even when the simulation macro disables the static distortion input.

## Input Configuration

The simulation reconstruction macro disables the static distortion input:

- `G4TPC::ENABLE_STATIC_DISTORTIONS = false`
- `G4TPC::ENABLE_STATIC_CORRECTIONS = false`

The closure expectation is therefore that the recovered average CPM correction
should be close to zero, up to reconstruction bias, finite statistics, magnetic
field approximation, pair selection, and residual solver effects.

The CPM sign convention used here is:

```text
distortion = voxel center - crossing point
```

For the final B3 histograms, `hIntDistortionP` is `delta_phi` in radians, not
`r * delta_phi`.

## Fixed Slice

The closure slice diagnostic follows the same style as
`jobB/plot/draw1D_r_from3D.C`:

- fixed phi near `4.7 rad`;
- fixed positive and negative z bins near `+/-10 cm`;
- scan the R bins.

This makes the new diagnostic plots directly comparable with the existing
`resid_vsR_from3D_atZ...` plots.

## Diagnostic Quantities

The B1 pair tree stores the two closest helix points and their midpoint:

- `point_a`: PoCA point on the positive-charge track;
- `point_b`: PoCA point on the negative-charge track;
- `midpoint`: pair crossing estimate used by CPM;
- `dca`: distance between `point_a` and `point_b`.

The closure slice plots compare all three against the same voxel center:

```text
delta(point_a) = voxel center - point_a
delta(point_b) = voxel center - point_b
delta(midpoint) = voxel center - midpoint
```

If `point_a` and `point_b` have opposite biases that mostly cancel, the midpoint
can still close well. If both points share the same bias, the midpoint will keep
that bias. If the midpoint changes substantially under tighter DCA cuts, the
pair-quality selection is likely contributing.

## Running The Test

The QA chain must have been run with `--write-pair-tree`, because this test reads
`cpm_poca_pairs`.

Standalone command for GenFit:

```sh
cd /sphenix/u/xyu3/workarea/cpm/jobB

./run_cpm_closure_slice.sh \
  --qa-dir output/sim_genfit_unweighted_qa \
  --prefix sim_genfit_unweighted \
  --plot-dir output/sim_genfit_unweighted_qa/closure_slice_phi4p70_z10 \
  --phi 4.7 \
  --z 10 \
  --dca-thresholds 2.0,1.0,0.5,0.2,0.1,0.05,0.02
```

Standalone command for ACTS:

```sh
cd /sphenix/u/xyu3/workarea/cpm/jobB

./run_cpm_closure_slice.sh \
  --qa-dir output/sim_acts_unweighted_qa \
  --prefix sim_acts_unweighted \
  --plot-dir output/sim_acts_unweighted_qa/closure_slice_phi4p70_z10 \
  --phi 4.7 \
  --z 10 \
  --dca-thresholds 2.0,1.0,0.5,0.2,0.1,0.05,0.02
```

Condor command:

```sh
cd /sphenix/u/xyu3/workarea/cpm/jobB
condor_submit condor-cpm-closure.job
```

The default Condor queue is listed in `closure_jobs.txt`.

## Output Files

The test writes only to closure-specific directories, so the original QA plots
and B3 files are not overwritten.

The delta-plot y-axis ranges are fixed to match
`jobB/plot/draw1D_r_from3D.C` for simulation R-scan comparisons:

- `delta_r`: `[-0.5, 0.5] cm`;
- `delta_phi`: `[-0.01, 0.03] rad`;
- `delta_z`: `[-0.1, 0.1] cm`.

For each sample, the output directory contains:

- `PREFIX_closure_slice_profiles.root`: all profiles plus a
  `cpm_closure_slice_summary` tree;
- `PREFIX_closure_points_posz_dca2p00.pdf`: positive-z point/midpoint
  comparison using the loose `2 cm` DCA cut;
- `PREFIX_closure_points_negz_dca2p00.pdf`: negative-z point/midpoint
  comparison using the loose `2 cm` DCA cut;
- `PREFIX_closure_midpoint_delta_r_posz_dca_scan.pdf`: positive-z midpoint
  radial closure for each DCA threshold, with the original B3 slice overlaid;
- `PREFIX_closure_midpoint_delta_phi_posz_dca_scan.pdf`: positive-z midpoint
  phi closure for each DCA threshold, with the original B3 slice overlaid;
- `PREFIX_closure_midpoint_delta_z_posz_dca_scan.pdf`: positive-z midpoint
  z closure for each DCA threshold, with the original B3 slice overlaid;
- corresponding `negz` DCA-scan plots;
- `PREFIX_closure_entries_posz_dca_scan.pdf` and
  `PREFIX_closure_entries_negz_dca_scan.pdf`: accepted pair counts after each
  DCA refilter.

Each closure PDF has a companion `*_legend.pdf` file. The data PDFs do not draw
legends inside the plot area, which avoids hiding the R-scan curves when many
DCA thresholds are shown. The entries/statistics PDFs use log-y with a minimum
of `1e3`.

The latest legend-separated tight-DCA outputs are:

- `output/sim_acts_unweighted_qa/closure_slice_phi4p70_z10_tightdca_nolegend`
- `output/sim_genfit_unweighted_qa/closure_slice_phi4p70_z10_tightdca_nolegend`

Each latest directory has 10 data PDFs, 10 companion legend PDFs, and one
profile ROOT file.

## Test Log

| Date | Sample | Command | Status | Notes |
| --- | --- | --- | --- | --- |
| 2026-06-02 | sim_acts_unweighted | `condor-cpm-closure.job`, cluster `2768145.0` | done | normal termination, return value 0 |
| 2026-06-02 | sim_genfit_unweighted | `condor-cpm-closure.job`, cluster `2768145.1` | done | normal termination, return value 0 |
| 2026-06-02 | sim_acts_unweighted | `condor-cpm-closure.job`, cluster `2768147.0` | done | tight DCA scan, plot-matched y ranges |
| 2026-06-02 | sim_genfit_unweighted | `condor-cpm-closure.job`, cluster `2768147.1` | done | tight DCA scan, plot-matched y ranges |
| 2026-06-02 | sim_acts_unweighted | `condor-cpm-closure.job`, cluster `2768148.0` | done | legends moved to companion PDFs |
| 2026-06-02 | sim_genfit_unweighted | `condor-cpm-closure.job`, cluster `2768148.1` | done | legends moved to companion PDFs |

## Results From The First Closure Slice

Slice selection:

- requested `phi = 4.7 rad`, selected phi bin center `4.62512 rad`;
- requested `|z| = 10 cm`, selected z bin centers `+8.97794 cm` and
  `-8.97794 cm`.

Pair-row statistics:

| Sample | Scanned B1 pair rows | Selected pair rows before DCA scan |
| --- | ---: | ---: |
| sim_acts_unweighted | 84,854,056 | 734,861 |
| sim_genfit_unweighted | 60,904,087 | 547,463 |

Pair-entry-weighted slice means for the midpoint are:

| Sample | Side | DCA cut | Pairs | mean delta_r [cm] | mean delta_phi [rad] | mean delta_z [cm] |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| sim_acts_unweighted | +Z | `<= 2.0 cm` | 369,950 | 0.177351 | -0.000607912 | 0.0315611 |
| sim_acts_unweighted | +Z | `<= 0.2 cm` | 181,057 | 0.139942 | -0.000117550 | 0.0275570 |
| sim_acts_unweighted | -Z | `<= 2.0 cm` | 364,911 | 0.197194 | -0.000511747 | -0.0387139 |
| sim_acts_unweighted | -Z | `<= 0.2 cm` | 179,688 | 0.158782 | -0.000101995 | -0.0334790 |
| sim_genfit_unweighted | +Z | `<= 2.0 cm` | 278,340 | 0.0871131 | -0.000910439 | 0.0160473 |
| sim_genfit_unweighted | +Z | `<= 0.2 cm` | 237,296 | 0.0295937 | -0.000124754 | 0.00671098 |
| sim_genfit_unweighted | -Z | `<= 2.0 cm` | 269,123 | 0.0975554 | -0.000517532 | -0.0190846 |
| sim_genfit_unweighted | -Z | `<= 0.2 cm` | 229,853 | 0.0313060 | -0.0000998302 | -0.0077664 |

For the loose `DCA <= 2.0 cm` selection, the B3 overlay extracted from
`PREFIX_B3_average_correction_histograms.root` matches the B1 midpoint profile
means exactly within print precision. This confirms that the observed offset is
already present in the B1 crossing-point estimate and is not introduced by B2
accumulation or B3 histogram splitting/guard-bin handling.

The `point_a` and `point_b` means are very close to each other. For example, in
the ACTS +Z slice with `DCA <= 2.0 cm`, `point_a` and `point_b` have mean
`delta_r = 0.177287 cm` and `0.177255 cm`, while the midpoint has
`delta_r = 0.177351 cm`. This indicates that the midpoint is mostly inheriting a
common same-direction bias of the two helix PoCA points rather than averaging
away opposite point biases.

Tightening the pair DCA cut improves the closure but does not remove the bias,
especially for ACTS. GenFit improves more strongly under the same `0.2 cm` DCA
refilter. The current evidence therefore points to B1-level pair quality and/or
track/helix-extrapolation bias as the next place to investigate, rather than
B2/B3 averaging.

## Tight DCA Scan With Plot-Matched Axes

After matching the delta-plot y-axis ranges to `draw1D_r_from3D.C`, the tight
DCA scan was rerun into new output directories:

- `output/sim_acts_unweighted_qa/closure_slice_phi4p70_z10_tightdca_plotrange`
- `output/sim_genfit_unweighted_qa/closure_slice_phi4p70_z10_tightdca_plotrange`

Each directory contains 10 PDF plots plus one
`PREFIX_closure_slice_profiles.root` file. The earlier
`closure_slice_phi4p70_z10` plots are unchanged.

Pair-entry-weighted midpoint means for the tight DCA thresholds:

| Sample | Side | DCA cut | Pairs | mean delta_r [cm] | mean delta_phi [rad] | mean delta_z [cm] |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| sim_acts_unweighted | +Z | `<= 0.20 cm` | 181,057 | 0.139942 | -0.000117550 | 0.0275570 |
| sim_acts_unweighted | +Z | `<= 0.10 cm` | 103,703 | 0.144834 | -0.000108427 | 0.0291459 |
| sim_acts_unweighted | +Z | `<= 0.05 cm` | 54,249 | 0.144757 | -0.000131250 | 0.0312075 |
| sim_acts_unweighted | +Z | `<= 0.02 cm` | 22,340 | 0.144987 | -0.000107280 | 0.0317233 |
| sim_acts_unweighted | -Z | `<= 0.20 cm` | 179,688 | 0.158782 | -0.000101995 | -0.0334790 |
| sim_acts_unweighted | -Z | `<= 0.10 cm` | 102,374 | 0.153138 | -0.0000931793 | -0.0328488 |
| sim_acts_unweighted | -Z | `<= 0.05 cm` | 53,703 | 0.155859 | -0.000112890 | -0.0334964 |
| sim_acts_unweighted | -Z | `<= 0.02 cm` | 21,876 | 0.157899 | 0.000112974 | -0.0354721 |
| sim_genfit_unweighted | +Z | `<= 0.20 cm` | 237,296 | 0.0295937 | -0.000124754 | 0.00671098 |
| sim_genfit_unweighted | +Z | `<= 0.10 cm` | 192,007 | 0.0288544 | -0.000112073 | 0.00669465 |
| sim_genfit_unweighted | +Z | `<= 0.05 cm` | 121,219 | 0.0290694 | -0.0000831481 | 0.00658041 |
| sim_genfit_unweighted | +Z | `<= 0.02 cm` | 52,672 | 0.0249674 | -0.0000855785 | 0.00575622 |
| sim_genfit_unweighted | -Z | `<= 0.20 cm` | 229,853 | 0.0313060 | -0.0000998302 | -0.0077664 |
| sim_genfit_unweighted | -Z | `<= 0.10 cm` | 186,298 | 0.0329873 | -0.0000355916 | -0.00774196 |
| sim_genfit_unweighted | -Z | `<= 0.05 cm` | 117,897 | 0.0294781 | -0.0000618612 | -0.00680631 |
| sim_genfit_unweighted | -Z | `<= 0.02 cm` | 51,643 | 0.0273260 | -0.0000493642 | -0.00632211 |

The tight DCA cuts preserve enough statistics in this slice. At
`DCA <= 0.02 cm`, ACTS still has about `2.2e4` accepted pairs per z side, while
GenFit has about `5.2e4`. GenFit continues to move modestly toward zero under
the tightest DCA cut. ACTS improves strongly when moving from `2.0 cm` to
`0.2 cm`, but the `0.2 -> 0.02 cm` scan mostly plateaus and does not eliminate
the remaining radial/z offset. This suggests that the strict pair-crossing
requirement helps, but the ACTS residual closure issue likely also contains a
track-extrapolation or solver-model component beyond large-DCA pairs.

## Interpretation Checklist

Use the following checks when reading the plots:

1. Compare `point_a`, `point_b`, and `midpoint` in the `closure_points_*`
   plots. A midpoint-only bias points to the crossing/midpoint construction; a
   shared point bias points upstream to track extrapolation or voxel association.
2. Compare `DCA <= 2.0, 1.0, 0.5, 0.2 cm` in the midpoint DCA-scan plots. If
   tightening the DCA cut moves the curve toward zero, the pair-quality cut is a
   major part of the closure issue.
3. Compare the B1 midpoint profile with the overlaid B3 slice. Agreement means
   B2/B3 averaging and histogram writing are not introducing the offset. A
   mismatch would point to accumulation or split/guard-bin handling.
4. Watch the entries plots while tightening DCA. A cleaner but statistically
   thin curve should not be overinterpreted.
