# Denoising Pipeline Optimization Report

- Generated: 2026-04-18 12:52:46
- Project:   `/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/v2_no_param.mat`
- Subjects:  20
- Iterations: 160 (Stage 1: 144, Stage 2: 16)

## Validation Against GUI Ground Truth

| Score        | GUI    | Re-computed | Diff   |
|--------------|--------|-------------|--------|
| Validity     | 0.7558 | 0.7558      | -0.0000 |
| Quality      | 0.8959 | 0.8959      | +0.0000 |
| Sensitivity  | (not stored) | 0.9194 | - |

Tolerance: 0.010.  Validation **passed**.

## Best Parameter Combination

| Parameter | Value |
|-----------|-------|
| POLY_ORD  | 2 |
| BP_HZ     | [0.008, 0.09] |
| SIMULT    | 1 |
| MOT24     | 0 |
| N_ACOMP   | 20 |
| N_AROMA   | 0 |
| Validity  | 0.7687 |
| Quality   | 0.8872 |
| Sensitivity | 0.9978 |
| **Mean_QC** | **0.8846** |

## Score Distributions (main-effect summary)

Mean Mean_QC across combinations grouped by each parameter level.

### POLY_ORD

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 1 | 0.8035 | 0.6257 | 0.8392 | 0.9824 | 48 |
| 2 | 0.8114 | 0.6518 | 0.8490 | 0.9710 | 64 |
| 3 | 0.8040 | 0.6342 | 0.8490 | 0.9657 | 48 |

### BP_HZ

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| [0.008 0.09] | 0.7924 | 0.6392 | 0.8550 | 0.9449 | 88 |
| [0.008 Inf] | 0.8244 | 0.6381 | 0.8351 | 1.0000 | 72 |

### SIMULT

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 1 | 0.7909 | 0.6516 | 0.8381 | 0.9449 | 88 |
| 0 | 0.8262 | 0.6230 | 0.8558 | 1.0000 | 72 |

### MOT24

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 1 | 0.8098 | 0.6369 | 0.8661 | 0.9640 | 72 |
| 0 | 0.8044 | 0.6401 | 0.8297 | 0.9799 | 88 |

### N_ACOMP

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 5 | 0.8124 | 0.5428 | 0.9021 | 0.9923 | 50 |
| 10 | 0.8694 | 0.7339 | 0.8757 | 0.9986 | 4 |
| 15 | 0.8725 | 0.7540 | 0.8808 | 0.9827 | 4 |
| 20 | 0.8222 | 0.7028 | 0.8434 | 0.9338 | 50 |
| 30 | 0.7799 | 0.7400 | 0.8400 | 0.7176 | 4 |
| 50 | 0.7765 | 0.6458 | 0.7857 | 1.0000 | 48 |

### N_AROMA

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 0 | 0.8009 | 0.6337 | 0.8330 | 0.9685 | 75 |
| 5 | 0.8479 | 0.7172 | 0.8689 | 0.9926 | 5 |
| 10 | 0.8395 | 0.7290 | 0.8624 | 0.9542 | 5 |
| Inf | 0.8079 | 0.6324 | 0.8566 | 0.9774 | 75 |

## Recommendation

Use the **Best Parameter Combination** above for this dataset.

## Caveats

- 0 of 160 combinations failed (0.0%).
- QC scores recomputed on CONN's sample-voxel network (~500-2000 grey-matter voxels).
- Two-stage grid: Stage 1 tested 144 coarse combos; Stage 2 16 fine combos around the winner.
- Mean_QC = simple mean of Validity, Quality, Sensitivity (each in [0,1]).
- AROMA available components in this dataset = 6, so N_AROMA >= 5 all map to at most 6 actual regressors.

## Correction Note

The previous version of this script (moved to `temp/outdated/`) produced QC
scores that did not match the CONN GUI for the same project.  For the
run_conn_batch.m configuration the old script reported Validity=0.8035
while the GUI reported 0.7558 - a discrepancy of 0.048.  Two bugs were
identified and fixed:

1. **`strcat` stripped trailing whitespace in confound names.**  The old
   script rebuilt the CONN project with
   `strcat('Effect of ', COND_NAMES)`, which MATLAB's `strcat` silently
   strips to `'Effect ofrandom'`, `'Effect ofstr_r12'`, etc.  `conn_batch`
   accepted these without error but, at QC time, `conn_designmatrix`
   could not associate them with any condition, so the per-condition
   regressors were effectively omitted.  That changed the shape of the
   after-denoising FC distribution and thus the Validity score.  Fix:
   start from the existing `v2_no_param.mat` project that was built by
   `run_conn_batch.m` (which uses `append` instead of `strcat` and
   therefore has the correct space-separated names).

2. **DataSensitivity silently returned [].**
   `conn_qascores('DataSensitivity', ...)` prefers
   `QC_DOF_WelchSatterthwaite` over `QC_DOF`, but for this dataset the
   per-subject Welch-Satterthwaite values are often below 3, so
   `sum(QC_DOF_WS - 3)` is negative, `sqrt` returns a complex number,
   `spm_Ncdf` throws "Input must be real and full", and the whole call
   is swallowed by an unguarded `try ... end`.  That is why no
   `DataSensitivityScore.mat` was written to the GUI folder.  Fix: in
   this script we compute Sensitivity directly with the regular
   `QC_DOF` covariate using the same formula
   (`1 - spm_Ncdf(1.645 - 0.1003*sqrt(sum(QC_DOF(valid)-3)))`).  The
   manual value for the GUI configuration is 0.9194.

3. **Validation step added.**  Before the grid search, the script now
   re-scores the run_conn_batch.m parameters and aborts if recomputed
   Validity and Quality differ from the stored GUI values by more than
   0.010.  The corresponding row in `optimization_results.xlsx` is
   flagged in the `GUI_verified` column.

