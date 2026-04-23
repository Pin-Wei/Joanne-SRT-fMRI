# Denoising Pipeline Optimization Report

- Generated: 2026-04-19 14:41:01
- Project:   `/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/v3_no_param.mat`
- Subjects:  20
- Iterations: 160 (Stage 1: 72, Stage 2: 13)

## Validation Against GUI Ground Truth

| Score        | GUI    | Re-computed | Diff    |
|--------------|--------|-------------|---------|
| Validity     | 0.7687 | 0.7687      | +0.0000 |
| Quality      | 0.8872 | 0.8872      | +0.0000 |
| Sensitivity  | 0.5421 | 0.5421      | +0.0000 |

Tolerance: 0.010.  Validation **passed**.

## Best Parameter Combination

| Parameter | Value |
|-----------|-------|
| POLY_ORD  | 3 |
| BP_HZ     | [0.008, 0.09] |
| SIMULT    | 0 |
| MOT24     | 1 |
| N_ACOMP   | 20 |
| N_AROMA   | Inf |
| Validity  | 0.7133 |
| Quality   | 0.9120 |
| Sensitivity | 1.0000 |
| **Mean_QC** | **0.8751** |

## Score Distributions (main-effect summary)

Mean Mean_QC across combinations grouped by each parameter level.

### POLY_ORD

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 1 | 0.8017 | 0.6262 | 0.8383 | 0.9889 | 48 |
| 2 | 0.8016 | 0.6273 | 0.8444 | 0.9815 | 48 |
| 3 | 0.8149 | 0.6426 | 0.8594 | 0.9790 | 64 |

### BP_HZ

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| [0.008 0.09] | 0.7926 | 0.6290 | 0.8595 | 0.9631 | 88 |
| [0.008 Inf] | 0.8244 | 0.6381 | 0.8351 | 1.0000 | 72 |

### SIMULT

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 1 | 0.7773 | 0.6339 | 0.8332 | 0.9508 | 72 |
| 0 | 0.8312 | 0.6324 | 0.8611 | 1.0000 | 88 |

### MOT24

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 1 | 0.8200 | 0.6426 | 0.8709 | 0.9929 | 88 |
| 0 | 0.7909 | 0.6215 | 0.8212 | 0.9707 | 72 |

### N_ACOMP

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 5 | 0.8040 | 0.5404 | 0.9036 | 0.9801 | 50 |
| 10 | 0.8501 | 0.6783 | 0.8719 | 1.0000 | 4 |
| 15 | 0.8640 | 0.7013 | 0.8906 | 1.0000 | 4 |
| 20 | 0.8261 | 0.7000 | 0.8438 | 0.9652 | 50 |
| 30 | 0.8620 | 0.7097 | 0.8763 | 1.0000 | 4 |
| 50 | 0.7771 | 0.6441 | 0.7885 | 1.0000 | 48 |

### N_AROMA

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 0 | 0.7988 | 0.6314 | 0.8326 | 0.9706 | 75 |
| 5 | 0.8474 | 0.6584 | 0.8839 | 1.0000 | 5 |
| 10 | 0.8489 | 0.6636 | 0.8830 | 1.0000 | 5 |
| Inf | 0.8095 | 0.6311 | 0.8599 | 0.9930 | 75 |

## Recommendation

Use the **Best Parameter Combination** above for this dataset.

## Caveats

- 0 of 85 new combinations failed (0.0%).
- QC scores recomputed on CONN's sample-voxel network (~500-2000 grey-matter voxels).
- Two-stage grid: Stage 1 tested 72 coarse combos; Stage 2 13 fine combos around the winner.
- Mean_QC = simple mean of Validity, Quality, Sensitivity (each in [0,1]).
- AROMA available components in this dataset = 6, so N_AROMA >= 5 all map to at most 6 actual regressors.

## Correction Note

This script has been corrected twice against CONN GUI ground truth.

### Fix 1 (initial) - confound-name whitespace

The original `temp/outdated/` version rebuilt the CONN project via
`strcat('Effect of ', COND_NAMES)`.  MATLAB's `strcat` silently strips
trailing whitespace, producing `'Effect ofrandom'` etc. that
`conn_designmatrix` could not match to any condition, so the
per-condition regressors were never applied during QC denoising.
Fix: load the existing project built by `run_conn_batch.m`
(which uses `append`, preserving the space) as the template.

### Fix 2 (this pass) - DataSensitivity used the wrong DOF covariate

The `temp/bad_2/` version computed Sensitivity from `QC_DOF` (the
integer pre-denoising DOF, ~108 per subject).  That made
`sum(QC_DOF(valid)-3)` ~2000 for every combination, driving the
score to ~0.998 regardless of the denoising parameters, and
disagreeing with the GUI value of **0.5421** stored in
`v3_no_param/results/qa/.../DataSensitivityScore.mat`.

The CONN GUI source (`conn_qascores.m` line 140-143) prefers
`QC_DOF_WelchSatterthwaite` - the fractional residual DOF populated
by the `QA_COV` step (typically 10-30 per subject).  Summing
`(QC_DOF_WS(valid) - 3)` for this dataset gives ~304.7, and
`1 - spm_Ncdf(1.645 - 0.1003*sqrt(304.7)) = 0.5421` - matches GUI.

`compute_sensitivity_manual` now pulls `QC_DOF_WelchSatterthwaite`
first and only falls back to `QC_DOF` if that covariate is missing.

### Validation step

Before the grid search the script re-scores the run_conn_batch.m
configuration and aborts if recomputed Validity/Quality/Sensitivity
differ from the stored GUI values by more than 0.010.  The matching
row in `compare_QC_scores_optimized.xlsx` is flagged `GUI_verified`.

## Second Optimization Run - Gap Fill

### Audit Results

- Prior xlsx files merged: 1 (rows: 75).
- Stage 1 intended (full coarse grid): 144 combinations.
- Stage 1 newly run this pass: 72 combinations.
- Stage 2 newly run this pass: 13 combinations.
- Total new combinations: 85.

### Root Cause of Missing N_AROMA > 0 Runs

The previous run collapsed AROMA -> {0} because
`Preproc.confounds.names` lacked an `aroma` entry (the v3 project was
built by `run_conn_batch.m` with `N_AROMA = 0`, so line 47 of that
script never added `aroma` to `CONFOUND_NAMES`).  The auto-collapse
guard at line 100-106 of the prior `run_conn_optimize.m` saw the
absence and disabled all AROMA>0 testing.

However `aroma` IS registered in `CONN_x.Setup.l1covariates.names`
(index 5), with all 20 subjects x 8 sessions of TSV timeseries
already loaded by the run_conn_batch.m Setup step.  The fix:
`update_confounds()` now ADDS an `aroma` entry to
`Preproc.confounds` at runtime when `n_aroma>0` and aroma is missing
from the names list.  The aroma availability check now consults
`Setup.l1covariates.names` instead of `Preproc.confounds.names`.

### Skip Logic

A new `EXISTING_RESULT_FILES` cell array at the top of the script
lists prior xlsx files.  At startup their parameter columns are
hashed into a set; combinations already in the set are logged and
skipped.  This run loaded 75 prior rows from 1 file(s) and skipped
the corresponding combinations during the coarse grid pass.

### New Best Parameter Combination (merged set)

| Parameter | Value |
|-----------|-------|
| Origin    | new |
| POLY_ORD  | 3 |
| BP_HZ     | [0.008, 0.09] |
| SIMULT    | 0 |
| MOT24     | 1 |
| N_ACOMP   | 20 |
| N_AROMA   | Inf |
| Validity  | 0.7133 |
| Quality   | 0.9120 |
| Sensitivity | 1.0000 |
| **Mean_QC** | **0.8751** |

### Comparison: Top 10 from Prior vs. New Top 10

#### Prior Top 10

| POLY | BP | SIMULT | MOT24 | N_ACOMP | N_AROMA | Validity | Quality | Sens | Mean_QC |
|------|----|--------|-------|---------|---------|----------|---------|------|---------|
| 3 | [0.008, 0.09] | 0 | 1 | 20 | 0 | 0.7083 | 0.8857 | 1.0000 | 0.8647 |
| 3 | [0.008, 0.09] | 0 | 1 | 15 | 0 | 0.7019 | 0.8866 | 1.0000 | 0.8628 |
| 3 | [0.008, 0.09] | 0 | 1 | 50 | 0 | 0.7251 | 0.8607 | 1.0000 | 0.8619 |
| 3 | [0.008, 0.09] | 0 | 1 | 30 | 0 | 0.7119 | 0.8728 | 1.0000 | 0.8616 |
| 3 | [0.008, 0.09] | 0 | 0 | 20 | 0 | 0.6967 | 0.8815 | 1.0000 | 0.8594 |
| 3 | [0.008, 0.09] | 0 | 0 | 50 | 0 | 0.7094 | 0.8576 | 1.0000 | 0.8557 |
| 2 | [0.008, 0.09] | 0 | 1 | 20 | 0 | 0.6974 | 0.8695 | 1.0000 | 0.8556 |
| 2 | [0.008, 0.09] | 0 | 1 | 50 | 0 | 0.7143 | 0.8464 | 1.0000 | 0.8536 |
| 1 | [0.008, 0.09] | 0 | 1 | 20 | 0 | 0.6871 | 0.8704 | 1.0000 | 0.8525 |
| 1 | [0.008, 0.09] | 0 | 1 | 50 | 0 | 0.7028 | 0.8502 | 1.0000 | 0.8510 |

#### New Top 10

| POLY | BP | SIMULT | MOT24 | N_ACOMP | N_AROMA | Validity | Quality | Sens | Mean_QC |
|------|----|--------|-------|---------|---------|----------|---------|------|---------|
| 3 | [0.008, 0.09] | 0 | 1 | 20 | Inf | 0.7133 | 0.9120 | 1.0000 | 0.8751 |
| 2 | [0.008, Inf] | 1 | 1 | 50 | Inf | 0.7142 | 0.9025 | 1.0000 | 0.8722 |
| 1 | [0.008, Inf] | 1 | 1 | 50 | Inf | 0.7142 | 0.9014 | 1.0000 | 0.8719 |
| 3 | [0.008, 0.09] | 0 | 1 | 15 | Inf | 0.7043 | 0.9098 | 1.0000 | 0.8714 |
| 3 | [0.008, Inf] | 1 | 1 | 50 | Inf | 0.7115 | 0.9013 | 1.0000 | 0.8710 |
| 2 | [0.008, Inf] | 1 | 1 | 20 | Inf | 0.7028 | 0.9084 | 1.0000 | 0.8704 |
| 3 | [0.008, Inf] | 1 | 1 | 20 | Inf | 0.7035 | 0.9066 | 1.0000 | 0.8700 |
| 1 | [0.008, Inf] | 1 | 1 | 20 | Inf | 0.7017 | 0.9078 | 1.0000 | 0.8698 |
| 3 | [0.008, 0.09] | 0 | 1 | 30 | Inf | 0.7103 | 0.8970 | 1.0000 | 0.8691 |
| 3 | [0.008, 0.09] | 0 | 1 | 50 | Inf | 0.7183 | 0.8885 | 1.0000 | 0.8689 |

