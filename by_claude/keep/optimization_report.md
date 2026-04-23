# Denoising Pipeline Optimization Report

- Generated: 2026-04-18 21:45:25
- Project:   `/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/v3_no_param.mat`
- Subjects:  20
- Iterations: 75 (Stage 1: 72, Stage 2: 3)

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
| N_AROMA   | 0 |
| Validity  | 0.7083 |
| Quality   | 0.8857 |
| Sensitivity | 1.0000 |
| **Mean_QC** | **0.8647** |

## Score Distributions (main-effect summary)

Mean Mean_QC across combinations grouped by each parameter level.

### POLY_ORD

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 1 | 0.7981 | 0.6257 | 0.8258 | 0.9826 | 24 |
| 2 | 0.7961 | 0.6267 | 0.8312 | 0.9709 | 24 |
| 3 | 0.8019 | 0.6406 | 0.8399 | 0.9598 | 27 |

### BP_HZ

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| [0.008 0.09] | 0.7799 | 0.6262 | 0.8434 | 0.9353 | 39 |
| [0.008 Inf] | 0.8193 | 0.6371 | 0.8208 | 1.0000 | 36 |

### SIMULT

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 1 | 0.7690 | 0.6381 | 0.8096 | 0.9281 | 36 |
| 0 | 0.8264 | 0.6252 | 0.8538 | 1.0000 | 39 |

### MOT24

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 1 | 0.8139 | 0.6443 | 0.8522 | 0.9844 | 39 |
| 0 | 0.7825 | 0.6174 | 0.8114 | 0.9568 | 36 |

### N_ACOMP

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 5 | 0.8001 | 0.5260 | 0.8959 | 0.9785 | 24 |
| 10 | 0.8503 | 0.6761 | 0.8749 | 1.0000 | 1 |
| 15 | 0.8628 | 0.7019 | 0.8866 | 1.0000 | 1 |
| 20 | 0.8225 | 0.7092 | 0.8382 | 0.9321 | 24 |
| 30 | 0.8616 | 0.7119 | 0.8728 | 1.0000 | 1 |
| 50 | 0.7664 | 0.6509 | 0.7579 | 1.0000 | 24 |

### N_AROMA

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|-------|---------|----------|---------|-------------|---|
| 0 | 0.7988 | 0.6314 | 0.8326 | 0.9706 | 75 |

## Recommendation

Use the **Best Parameter Combination** above for this dataset.

## Caveats

- 0 of 75 combinations failed (0.0%).
- QC scores recomputed on CONN's sample-voxel network (~500-2000 grey-matter voxels).
- Two-stage grid: Stage 1 tested 72 coarse combos; Stage 2 3 fine combos around the winner.
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

