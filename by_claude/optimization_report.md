# Denoising Pipeline Optimization Report (with RM_GMR / GSR)

- Generated: 2026-04-24 15:37:24
- Project:   `/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/no_PM_260423.mat`
- Subjects:  20
- Iterations: 216 (Stage 1: 192, Stage 2: 24, Prior: 0)

## New Parameter: RM_GMR

Binary flag. RM_GMR=1 adds a `Grey Matter` entry to `Preproc.confounds` (dim=1, deriv=0, power=1), equivalent to GSR. Motivated by Ciric 2017 / Mascali 2021 benchmarks.

## Validation

| Score | GUI | Re-computed | Diff |
|---|---|---|---|
| Validity    | 0.7133 | 0.7133 | -0.0000 |
| Quality     | 0.9120 | 0.9120 | +0.0000 |
| Sensitivity | 1.0000 | 1.0000 | -0.0000 |

## Best Parameter Combination

| Parameter | Value |
|---|---|
| POLY_ORD  | 3 |
| BP_HZ     | [0.008, 0.09] |
| SIMULT    | 0 |
| MOT24     | 1 |
| N_ACOMP   | 5 |
| N_AROMA   | 10 |
| **RM_GMR** | **1** |
| Validity  | 0.9768 |
| Quality   | 0.9694 |
| Sensitivity | 1.0000 |
| **Mean_QC** | **0.9821** |

## Score Distributions (main-effect summary)

### POLY_ORD

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|---|---|---|---|---|---|
| 2 | 0.8521 | 0.7743 | 0.8444 | 0.9789 | 96 |
| 3 | 0.8637 | 0.7850 | 0.8631 | 0.9751 | 120 |

### BP_HZ

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|---|---|---|---|---|---|
| [0.008 0.09] | 0.8460 | 0.7667 | 0.8783 | 0.9514 | 120 |
| [0.008 Inf] | 0.8742 | 0.7972 | 0.8254 | 1.0000 | 96 |

### SIMULT

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|---|---|---|---|---|---|
| 1 | 0.8204 | 0.7587 | 0.8363 | 0.9331 | 96 |
| 0 | 0.8890 | 0.7975 | 0.8696 | 1.0000 | 120 |

### MOT24

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|---|---|---|---|---|---|
| 1 | 0.8748 | 0.7853 | 0.8865 | 0.9903 | 120 |
| 0 | 0.8382 | 0.7739 | 0.8152 | 0.9606 | 96 |

### N_ACOMP

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|---|---|---|---|---|---|
| 5 | 0.8730 | 0.7525 | 0.8985 | 0.9739 | 68 |
| 10 | 0.9117 | 0.8210 | 0.9141 | 1.0000 | 8 |
| 15 | 0.9168 | 0.8310 | 0.9194 | 1.0000 | 8 |
| 20 | 0.8722 | 0.8254 | 0.8555 | 0.9535 | 68 |
| 50 | 0.8148 | 0.7504 | 0.7920 | 1.0000 | 64 |

### N_AROMA

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|---|---|---|---|---|---|
| 0 | 0.8455 | 0.7745 | 0.8366 | 0.9594 | 100 |
| 5 | 0.9091 | 0.8050 | 0.9225 | 1.0000 | 8 |
| 10 | 0.9119 | 0.8113 | 0.9245 | 1.0000 | 8 |
| Inf | 0.8633 | 0.7816 | 0.8619 | 0.9912 | 100 |

### RM_GMR

| Level | Mean_QC | Validity | Quality | Sensitivity | n |
|---|---|---|---|---|---|
| 0 | 0.8070 | 0.6345 | 0.8505 | 0.9792 | 108 |
| 1 | 0.9101 | 0.9260 | 0.8591 | 0.9743 | 108 |

## RM_GMR paired comparison

- Paired cells: 108
- Mean Mean_QC lift (GMR=1 minus GMR=0): +0.1031 (std 0.0373)
- % of cells where GMR=1 is better: 98.1%

## Top 10 (RM_GMR=1 arm)

| POLY | BP | SIMULT | MOT24 | N_ACOMP | N_AROMA | Validity | Quality | Sens | Mean_QC |
|---|---|---|---|---|---|---|---|---|---|
| 3 | [0.008, 0.09] | 0 | 1 | 5 | 10 | 0.9768 | 0.9694 | 1.0000 | 0.9821 |
| 3 | [0.008, 0.09] | 0 | 1 | 5 | Inf | 0.9767 | 0.9626 | 1.0000 | 0.9798 |
| 2 | [0.008, 0.09] | 0 | 1 | 5 | Inf | 0.9831 | 0.9562 | 1.0000 | 0.9798 |
| 3 | [0.008, 0.09] | 1 | 1 | 20 | 0 | 0.9813 | 0.9781 | NaN | 0.9797 |
| 3 | [0.008, 0.09] | 0 | 1 | 5 | 5 | 0.9717 | 0.9640 | 1.0000 | 0.9786 |
| 3 | [0.008, 0.09] | 0 | 1 | 10 | 10 | 0.9670 | 0.9658 | 1.0000 | 0.9776 |
| 3 | [0.008, 0.09] | 1 | 0 | 20 | Inf | 0.9773 | 0.9752 | NaN | 0.9763 |
| 3 | [0.008, 0.09] | 0 | 1 | 10 | Inf | 0.9715 | 0.9556 | 1.0000 | 0.9757 |
| 3 | [0.008, 0.09] | 0 | 1 | 15 | 10 | 0.9670 | 0.9591 | 1.0000 | 0.9754 |
| 3 | [0.008, 0.09] | 0 | 1 | 10 | 5 | 0.9610 | 0.9614 | 1.0000 | 0.9741 |

## Top 10 (RM_GMR=0 arm, for comparison)

| POLY | BP | SIMULT | MOT24 | N_ACOMP | N_AROMA | Validity | Quality | Sens | Mean_QC |
|---|---|---|---|---|---|---|---|---|---|
| 3 | [0.008, 0.09] | 0 | 1 | 20 | Inf | 0.7133 | 0.9120 | 1.0000 | 0.8751 |
| 2 | [0.008, Inf] | 1 | 1 | 50 | Inf | 0.7142 | 0.9025 | 1.0000 | 0.8722 |
| 3 | [0.008, 0.09] | 0 | 1 | 15 | Inf | 0.7043 | 0.9098 | 1.0000 | 0.8714 |
| 3 | [0.008, Inf] | 1 | 1 | 50 | Inf | 0.7115 | 0.9013 | 1.0000 | 0.8710 |
| 2 | [0.008, Inf] | 1 | 1 | 20 | Inf | 0.7028 | 0.9084 | 1.0000 | 0.8704 |
| 3 | [0.008, Inf] | 1 | 1 | 20 | Inf | 0.7035 | 0.9066 | 1.0000 | 0.8700 |
| 3 | [0.008, 0.09] | 0 | 1 | 50 | Inf | 0.7183 | 0.8885 | 1.0000 | 0.8689 |
| 3 | [0.008, 0.09] | 0 | 1 | 20 | 0 | 0.7083 | 0.8857 | 1.0000 | 0.8647 |
| 3 | [0.008, 0.09] | 0 | 1 | 20 | 5 | 0.7064 | 0.8861 | 1.0000 | 0.8642 |
| 2 | [0.008, 0.09] | 0 | 1 | 20 | Inf | 0.7009 | 0.8901 | 1.0000 | 0.8637 |

## Caveats

- 0 of 216 new combinations failed (0.0%).
- Mean_QC = simple mean of Validity, Quality, Sensitivity.
- RM_GMR implementation: adds Grey Matter confound (dim=1, deriv=0, power=1).
