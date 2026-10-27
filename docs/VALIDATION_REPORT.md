# Final validation report

Release-candidate validation for the TER reproduction. All runs were executed in a sandboxed mirror (`runs/mirror`) of the frozen baseline project; the baseline itself was never modified (hash-verified).

## 1. Steps executed

### run_summary

| step | status | seconds |
|---|---|---|
| fish01 | TIMEOUT | 0,5 |
| fish02 | TIMEOUT | 0 |
| fish03 | TIMEOUT | 0 |
| fish04 | TIMEOUT | 0 |
| fish06 | TIMEOUT | 0 |
| fish_export_full | TIMEOUT | 0 |
| fish07 | EXIT_ | 340,9 |
| fish07b | EXIT_ | 4,6 |
| fish07c | EXIT_ | 134,7 |
| fish07d | EXIT_ | 2,9 |
| fish08 | EXIT_ | 251 |
| fish08b | EXIT_ | 6,5 |
| fish_stage_foldsafe | TIMEOUT | 0 |
| vacc10 | EXIT_ | 47,8 |
| vacc11 | EXIT_ | 10,5 |
| vacc20 | EXIT_ | 1032,5 |
| vacc21 | EXIT_ | 1114,5 |
| vacc21b | EXIT_ | 4,8 |
| vacc21b2 | EXIT_ | 4,4 |
| vacc21c | EXIT_ | 9,6 |
| vacc21d | EXIT_ | 1,4 |
| vacc22 | EXIT_ | 8,4 |
| vacc22b | EXIT_ | 6 |
| vacc30 | EXIT_ | 1,1 |
| vacc30b | EXIT_ | 4,9 |
| vacc25fdr | EXIT_ | 0,8 |
| fig40 | EXIT_ | 9,8 |
| fig41 | EXIT_ | 5,2 |

### run_summary_part2

| step | status | seconds |
|---|---|---|
| fish01 | EXIT_ | 0,3 |
| fish02 | EXIT_ | 132,1 |
| fish03 | EXIT_ | 73,4 |
| fish04 | EXIT_ | 2,8 |
| fish06 | EXIT_ | 2,9 |
| fish_export_full | EXIT_ | 21,8 |
| fish07d | EXIT_ | 448,4 |
| fish_stage_foldsafe | EXIT_ | 33,8 |
| vacc_within_foldsafe | EXIT_ | 293,2 |
| fig41_rerun | EXIT_ | 4,6 |
| compare | EXIT_ | 1,4 |

### run_summary_part3a

| step | status | seconds |
|---|---|---|
| r21b2 | DONE(no-code) | 4,9 |
| r21c | DONE(no-code) | 634,3 |
| r21d | DONE(no-code) | 1,7 |
| r22 | DONE(no-code) | 19,3 |
| r22b | DONE(no-code) | 6,3 |
| f08b | DONE(no-code) | 5,8 |

### run_summary_part3b

| step | status | seconds |
|---|---|---|
| r30 | DONE(no-code) | 644,9 |
| r30b | DONE(no-code) | 4,7 |
| r25fdr | DONE(no-code) | 0,7 |
| rfig40 | DONE(no-code) | 8 |
| compare_final | DONE(no-code) | 1,8 |
| validation_report | DONE(no-code) | 1,3 |

### run_summary_part3c

| step | status | seconds |
|---|---|---|
| g01 | EXIT_ | 3,3 |
| g06 | EXIT_ | 3 |
| g08b | EXIT_ | 277 |
| g08c | EXIT_ | 306 |
| gfig41 | EXIT_ | 4,4 |

## 2. Mirror vs frozen baseline (value-level comparison)

- files compared: **63**; identical: **61**; changed/missing: **2**

| group | file | status | max_abs_diff | char_mismatch | notes |
|---|---|---|---|---|---|
| fishStage | TA_gene_fdr.csv.gz | VALUES_EQUAL | 0.0 | 0.0 |  |
| fishStage | TA_gsea_enrichr.csv | MISSING_IN_MIRROR | nan | nan |  |

## 3. Fish fold-safe re-derivation (grouped LOSO)

| metric | committed (leaky) | fold-safe | delta |
|---|---|---|---|
| matrix PLS-DA grouped LOSO | 0.620056 | 0.620056 | +0.000000 |
| tPLS-DA grouped LOSO | 0.681576 | 0.681576 | +0.000000 |

### Stage-model LOSO (audit / PCA / PLS-age / ridge)

| model | foldsafe_loso_r2_age | foldsafe_loso_auc_group | committed_loso_r2_age | committed_loso_auc_group | delta_r2 | delta_auc |
|---|---|---|---|---|---|---|
| audit | 0.308388 | 0.792216 | 0.309697 | 0.792216 | -0.001309 | 0.000000 |
| pca | 0.291660 | 0.782800 | 0.292781 | 0.782800 | -0.001121 | 0.000000 |
| pls_age | 0.327837 | 0.779975 | 0.328326 | 0.780917 | -0.000489 | -0.000942 |
| ridge | 0.380328 | 0.796924 | 0.380886 | 0.796610 | -0.000558 | 0.000314 |

### tPLS-regression LOSO (fold-safe)

| metric | leaky | foldsafe | delta | mean_sel_overlap |
|---|---|---|---|---|
| tPLS regression LOSO R2 | -0.085510 | -0.083423 | 0.002087 | 991.035398 |

## 4. Within-cohort LOSO fold-safe check

| cohort | n | leaky_matrix_acc | foldsafe_matrix_acc | leaky_tensor_acc | foldsafe_tensor_acc | delta_matrix | delta_tensor |
|---|---|---|---|---|---|---|---|
| Atlanta | 37 | 0.379464 | 0.519345 | 0.410714 | 0.410714 | 0.139881 | 0 |
| Denver | 35 | 0.866667 | 0.775000 | 0.791667 | 0.791667 | -0.091667 | 0 |

## 4b. Old vs new (fold-safe) values — consolidated

| metric | old (committed) | new (fold-safe) | delta |
|---|---|---|---|
| matrix PLS-DA grouped LOSO | 0.620056 | 0.620056 | +0.000000 |
| tPLS-DA grouped LOSO | 0.681576 | 0.681576 | +0.000000 |
| tPLS regression LOSO R2 | -0.085510 | -0.083423 | +0.002087 |
| stage audit LOSO R2 | 0.309697 | 0.308388 | -0.001309 |
| stage pca LOSO R2 | 0.292781 | 0.291660 | -0.001121 |
| stage pls_age LOSO R2 | 0.328326 | 0.327837 | -0.000489 |
| stage ridge LOSO R2 | 0.380886 | 0.380328 | -0.000558 |
| within-cohort Atlanta matrix LOSO | 0.379464 | 0.519345 | +0.139881 |
| within-cohort Atlanta tensor LOSO | 0.410714 | 0.410714 | +0.000000 |
| within-cohort Denver matrix LOSO | 0.866667 | 0.775000 | -0.091667 |
| within-cohort Denver tensor LOSO | 0.791667 | 0.791667 | +0.000000 |

Report-text updates driven by these values are listed in the report's change log (Table 1 and Section 5.3 cite the fold-safe versions where they differ).

## 5. FDR corrections

- families: 15; variants: 36; variants where BH changes a 0.05-decision vs raw: 1

## 6. Known deviations and limits

- V02 sandbox run: 11 min 25 s (above the 6–10 min expectation; disclosed in the phase record).
- Fold-safe re-derivations changed four reported quantities (all updated in the final report; old→new in section 4b): audit LOSO R² 0.310→0.308; ridge LOSO R² 0.381→0.380; tPLS-regression LOSO R² −0.086→−0.083; within-cohort matrix-arm LOSO Denver 0.87→0.78 and Atlanta 0.38→0.52 (tensor arms unchanged). No conclusion changed.
- MB06 metabolome audit p-values are now persisted (`vacc_metab_audit.csv`) and included as an FDR reference family (final: 15 families / 36 variants; family corrections leave every reported conclusion unchanged; one exploratory joint reference family flips two uncited gene-level calls at the 0.05 threshold).
- `TA_gsea_enrichr.csv` (script 05) was not regenerated (requires network access); enrichment figures remain the baseline-derived (P1) versions.
- Assembly fixes found and resolved during the reproduction: 21c GSM quote-stripping (affected 21d/22/22b/30/30b chain; re-run successfully), 21b2 validation column suffixes, and a comparison-tool dtype bug.
- `08b_tpls_reg_check.R` had never executed due to a filename typo in earlier runners; corrected and now verified (projection cor = 1.000000, max|diff| = 0; LOSO R² = −0.086 reproduced).
- Figures FA1/FA3/FA5 are the corrected (P1) versions; all other figures were regenerated by the pipeline and are byte-identical to the P1 set.
- Deviations disclosed at earlier phases (matrix-probability supplement runtime above its declared estimate; BH-checker v1 implementation bug, self-caught) remain recorded in the project audit trail.

## 7. Verdict

Core results are generated by the committed code; value-level comparison against the frozen baseline is recorded above; remaining non-blocking items are tracked in BACKLOG.md.
