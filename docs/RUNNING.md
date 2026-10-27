# Running guide (order, inputs, outputs)

Set the repository root once; all scripts resolve paths relative to it:

```bash
export TER_ROOT=/path/to/repo          # Linux/macOS
setx TER_ROOT C:\path\to\repo          # Windows PowerShell
```

(If unset, scripts fall back to the current working directory — run from the repo root.)

Optional path overrides:

- `GSE79396_SERIES_MATRIX` - GSE79396 series matrix; defaults to
  `data/GSE79396/GSE79396_series_matrix.txt.gz`.
- `TER_V1_DIR` - directory with the v1 result tables (`matrix_pca_scores.csv`,
  `tensor_pca_scores.csv`, `top_loadings_pc*.csv`); required only by `11_reproduce_v1.R`.
- `TER_V1_METADATA` - optional v1 sample metadata for the comparison block in
  `10_prepare_gse79396.R`; defaults to `data/v1_tables/sample_metadata_complete.csv`.

## Fish pipeline (Python) — `code/fish_pipeline/`

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `01_parse_gff.py` | `data/GSE150318/*.gff3.gz` | gene annotation tables |
| 2 | `02_build_fish_dataset.py` | counts + series matrix | `fish_data.npz`, `fish_metadata.tsv` |
| 3 | `03_fish_analysis.py` | `fish_data.npz` | `results/…/TA_summary.json`, permutations, LOSO, FDR, sensitivity |
| 4 | `04_annotate_hits.py` | TA tables | annotated top tables, enrichment summary |
| 5 | `05_fish_figures_gsea.py` | TA tables (+ network for Enrichr) | FA5 input tables/figures |
| 6 | `06_ablation.py` | `fish_data.npz` | ablation table |

## Fish comparison (R) — `code/pipeline/`

`07_fish_comparison.R` → `fish_manova.csv`, `fish_perm_inflation.csv`,
`fish_loso_naive_vs_grouped.csv`, scores;
`07b_fish_scores_fix.R` → score tables; `07c_fish_tplsda_loo.R` → tensor LOSO row;
`07d_fish_foldsafe.R` → `fish_foldsafe.csv` (fold-safe re-derivation; requires the full-gene
export produced by `tools/export_fish_full.py` in the reproduction workflow);
`08_fish_tpls_regression.R`, `08b_fish_tpls_reg_check.R` → tPLS regression + LOSO.

## Vaccine pipeline (R) — `code/pipeline/`

Run in this order (dependencies in parentheses):

```
10_prepare_gse79396.R            → data/GSE79396/*.rds, samples_all.csv
11_reproduce_v1.R                → preprocessed_A.rds (selection, z-scores, tensors)
20_vaccine_valueprop.R           → manova / permutation / LOSO / scores      (needs 11)
21_vaccine_batch.R               → batch ladder metrics                       (needs 11)
21b_vaccine_heterogeneity.R      → within-cohort summaries                    (needs 11)
21b2_comps_full.R                → full per-component p values                (needs 11)
21c_vaccine_chip_batch.R         → plate serials + plate correction           (needs 11; before 21d/22/22b/30/30b)
21d_vaccine_chip_adjust.R        → covariate-adjusted age stats               (needs 21c plate column)
22_vaccine_deep.R                → v1 genes, top genes, BTM                     (needs 21c)
22b_baseline_dynamics.R          → day-resolved axis model                      (needs 22)
30_vaccine_multiblock.R          → multi-block integration (incl. metabolome audit) (needs 11)
30b_multiblock_corrected_align.R → plate-corrected alignment (idempotent)      (needs 30, 21c)
25_fdr_corrections.R             → family-wise BH corrections over the stored p values
40_vaccine_figures.R             → FV1–FV4 figures
41_fish_figures.R                → FF1–FF2, FA1–FA5 figures
```

## Diagnostic scripts

`97_prep_checks.R`, `98_test_vaccine_pred.R`, `99_test_predict.R` — smoke/diagnostic checks
(they must exit 0). `14_chip_batch_check.py` — plate crosstabs (path configurable).

## Determinism notes

- All permutation/bootstrap heavy scripts set explicit seeds; R scripts are deterministic given
  seeds and package versions (recorded in the report).
- Table writers use idempotent/atomic updates where re-runs previously duplicated rows
  (`30b`, `07c`) or rewrote in place (`21c`).
