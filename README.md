# TER — Tensor reanalysis and validation of two longitudinal omics datasets

**A reproducible reanalysis of killifish lifespan transcriptomics (GSE150318) and human
vaccine-response transcriptomics (GSE79396), with an explicit statistical-validation layer
for tensor methods on longitudinal data.**

This repository contains the analysis code, the final report (Markdown / LaTeX / PDF / DOCX),
figures, and reproduction records for the project. It was developed as a research portfolio
project and is published as a *release candidate* — no external data are redistributed here.

---

## What the project found

**Killifish (GSE150318; 113 fish, fin transcriptomes at 10 and 20 weeks, age at death).**
- Early transcriptome **levels** carry a reproducible association with lifespan: in-sample
  R² ≈ 0.5 (bootstrap 95% CI [0.43, 0.66]); subject-level label permutation p = 0.002.
- A structure-preserving **time reshuffle** gives p = 0.150: the signal is level-based rather
  than change-based (levels-only R² 0.541 vs 0.488 for 20-week-only and 0.329 for 10-week-only).
- The signal is organised in a **circadian-clock module** (7/8 core clock genes FDR-significant,
  direction-consistent) and a **proteasome/ubiquitin** axis; the original study's mitochondrial
  complex-I mechanism is recovered only as a diffuse trend (1.60× enrichment).
- **No prediction claim**: grouped LOSO R² ≈ 0.31 for the audit model; a ridge baseline matches
  it; a two-component tensor regression does not generalise (negative LOSO R²).

**Vaccine (GSE79396; 72 subjects, PBMC transcriptomes at days 0/1/3/7).**
- Two classical matrix-analysis failure modes are quantified: **pseudo-replication**
  (matrix PLS-DA permutation p = 0.002 with naive row-level shuffling vs 0.094 with
  subject-level shuffling) and **cross-validation leakage** (row-level LOSO balanced accuracy
  0.722 vs 0.644 subject-grouped; the grouped tensor model reaches 0.685).
- The **processing plate** (recovered from raw CEL file names) is the finest batch level and is
  confounded with age. After plate-level ComBat the age association collapses (permutation
  p ≈ 0.003 → 0.177), and the original "inflammaging" panel loses significance under
  plate- and sex-adjusted models. The original axis is not supported.
- PLSDA-batch is applied **on top of** ComBat: a layered correction, not an independent
  replication.

**Methodological contribution.** A transferable validation protocol — subject-level permutation,
structure-preserving time reshuffle, subject-grouped and fold-safe cross-validation, family-wise
FDR, and nested batch auditing — together with a machine-precision-verified manual projection
for tensor models (needed because `tensorOmics` has no `predict()` method).

Full details, tables and limitations: [`report/TER_report_final.md`](report/TER_report_final.md)
(also `.tex` / `.pdf` / `.docx`).

---

## Repository layout

```
code/
  fish_pipeline/    Python: dataset build + main killifish analysis (01–06)
  pipeline/         R: preprocessing, tensor analyses, batch audits, FDR, figures (07–41, 97–99)
  R/                helpers: lib_tensor_helpers.R (projection etc.), lib_p2a_helpers.R
data/               acquisition instructions + small metadata files (no raw data)
report/             final report (md/tex/pdf/docx) + figures/
docs/               running guide, AI disclosure, contributions, validation report
reproduction/       run summaries and comparison records from the final reproduction
```

## Quick start

1. **Get the data** — follow [`data/README.md`](data/README.md) (GEO GSE150318 and GSE79396;
   PDB/GEO terms apply; Li et al. 2017 supplement for the metabolomics block).
2. **Environment** — R 4.6 with `tensorOmics` (pinned commit `61afa63`), `mixOmics`, `MASS`,
   `data.table`, `lme4`, `lmerTest`, `sva`, `PLSDAbatch`, `ggplot2`, `patchwork`;
   Python 3.14 with `numpy`, `pandas`, `scipy`, `scikit-learn`.
3. **Run** — execute from the repository root (see [`docs/RUNNING.md`](docs/RUNNING.md) for the
   exact order and per-script inputs/outputs):

```bash
python code/fish_pipeline/02_build_fish_dataset.py
python code/fish_pipeline/03_fish_analysis.py
python code/fish_pipeline/04_annotate_hits.py
Rscript code/pipeline/07_fish_comparison.R
Rscript code/pipeline/07c_fish_tplsda_loo.R
Rscript code/pipeline/08_fish_tpls_regression.R
Rscript code/pipeline/10_prepare_gse79396.R
Rscript code/pipeline/11_reproduce_v1.R
Rscript code/pipeline/20_vaccine_valueprop.R
Rscript code/pipeline/21_vaccine_batch.R
Rscript code/pipeline/21b_vaccine_heterogeneity.R
Rscript code/pipeline/21c_vaccine_chip_batch.R
Rscript code/pipeline/21d_vaccine_chip_adjust.R
Rscript code/pipeline/22_vaccine_deep.R
Rscript code/pipeline/22b_baseline_dynamics.R
Rscript code/pipeline/30_vaccine_multiblock.R
Rscript code/pipeline/30b_multiblock_corrected_align.R
Rscript code/pipeline/25_fdr_corrections.R
Rscript code/pipeline/40_vaccine_figures.R
Rscript code/pipeline/41_fish_figures.R
```

## Reproducibility

- Fixed seeds and recorded environments; every reported number traces to a script and an output
  artifact (see `reproduction/` and the report's Appendix traceability table).
- Analyses were re-run in a sandboxed mirror of the full project; outputs were compared against
  the frozen baseline at value level. Where preprocessing choices could leak information across
  folds (feature selection, scaling), fold-safe re-derivations were run and are reported.
- Engineering fixes verified by regression tests: idempotent/atomic table updates, GSM-keyed
  plate mapping (reproduces the stored mapping exactly), deterministic diagnostic scripts.

## Data and licensing

- Raw data are **not** redistributed; `data/README.md` lists sources, and `data/CHECKSUMS.sha256`
  the expected file hashes.
- Our code: **MIT** (see `LICENSE`). Third-party notices and attributions: `NOTICE`.
- `tensorOmics` (Brendan Lu, AGPL-3.0+) is used as an installed R package; it is not vendored here.

## AI-assistance disclosure

Code in this project was developed with AI assistance (a science harness for the killifish
optimisation stage; a local coding agent for the R pipeline, validation layer, tests and
figures), under the author's direction; scientific decisions and verification are the author's.
Two unused files that were verbatim copies of upstream `tensorOmics` source were excluded from
this release. See [`docs/AI_DISCLOSURE.md`](docs/AI_DISCLOSURE.md) and
[`docs/CONTRIBUTIONS.md`](docs/CONTRIBUTIONS.md).

## References

1. tensorOmics framework: Lu, B. et al. — *tensorOmics* R package and preprint
   (github.com/brendanlu/tensorOmics; pinned commit `61afa63`).
2. GSE150318: Baumgart et al. (killifish lifespan transcriptomics).
3. GSE79396: Li et al. (2017) — vaccine response transcriptomics + supplements.
4. Mor et al. (2022) — tensor PCA for longitudinal data.
5. Kilmer et al. (2021) — tensor SVD (t-SVDM).
