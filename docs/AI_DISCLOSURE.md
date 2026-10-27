# AI-assistance disclosure

Code in this repository was developed with **AI assistance under the author's direction**:

- The killifish optimisation stage was drafted with a science harness (ClawsGO);
- the R analysis pipeline, the validation layer (fold-safe CV, FDR layer, tests), and the
  figures were developed with a local coding agent, with the author specifying and reviewing
  every scientific decision;
- the final report was assembled by the author from verified outputs; no reported number is
  accepted without an executable path to its producing script and artifact.

**What this means for reuse.** The code should be read as "AI-assisted implementation of the
author's analysis design", not as hand-written reference software. The author verified
numerical behaviour by re-running the pipeline in a sandbox, comparing outputs against a frozen
baseline, and regression-testing the fixes.

**Excluded upstream material.** Two unused files from the original working directory
(`_tpls_source.R`, `_tplsda_source.R`) were verbatim deparse copies of upstream `tensorOmics`
functions. They are excluded from this release; the package itself is used as an installed
dependency (AGPL-3.0+, not vendored). See `NOTICE`.

**Provenance statement.** Third-party methodological credit: `tensorOmics` (Brendan Lu) for
the tensor methods; `mixOmics` for matrix comparators; `sva`/`PLSDAbatch` for batch correction;
GEO datasets and the Li et al. (2017) supplement for the data.
