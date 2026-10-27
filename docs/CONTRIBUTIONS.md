# Contributions

**What the author did (this project):**
- Defined the reanalysis questions and the validation protocol (subject-level permutation,
  structure-preserving time reshuffle, subject-grouped cross-validation, fold-safe
  re-derivation, family-wise FDR, nested batch auditing).
- Executed and verified the analyses on the two public datasets; built the figures and the
  final report; audited all reported numbers end-to-end.
- Implemented the verification utilities: the validated manual projection for tensor models
  (required because `tensorOmics` has no `predict()`), the fold-safe cross-validation, the FDR
  correction layer, idempotent/safe table updates, and the reproduction/regression harnesses.

**What is reused with attribution:**
- `tensorOmics` (Brendan Lu; AGPL-3.0+; pinned commit `61afa63`) for the tensor methods
  (tubal-transform algebra, mean-deviation centring, t-SVDM, tPCA/tPLS/tPLS-DA/block.tPLS) —
  used as an installed dependency, not vendored.
- Matrix comparators via `mixOmics`; batch correction via `sva` and `PLSDAbatch`.

**AI assistance:** disclosed separately in `AI_DISCLOSURE.md` (AI-assisted implementation under
the author's direction; scientific decisions and verification by the author). Eight specific
decision-ownership items from the project's audit remain to be confirmed individually before
this statement is used in applications.
