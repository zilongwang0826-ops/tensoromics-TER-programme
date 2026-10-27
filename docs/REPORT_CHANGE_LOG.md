# Report change log — P2B candidate → final (release candidate)

Baseline: `audit/phase4-p2b/revised/TER_report_v2_p2b.md`.
Final: `repo/report/TER_report_final.md` (plus .tex/.pdf/.docx).
All changes follow from the completed reproduction; old→new numeric comparisons are recorded in
`VALIDATION_REPORT.md` §4b.

| # | Location | Old | New | Reason |
|---|---|---|---|---|
| 1 | §3.3 item 4 | "…unchanged at full precision." | Rewritten: pooled vaccine LOSO unchanged; killifish baseline R² shifts ≤ 0.001; within-cohort matrix-arm LOSO changes (0.87→0.78; 0.38→0.52); tensor arms unchanged; fold-safe values are cited | Fold-safe re-derivations completed (verified runs) |
| 2 | Table 1 | tPLS regression LOSO R² = −0.086 | −0.083 (labelled fold-safe) | Fold-safe re-derivation (08c) |
| 3 | Table 1 | Audit LOSO R²/AUC = 0.310/0.792 | 0.308/0.792 (labelled fold-safe) | Fold-safe re-derivation (stage) |
| 4 | Table 1 | Ridge LOSO R²/AUC = 0.381/0.797 | 0.380/0.797 (labelled fold-safe) | Fold-safe re-derivation (stage) |
| 5 | §4.2 | "…does not generalise at all (LOSO R² = −0.086)" | −0.083 | Same as #2 |
| 6 | §5.3 | "LOSO balanced accuracy 0.79–0.87" / "LOSO 0.38–0.41" | "fold-safe LOSO balanced accuracy 0.78–0.79" / "fold-safe LOSO 0.41–0.52" | Within-cohort fold-safe re-derivation |
| 7 | §7 | "…(LOSO R² = 0.381)" | 0.380 | Same as #4 |
| 8 | §7.6 | "…excluded from correction claims until the final re-execution stores them. Two audit-layer numbers are excluded… pending the final re-run…" | Final statement: BH leaves reported conclusions unchanged (exploratory joint families can shift uncited calls); fold-safe re-derivations completed; deviations documented in the validation record | MB06 persisted + all fold-safe runs completed |

**Unchanged by these edits** (verified): Abstract; pseudo-replication and leakage numbers
(0.002/0.094; 0.722/0.644/0.685); fish classification LOSO (0.620/0.682/0.660); batch ladder
values; F22 removal; layered-correction wording; FA5 caption note; all figures (FA1/FA3/FA5 are
the corrected P1 versions; all others regenerated identical).
