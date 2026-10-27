#!/usr/bin/env python3
"""Task A — Step 3: post-hoc annotation of FDR hits + complex-I enrichment vs background."""
import os
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact, binomtest

DATA = Path("os.path.join(os.environ.get("TER_ROOT", "."), "data", "GSE150318")")
OUTT = Path("os.path.join(os.environ.get("TER_ROOT", "."), "results", "tables")")

ann = pd.read_csv(DATA / "gene_annotation.tsv", sep="\t").set_index("gene_id")
fdr = pd.read_csv(OUTT / "TA_gene_fdr.csv.gz", index_col=0)
# the saved table is in the same gene order as fish_data.npz; attach ids by position
d = np.load(DATA / "fish_data.npz", allow_pickle=True)
assert len(d["genes"]) == len(fdr), "gene order mismatch"
fdr.index = [str(g) for g in d["genes"]]
fdr.index.name = "gene_id"
print("FDR table:", fdr.shape, "| q<0.05:", int((fdr["delta_age_q"] < 0.05).sum()))

# annotated top-50 FDR genes
top = fdr.sort_values("delta_age_p").head(50).copy()
top = top.join(ann[["short_name", "description", "orthologs"]])
top.round(5).to_csv(OUTT / "TA_topFDR_annotated.csv")
print("\nTop-15 FDR genes (annotated):")
print(top.head(15)[["delta_age_rho", "delta_age_q", "short_name", "description"]].to_string())

# where do the report's four bootstrap-stable genes land?
four = ["Nfu_g_1_010933", "Nfu_g_1_007097", "Nfu_g_1_001084", "Nfu_g_1_022498"]
fdr2 = fdr.copy()
fdr2["p_rank"] = fdr2["delta_age_p"].rank().astype(int)
print("\nreport's four bootstrap-stable genes, ranked by delta-age p:")
for g in four:
    if g in fdr2.index:
        r = fdr2.loc[g]
        print(f"  {g} | rho={r['delta_age_rho']:+.3f} q={r['delta_age_q']:.2e} "
              f"p-rank={int(r['p_rank'])}/26739 | {ann.loc[g,'short_name']} | {str(ann.loc[g,'description'])[:60]}")

# circadian candidates
circ = ["PER1", "NPAS2", "NR1D1", "PER2", "PER3", "CLOCK", "ARNTL", "CRY1", "CRY2", "DBP", "TEF", "HLF", "NFIL3"]
circ_ids = ann[ann["short_name"].astype(str).str.upper().isin([c for c in circ])]
rows = []
for g, r in circ_ids.iterrows():
    if g in fdr2.index:
        rr = fdr2.loc[g]
        rows.append({"gene_id": g, "short_name": r["short_name"], "rho": rr["delta_age_rho"],
                     "q": rr["delta_age_q"], "p_rank": int(rr["p_rank"]), "in_top50_fdr": int(rr["p_rank"]) <= 50})
circ_df = pd.DataFrame(rows).sort_values("p_rank")
circ_df.round(6).to_csv(OUTT / "TA_circadian_genes.csv", index=False)
print("\ncircadian clock genes in the fish data:")
print(circ_df.to_string())

# complex-I enrichment vs global background
ci = pd.read_csv(OUTT / "TA_complexI_members.csv")
n_bg = len(fdr); k_bg = int((fdr["delta_age_q"] < 0.05).sum())
n_ci = len(ci); k_ci = int((ci["delta_age_q"] < 0.05).sum())
print(f"\nbackground q<0.05 rate: {k_bg}/{n_bg} = {k_bg/n_bg:.3f}")
print(f"complex I q<0.05 rate:  {k_ci}/{n_ci} = {k_ci/n_ci:.3f}")
bt = binomtest(k_ci, n_ci, k_bg / n_bg, alternative="greater")
print(f"binomial test (vs background rate): p = {bt.pvalue:.4f}, enrichment = {(k_ci/n_ci)/(k_bg/n_bg):.2f}x")

# also: median |rho| of complex I vs all genes
all_abs = fdr["delta_age_rho"].abs()
ci_abs = np.abs(ci["delta_age_rho"])
print(f"median |rho|: complex I = {np.nanmedian(ci_abs):.3f} | all genes = {all_abs.median():.3f}")
from scipy.stats import mannwhitneyu
u = mannwhitneyu(ci_abs.dropna(), all_abs, alternative="greater")
print(f"Mann-Whitney (complex I |rho| > background): p = {u.pvalue:.4f}")

# save a compact summary for the report
summary = {
    "n_genes_total": int(n_bg), "n_fdr_sig": int(k_bg), "frac_sig": float(k_bg / n_bg),
    "complexI_n": int(n_ci), "complexI_sig": int(k_ci), "complexI_frac": float(k_ci / n_ci),
    "complexI_enrichment": float((k_ci / n_ci) / (k_bg / n_bg)),
    "complexI_binom_p": float(bt.pvalue),
    "complexI_median_absrho": float(np.nanmedian(ci_abs)),
    "all_median_absrho": float(all_abs.median()),
    "complexI_mw_p": float(u.pvalue),
}
import json
(OUTT / "TA_enrichment_summary.json").write_text(json.dumps(summary, indent=2))
print("\nsaved TA_enrichment_summary.json")
