#!/usr/bin/env python3
"""Task A — Step 4: figures for the upgraded fish report + g:Profiler GSEA (zebrafish orthologs)."""
import os
import json
import urllib.request
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

DATA = Path("os.path.join(os.environ.get("TER_ROOT", "."), "data", "GSE150318")")
OUTF = Path("os.path.join(os.environ.get("TER_ROOT", "."), "results", "figures")")
OUTT = Path("os.path.join(os.environ.get("TER_ROOT", "."), "results", "tables")")
OUTF.mkdir(exist_ok=True, parents=True)

d = np.load(DATA / "fish_data.npz", allow_pickle=True)
genes = np.array([str(g) for g in d["genes"]])
st = np.load(DATA / "TA_model_state.npz", allow_pickle=True)
ages = d["ages"]; groups = d["groups"]
T = st["T"]; load_gene = st["load_gene"]; sel1000 = st["sel1000"]
rho_g = st["rho_g"]; q_g = st["q_g"]
perm = np.load(OUTT / "TA_permutations.npz")
cv = pd.read_csv(OUTT / "TA_loso_cv.csv")
ann = pd.read_csv(DATA / "gene_annotation.tsv", sep="\t").set_index("gene_id")
sel1000_list = [str(genes[i]) for i in sel1000]

# ---------- Fig 1: score plot + age regression ----------
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
col = {"0": "#4C72B0", "1": "#DD8452"}
for g, lab in [(0, "short-lived"), (1, "long-lived")]:
    sel = groups == g
    axes[0].scatter(T[sel, 0], T[sel, 1], s=28, alpha=0.85, color=col[str(g)], label=lab)
axes[0].set_xlabel("audit component 1"); axes[0].set_ylabel("audit component 2")
axes[0].legend(fontsize=9); axes[0].set_title("(a) trajectory-PLS scores (median-lifespan split)")
axes[0].axhline(0, color="grey", lw=0.4); axes[0].axvline(0, color="grey", lw=0.4)
axes[1].scatter(T[:, 0], ages, s=28, alpha=0.85, c=[col[str(g)] for g in groups])
r1 = np.corrcoef(T[:, 0], ages)[0, 1]
axes[1].set_xlabel("audit component 1"); axes[1].set_ylabel("age at death (days)")
axes[1].set_title(f"(b) component 1 vs lifespan (r = {r1:.3f})")
fig.tight_layout()
fig.savefig(OUTF / "TA_fig1_scores.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- Fig 2: permutation nulls + LOSO ----------
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
axes[0].hist(perm["R2_null_time"], bins=30, color="#B0B0B0", edgecolor="white")
axes[0].axvline(0.505, color="#C44E52", lw=2)
axes[0].set_title("(a) within-fish time reshuffle\n(499 perms; p = 0.150)")
axes[0].set_xlabel("null R2 (age model)")
axes[1].hist(perm["R2_null_label"], bins=30, color="#B0B0B0", edgecolor="white")
axes[1].axvline(0.505, color="#C44E52", lw=2)
axes[1].set_title("(b) label permutation\n(499 perms; p = 0.002)")
axes[1].set_xlabel("null R2 (age model)")
x = np.arange(len(cv))
axes[2].bar(x - 0.2, cv["loso_r2_age"], 0.4, label="LOSO R2 (age)", color="#4C72B0")
axes[2].bar(x + 0.2, cv["loso_auc_group"], 0.4, label="LOSO AUC (group)", color="#DD8452")
axes[2].set_xticks(x); axes[2].set_xticklabels(["audit PLS", "PCA", "PLS(age)", "ridge"], fontsize=8)
axes[2].axhline(0, color="grey", lw=0.5)
axes[2].set_title("(c) leave-one-subject-out validation")
axes[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUTF / "TA_fig2_audit.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- Fig 3: annotated top loadings ----------
top_idx = np.argsort(load_gene)[::-1][:20]
names, vals, stable = [], [], []
four = {"Nfu_g_1_010933", "Nfu_g_1_007097", "Nfu_g_1_001084", "Nfu_g_1_022498"}
for i in top_idx:
    gid = sel1000_list[i]
    nm = str(ann.loc[gid, "short_name"]) if gid in ann.index else ""
    if nm in ("nan", ""):
        nm = gid.replace("Nfu_g_1_", "g_")
    names.append(nm); vals.append(load_gene[i]); stable.append(gid in four)
fig, ax = plt.subplots(figsize=(6.5, 6))
colors = ["#C44E52" if s else "#4C72B0" for s in stable]
ax.barh(range(len(names)), vals, color=colors)
ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8)
ax.invert_yaxis()
ax.set_xlabel("|loading| on audit component 1 (sum over 10w / 20w / delta blocks)")
ax.set_title("top-20 loading genes, annotated\n(red = the four bootstrap-stable genes from the original report)")
fig.tight_layout()
fig.savefig(OUTF / "TA_fig3_annotated_loadings.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- Fig 4: circadian panel + complex I summary ----------
circ = pd.read_csv(OUTT / "TA_circadian_genes.csv")
ci = pd.read_csv(OUTT / "TA_complexI_members.csv")
enr = json.loads((OUTT / "TA_enrichment_summary.json").read_text())
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
cc = circ.sort_values("p_rank")
sig_color = ["#C44E52" if q < 0.05 else "#B0B0B0" for q in cc["q"]]
axes[0].barh(range(len(cc)), cc["rho"], color=sig_color)
axes[0].set_yticks(range(len(cc))); axes[0].set_yticklabels(cc["short_name"], fontsize=9)
axes[0].invert_yaxis(); axes[0].axvline(0, color="grey", lw=0.5)
axes[0].set_xlabel("Spearman rho (10->20wk delta vs age at death)")
axes[0].set_title("(a) core circadian-clock genes\n(red = FDR q < 0.05)")
ci_abs = np.abs(ci["delta_age_rho"].dropna())
all_abs = np.abs(st["rho_g"])
axes[1].hist(all_abs, bins=60, color="#B0B0B0", density=True, label="all genes")
axes[1].hist(ci_abs, bins=15, color="#C44E52", density=True, alpha=0.7, label="complex I (n=41)")
axes[1].set_xlabel("|Spearman rho| (delta vs age)"); axes[1].set_ylabel("density")
axes[1].set_title(f"(b) complex I vs background\nmedian |rho| 0.235 vs 0.135; enrichment {enr['complexI_enrichment']:.2f}x, p = {enr['complexI_binom_p']:.3f}")
axes[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUTF / "TA_fig4_circadian_complexI.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("figures saved")

# ---------- GSEA via Enrichr (human libraries on uppercase short names)
# Caveat recorded: killifish short names follow ortholog-style symbols; mapping to human
# libraries is the standard practice for non-model fish, with the known limits.
import gseapy

fdr = pd.read_csv(OUTT / "TA_gene_fdr.csv.gz")
if "gene_id" not in fdr.columns:          # older-format file: rows are in gene order
    fdr["gene_id"] = genes
fdr = fdr.join(ann[["short_name"]], on="gene_id")
fdr["sym"] = fdr["short_name"].astype(str).str.upper()
ok = fdr["sym"].str.match(r"^[A-Z0-9\-]+$") & ~fdr["sym"].isin(["NAN", ""])
fdr = fdr.sort_values("delta_age_rho", ascending=False)
up = fdr[ok].drop_duplicates("sym").head(200)["sym"].tolist()
dn = fdr[ok].drop_duplicates("sym").tail(200)["sym"].tolist()
print(f"GSEA queries: up {len(up)}, down {len(dn)} symbols")

LIBRARIES = ["GO_Biological_Process_2023", "KEGG_2021_Human"]
gsea_frames = []
for name, lst in [("delta_up", up), ("delta_down", dn)]:
    enr = gseapy.enrichr(gene_list=lst, gene_sets=LIBRARIES, organism="human",
                         outdir=None, cutoff=1.0)
    df = enr.results.sort_values("Adjusted P-value")
    df.insert(0, "direction", name)
    gsea_frames.append(df)
    print(f"Enrichr {name}: {len(df)} terms; top3:",
          "; ".join(df["Term"].head(3).astype(str)))
gsea = pd.concat(gsea_frames, ignore_index=True)
gsea.to_csv(OUTT / "TA_gsea_enrichr.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
for ax, name in zip(axes, ["delta_up", "delta_down"]):
    sub = gsea[gsea["direction"] == name].head(12).iloc[::-1]
    labels = [t[:52] for t in sub["Term"]]
    ax.barh(range(len(sub)), -np.log10(sub["Adjusted P-value"].clip(1e-20)),
            color="#55A868" if name == "delta_up" else "#C44E52")
    ax.set_yticks(range(len(sub))); ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_xlabel("-log10 FDR")
    ax.set_title(f"Enrichr GSEA — {name} (top 200 delta-age genes)")
fig.tight_layout()
fig.savefig(OUTF / "TA_fig5_gsea.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("GSEA figure saved")

print("Task A figures + GSEA done")
