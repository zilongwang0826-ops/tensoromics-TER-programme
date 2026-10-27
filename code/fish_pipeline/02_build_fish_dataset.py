#!/usr/bin/env python3
"""Task A — Step 1: build the fish analysis dataset from GSE150318.

Inputs : GSE150318_counts.csv.gz, GSE150318_series_matrix.txt.gz, gene_annotation.tsv
Outputs: fish_data.npz with X[fish, 2, genes] (log2 CPM), ages, groups, gene ids,
         plus fish_metadata.tsv
"""
import os
import gzip
import re
import numpy as np
import pandas as pd

DATA = "os.path.join(os.environ.get("TER_ROOT", "."), "data", "GSE150318")"

# ---- 1. counts ----
print("loading counts...", flush=True)
counts = pd.read_csv(f"{DATA}/GSE150318_counts.csv.gz", index_col=0)
print("counts:", counts.shape)
sample_cols = list(counts.columns)
# sample names like FC1_10w / FC1_20w
meta = pd.DataFrame({
    "sample": sample_cols,
    "animal": [re.match(r"(FC\d+)_(\d+)w", s).group(1) for s in sample_cols],
    "week": [int(re.match(r"(FC\d+)_(\d+)w", s).group(2)) for s in sample_cols],
})
print("animals:", meta["animal"].nunique(), "| weeks:", sorted(meta["week"].unique()))

# ---- 2. age at death from series matrix ----
print("parsing series matrix...", flush=True)
titles, aad = [], {}
with gzip.open(f"{DATA}/GSE150318_series_matrix.txt.gz", "rt") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
        if line.startswith('!Sample_characteristics_ch1\t"age at death (days):'):
            vals = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            for t, v in zip(titles, vals):
                m = re.search(r"age at death \(days\): (\d+)", v)
                aad[t] = int(m.group(1)) if m else np.nan
print("age-at-death records:", len(aad))

# map (animal, week) -> title content-wise (order-independent)
title_key = {}
for t in titles:
    m = re.match(r"(FC\d+) fin clipped at (\d+) weeks", t)
    if m:
        title_key[(m.group(1), int(m.group(2)))] = t
print("parsed titles:", len(title_key))

# ---- 3. assemble per-fish matrix ----
meta["title"] = [title_key.get((a, w)) for a, w in zip(meta["animal"], meta["week"])]
meta["age_at_death"] = meta["title"].map(aad)
m10 = meta[meta["week"] == 10].set_index("animal")
m20 = meta[meta["week"] == 20].set_index("animal")
common = sorted(set(m10.index) & set(m20.index))
complete = [a for a in common if not np.isnan(m10.loc[a, "age_at_death"])
            and not np.isnan(m20.loc[a, "age_at_death"])]
same = [a for a in complete if m10.loc[a, "age_at_death"] == m20.loc[a, "age_at_death"]]
print(f"fish with both timepoints: {len(common)}; with age: {len(complete)}; consistent: {len(same)}")

fish = same
ages = np.array([m10.loc[a, "age_at_death"] for a in fish], dtype=float)
s10 = [m10.loc[a, "sample"] for a in fish]
s20 = [m20.loc[a, "sample"] for a in fish]
# counts is (genes x samples); build (fish, 2, genes)
X = np.stack([counts[s10].to_numpy(float).T, counts[s20].to_numpy(float).T], axis=1)
genes = counts.index.to_numpy()
print("X raw:", X.shape)

# ---- 4. CPM + log2 ----
cpm = X / X.sum(axis=2, keepdims=True) * 1e6
Xlog = np.log2(cpm + 1)
print("X:", Xlog.shape)

# ---- 5. gene selection: top-N most variable across the 2*fish observations ----
flat = Xlog.reshape(-1, Xlog.shape[2])
var = flat.var(axis=0)
order = np.argsort(var)[::-1]
sel = {n: order[:n] for n in [500, 1000, 2000]}
print("top-1000 var genes overlap with annotation:",
      sum(1 for g in genes[sel[1000]] if g in set(pd.read_csv(f"{DATA}/gene_annotation.tsv",
          sep="\t")["gene_id"])), "/1000")

# ---- 6. save ----
groups = (ages > np.median(ages)).astype(int)  # 1 = long-lived
print("median age:", np.median(ages), "| short/long:", np.bincount(groups))
np.savez_compressed(f"{DATA}/fish_data.npz",
                    X=Xlog, ages=ages, groups=groups, genes=genes,
                    sel500=sel[500], sel1000=sel[1000], sel2000=sel[2000],
                    fish=np.array(fish))
pd.DataFrame({"fish": fish, "age_at_death": ages, "group": groups}).to_csv(
    f"{DATA}/fish_metadata.tsv", sep="\t", index=False)
print("saved fish_data.npz")
print("done.")
