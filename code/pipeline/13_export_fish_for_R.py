# 13_export_fish_for_R.py
# Export the 113-fish tensor (1000 HVGs, 2 time points) as CSV for the R comparison script.
import os
import numpy as np
import pandas as pd
from pathlib import Path

SRC = Path(os.environ.get("TER_ROOT", ".")) / "data" / "GSE150318"
OUT = Path(os.environ.get("TER_ROOT", ".")) / "data" / "GSE150318"
OUT.mkdir(parents=True, exist_ok=True)

d = np.load(SRC / "fish_data.npz", allow_pickle=True)
print("keys:", list(d.keys()))
X = d["X"]                 # (113, 2, 26739) log2 CPM, order = [10w, 20w]
ages = d["ages"].astype(float)
groups = d["groups"]       # 1 = long-lived
genes = np.asarray(d["genes"])
sel = np.asarray(d["sel1000"])
print("X:", X.shape, "| sel:", len(sel), "| n_long:", int(groups.sum()))

Xsel = X[:, :, sel]        # (113, 2, 1000)
genes_sel = [str(g) for g in genes[sel]]

n, t, p = Xsel.shape
rows = []
for i in range(n):
    for k, tp in enumerate([10, 20]):
        rows.append([i + 1, tp] + Xsel[i, k, :].tolist())
cols = ["subj", "time"] + genes_sel
df = pd.DataFrame(rows, columns=cols)
df.to_csv(OUT / "fish_X1000_long.csv", index=False)
print("saved fish_X1000_long.csv:", df.shape)

meta = pd.DataFrame({
    "subj": np.arange(1, n + 1),
    "age_death": ages,
    "long_lived": groups.astype(int),
})
meta.to_csv(OUT / "fish_meta.csv", index=False)
print("saved fish_meta.csv:", meta.shape, "| median age:", np.median(ages))

# quick sanity: group sizes
print("short/long:", np.bincount(groups.astype(int)))
