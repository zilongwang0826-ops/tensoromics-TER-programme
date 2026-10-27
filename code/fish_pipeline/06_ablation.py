import os
import numpy as np, pandas as pd
from pathlib import Path
from sklearn.cross_decomposition import PLSRegression
DATA = Path("data")
d = np.load(DATA / "fish_data.npz", allow_pickle=True)
X = d["X"]; ages = d["ages"].astype(float); groups = d["groups"]; sel = d["sel1000"]
def r2_of(features):
    F = features.copy(); mu, sd = F.mean(0), F.std(0); sd[sd == 0] = 1; F = (F - mu) / sd
    pls = PLSRegression(n_components=2, scale=False).fit(F, groups.astype(float).reshape(-1, 1))
    T = pls.x_scores_.copy()
    for c in range(T.shape[1]):
        if np.corrcoef(T[:, c], ages)[0, 1] < 0: T[:, c] = -T[:, c]
    A = np.column_stack([np.ones(len(ages)), T]); coef, *_ = np.linalg.lstsq(A, ages, rcond=None)
    return 1 - ((ages - A @ coef) ** 2).sum() / ((ages - ages.mean()) ** 2).sum()
F10 = X[:, 0, sel]; F20 = X[:, 1, sel]
rows = [
  {"feature_block": "10w only", "r2_age": r2_of(F10)},
  {"feature_block": "20w only", "r2_age": r2_of(F20)},
  {"feature_block": "10w+20w", "r2_age": r2_of(np.hstack([F10, F20]))},
  {"feature_block": "10w+20w+delta (audit model)", "r2_age": r2_of(np.hstack([F10, F20, F20 - F10]))},
  {"feature_block": "delta only", "r2_age": r2_of(F20 - F10)},
]
pd.DataFrame(rows).round(4).to_csv("results/tables/TA_ablation.csv", index=False)
print(pd.DataFrame(rows).round(3).to_string())
