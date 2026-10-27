#!/usr/bin/env python3
"""Task A — Step 2: fish lifespan analysis with the P0 audit upgrades.

Reproduces the audit-style trajectory PLS (10w expression | 20w expression | 20w-10w delta,
supervised on the long-lived label) and adds the four P0 upgrades from the TER optimisation list:
  1. structure-preserving permutation (within-fish time reshuffle) [+ label permutation]
  2. leave-one-subject-out CV + matrix baselines + bootstrap CI
  3. rank / feature-selection sensitivity + FDR on gene-level tests
  4. mitochondrial complex I hypothesis test (annotation closure, Fisher enrichment)

Outputs: ../results/{tables,figures}/TA_*.csv|png
"""
import os
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.metrics import roc_auc_score

DATA = Path("os.path.join(os.environ.get("TER_ROOT", "."), "data", "GSE150318")")
OUTF = Path("os.path.join(os.environ.get("TER_ROOT", "."), "results", "figures")")
OUTT = Path("os.path.join(os.environ.get("TER_ROOT", "."), "results", "tables")")
OUTF.mkdir(parents=True, exist_ok=True)
OUTT.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(0)
SEED = 0
N_PERM = 499
N_BOOT = 1000

# ---------------------------------------------------------------- load
d = np.load(DATA / "fish_data.npz", allow_pickle=True)
X = d["X"]                    # (113, 2, 26739) log2 CPM
ages = d["ages"].astype(float)
groups = d["groups"]          # 1 = long-lived
genes = d["genes"]
sel1000 = d["sel1000"]
print("fish:", X.shape[0], "| genes:", X.shape[2], "| short/long:", np.bincount(groups))


def audit_features(X, sel):
    """[10w | 20w | delta] feature matrix, standardized; returns (F, meta)"""
    F10 = X[:, 0, sel]
    F20 = X[:, 1, sel]
    F = np.hstack([F10, F20, F20 - F10])
    mu, sd = F.mean(0), F.std(0)
    sd[sd == 0] = 1.0
    return (F - mu) / sd


def fit_audit_pls(F, y, ages_ref, n_comp=2):
    """group-supervised PLS; components sign-oriented so positive = longer lifespan.
    Orientation uses ages_ref only (train folds in CV, all fish otherwise)."""
    pls = PLSRegression(n_components=n_comp, scale=False).fit(F, y.reshape(-1, 1))
    T = pls.x_scores_.copy()
    signs = np.ones(n_comp)
    for c in range(T.shape[1]):
        if np.corrcoef(T[:, c], ages_ref)[0, 1] < 0:
            T[:, c] *= -1
            pls.x_loadings_[:, c] *= -1
            pls.x_weights_[:, c] *= -1
            signs[c] = -1
    return pls, T, signs


def age_r2(T, ages):
    """OLS R2 of age ~ T (with intercept); T: (n, k) already oriented"""
    A = np.column_stack([np.ones(len(ages)), T])
    coef, res, *_ = np.linalg.lstsq(A, ages, rcond=None)
    pred = A @ coef
    ss_res = ((ages - pred) ** 2).sum()
    ss_tot = ((ages - ages.mean()) ** 2).sum()
    return 1 - ss_res / ss_tot, pred


# ---------------------------------------------------------------- main model
F = audit_features(X, sel1000)
pls, T, _ = fit_audit_pls(F, groups.astype(float), ages, n_comp=2)
R2, age_pred = age_r2(T, ages)
r1 = np.corrcoef(T[:, 0], ages)[0, 1]
r2 = np.corrcoef(T[:, 1], ages)[0, 1]
print(f"[main] audit PLS(2) age R2 = {R2:.3f} | r(comp1,age) = {r1:.3f} | r(comp2,age) = {r2:.3f}")

# supplementary: how much does the temporal DELTA add on top of cross-sectional levels?
F10 = X[:, 0, sel1000]; F20 = X[:, 1, sel1000]
Flev = np.hstack([F10, F20]); mu, sd = Flev.mean(0), Flev.std(0); sd[sd == 0] = 1.0
Flev = (Flev - mu) / sd
_, Tlev, _ = fit_audit_pls(Flev, groups.astype(float), ages, n_comp=2)
R2_lev, _ = age_r2(Tlev, ages)
print(f"[levels-only] audit PLS(2) without delta features: R2 = {R2_lev:.3f} (vs {R2:.3f} with delta)")

# ---------------------------------------------------------------- P0.1 permutations
p_time, R2_null_time, r1_null_time = [], [], []
for _ in range(N_PERM):
    Xp = X.copy()
    swap = rng.random(X.shape[0]) < 0.5
    Xp[swap, 0, :], Xp[swap, 1, :] = Xp[swap, 1, :], Xp[swap, 0, :].copy()
    Fp = audit_features(Xp, sel1000)
    _, Tp, _ = fit_audit_pls(Fp, groups.astype(float), ages, n_comp=2)
    r2p, _ = age_r2(Tp, ages)
    R2_null_time.append(r2p)
    r1_null_time.append(np.corrcoef(Tp[:, 0], ages)[0, 1])
R2_null_time = np.array(R2_null_time); r1_null_time = np.array(r1_null_time)
p_time = (np.sum(R2_null_time >= R2) + 1) / (N_PERM + 1)

R2_null_label, r1_null_label = [], []
for _ in range(N_PERM):
    yp = rng.permutation(groups).astype(float)
    _, Tp, _ = fit_audit_pls(F, yp, ages, n_comp=2)
    r2p, _ = age_r2(Tp, ages)
    R2_null_label.append(r2p)
    r1_null_label.append(np.corrcoef(Tp[:, 0], ages)[0, 1])
R2_null_label = np.array(R2_null_label); r1_null_label = np.array(r1_null_label)
p_label = (np.sum(R2_null_label >= R2) + 1) / (N_PERM + 1)
print(f"[perm] time-reshuffle p(R2) = {p_time:.3f} | label-permutation p(R2) = {p_label:.3f}")

# ---------------------------------------------------------------- P0.2 LOSO CV + baselines
def loso_predict(Fm, ages, groups, model="audit", n_comp=2):
    """out-of-sample predictions over fish; every orientation decision uses train folds only"""
    n = len(ages)
    pred_age = np.zeros(n); pred_score = np.zeros(n)
    for i in range(n):
        tr = np.arange(n) != i
        Ftr, Fte = Fm[tr], Fm[[i]]
        if model == "audit":
            pls_i, Ttr, signs = fit_audit_pls(Ftr, groups[tr].astype(float), ages[tr], n_comp=n_comp)
            Tte = pls_i.transform(Fte) * signs
            A = np.column_stack([np.ones(tr.sum()), Ttr])
            coef, *_ = np.linalg.lstsq(A, ages[tr], rcond=None)
            pred_age[i] = np.r_[1, Tte[0]] @ coef
            pred_score[i] = (Tte[0, 0] - Ttr[:, 0].mean()) / (Ttr[:, 0].std() + 1e-12)
        elif model == "pca":
            pca_i = PCA(n_components=n_comp, random_state=SEED).fit(Ftr)
            Ttr = pca_i.transform(Ftr); Tte = pca_i.transform(Fte)
            for c in range(n_comp):
                if np.corrcoef(Ttr[:, c], ages[tr])[0, 1] < 0:
                    Ttr[:, c] *= -1; Tte[:, c] *= -1
            A = np.column_stack([np.ones(tr.sum()), Ttr])
            coef, *_ = np.linalg.lstsq(A, ages[tr], rcond=None)
            pred_age[i] = np.r_[1, Tte[0]] @ coef
            pred_score[i] = (Tte[0, 0] - Ttr[:, 0].mean()) / (Ttr[:, 0].std() + 1e-12)
        elif model == "pls_age":
            pls_i = PLSRegression(n_components=n_comp, scale=False).fit(Ftr, ages[tr])
            pred_age[i] = float(np.ravel(pls_i.predict(Fte))[0])
            pred_score[i] = pred_age[i]  # higher predicted age = longer-lived
        elif model == "ridge":
            r_i = RidgeCV(alphas=np.logspace(-2, 4, 13)).fit(Ftr, ages[tr])
            pred_age[i] = float(np.ravel(r_i.predict(Fte))[0])
            pred_score[i] = pred_age[i]
    from sklearn.metrics import r2_score
    r2_loso = r2_score(ages, pred_age)  # 1 - SSres/SStot vs mean baseline
    auc = roc_auc_score(groups, pred_score)
    return r2_loso, auc, pred_age

cv_table = []
for name in ["audit", "pca", "pls_age", "ridge"]:
    r2_cv, auc_cv, _ = loso_predict(F, ages, groups, model=name)
    cv_table.append({"model": name, "loso_r2_age": r2_cv, "loso_auc_group": auc_cv})
    print(f"[LOSO] {name:8s}: R2(age) = {r2_cv:+.3f} | AUC(group) = {auc_cv:.3f}")

# bootstrap CI (over fish) for in-sample R2 & r1
boot_r2, boot_r1 = [], []
n = len(ages)
for _ in range(N_BOOT):
    idx = rng.integers(0, n, n)
    Fb = F[idx]
    mu, sd = Fb.mean(0), Fb.std(0); sd[sd == 0] = 1.0
    _, Tb, _ = fit_audit_pls((Fb - mu) / sd, groups[idx].astype(float), ages[idx], n_comp=2)
    r2b, _ = age_r2(Tb, ages[idx])
    boot_r2.append(r2b); boot_r1.append(np.corrcoef(Tb[:, 0], ages[idx])[0, 1])
boot_r2 = np.array(boot_r2); boot_r1 = np.array(boot_r1)
ci_r2 = np.percentile(boot_r2, [2.5, 97.5]); ci_r1 = np.percentile(boot_r1, [2.5, 97.5])
print(f"[boot] R2 95% CI = [{ci_r2[0]:.3f}, {ci_r2[1]:.3f}] | r1 CI = [{ci_r1[0]:.3f}, {ci_r1[1]:.3f}]")

# ---------------------------------------------------------------- P0.3 sensitivity + FDR
sens = []
for n_genes in [500, 1000, 2000]:
    sel = d[f"sel{n_genes}"]
    Fg = audit_features(X, sel)
    for nc in [1, 2, 3, 4]:
        _, Tg, _ = fit_audit_pls(Fg, groups.astype(float), ages, n_comp=nc)
        r2g, _ = age_r2(Tg, ages)
        sens.append({"n_genes": n_genes, "n_comp": nc, "in_sample_r2": r2g})
sens = pd.DataFrame(sens)
print("[sens] in-sample R2 grid:")
print(sens.pivot(index="n_comp", columns="n_genes", values="in_sample_r2").round(3).to_string())

# FDR on gene-level delta-vs-age correlations (all genes); NaN p (constant genes) -> 1
import warnings
delta = X[:, 1, :] - X[:, 0, :]
rho_g = np.zeros(X.shape[2]); p_g = np.ones(X.shape[2])
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for g in range(X.shape[2]):
        rho_g[g], p_g[g] = spearmanr(delta[:, g], ages)
p_g = np.where(np.isnan(p_g), 1.0, p_g)
rho_g = np.where(np.isnan(rho_g), 0.0, rho_g)
m = len(p_g)
order = np.argsort(p_g)
p_sorted = p_g[order]
bh = p_sorted * m / np.arange(1, m + 1)
bh = np.minimum.accumulate(bh[::-1])[::-1]   # BH step-up: q_i = min_{j>=i} (m/j) p_(j)
q_g = np.empty(m)
q_g[order] = np.minimum(bh, 1.0)
n_sig = int((q_g < 0.05).sum())
print(f"[FDR] genes with q<0.05 (delta vs age): {n_sig}/{m}")
if n_sig:
    top_hits = pd.DataFrame({"gene_id": genes[order[:20]], "rho": rho_g[order[:20]],
                             "p": p_g[order[:20]], "q": q_g[order[:20]]})
    print(top_hits.to_string())

# ---------------------------------------------------------------- P0.4 complex I
ann = pd.read_csv(DATA / "gene_annotation.tsv", sep="\t").set_index("gene_id")
ci_ids = [g for g in ann.index
          if str(ann.loc[g, "short_name"]).upper().startswith("NDUF")
          and not str(ann.loc[g, "short_name"]).upper().startswith("NDUFAF")]
ci_ids += [g for g in ann.index if str(ann.loc[g, "description"]).lower().startswith("nadh dehydrogenase")]
ci_ids = list(dict.fromkeys(ci_ids))
gene_pos = {g: k for k, g in enumerate(genes)}
ci_present = [g for g in ci_ids if g in gene_pos]
print(f"[complexI] annotated members: {len(ci_ids)} | in expression matrix: {len(ci_present)}")

# loading of each complex I gene in the main model (absolute loading on comp1, folded across the 3 feature blocks)
load = np.abs(pls.x_loadings_[:, 0])
load_gene = load[:1000] + load[1000:2000] + load[2000:3000]  # by selected gene index
sel1000_list = [str(genes[i]) for i in sel1000]
load_map = dict(zip(sel1000_list, load_gene))
ci_in_hvg = [g for g in ci_present if g in load_map]
print(f"[complexI] members among the 1000 HVGs: {len(ci_in_hvg)}/{len(ci_present)}")

ci_rows = []
for g in ci_present:
    k = gene_pos[g]
    ci_rows.append({
        "gene_id": g,
        "short_name": ann.loc[g, "short_name"],
        "description": ann.loc[g, "description"],
        "in_hvg1000": g in load_map,
        "delta_age_rho": rho_g[k], "delta_age_q": q_g[k],
        "loading_comp1": load_map.get(g, np.nan),
    })
ci_df = pd.DataFrame(ci_rows).sort_values("loading_comp1", ascending=False)
ci_df.round(4).to_csv(OUTT / "TA_complexI_members.csv", index=False)
print(ci_df.head(12).to_string())

# Fisher enrichment: complex I among top-100 |loading| genes, within the 1000-HVG universe
top100 = set(np.array(sel1000_list)[np.argsort(load_gene)[::-1][:100]])
ci_in_top = len([g for g in ci_in_hvg if g in top100])
from scipy.stats import fisher_exact
table = [[ci_in_top, len(ci_in_hvg) - ci_in_top],
         [len(top100) - ci_in_top, 1000 - len(top100) - (len(ci_in_hvg) - ci_in_top)]]
odds, p_fisher = fisher_exact(table, alternative="greater")
print(f"[complexI] top-100 loadings contain {ci_in_top}/{len(ci_in_hvg)} complex I genes (of the HVG1000); Fisher OR={odds:.2f}, p={p_fisher:.4f}")
ci_sig = int((ci_df["delta_age_q"] < 0.05).sum())
print(f"[complexI] complex-I genes with delta-age q<0.05: {ci_sig}/{len(ci_present)}")
if ci_sig:
    print(ci_df[ci_df["delta_age_q"] < 0.05].to_string())

# top-30 loading genes of the main model, annotated (the annotation-closure deliverable)
top_idx = np.argsort(load_gene)[::-1][:30]
top_rows = []
for r, i in enumerate(top_idx):
    gid = sel1000_list[i]
    top_rows.append({
        "rank": r + 1, "gene_id": gid, "loading_comp1": float(load_gene[i]),
        "short_name": ann.loc[gid, "short_name"] if gid in ann.index else "",
        "description": ann.loc[gid, "description"] if gid in ann.index else "",
        "delta_age_rho": float(rho_g[gene_pos[gid]]), "delta_age_q": float(q_g[gene_pos[gid]]),
    })
pd.DataFrame(top_rows).round(4).to_csv(OUTT / "TA_top_loadings_annotated.csv", index=False)
print("\n[top loadings annotated]")
print(pd.DataFrame(top_rows).head(15)[["rank", "gene_id", "short_name", "description", "delta_age_rho", "delta_age_q"]].to_string())

# ---------------------------------------------------------------- save everything
pd.DataFrame(cv_table).round(4).to_csv(OUTT / "TA_loso_cv.csv", index=False)
sens.round(4).to_csv(OUTT / "TA_sensitivity.csv", index=False)
pd.DataFrame({
    "gene_id": [str(g) for g in genes], "delta_age_rho": rho_g,
    "delta_age_p": p_g, "delta_age_q": q_g,
}).to_csv(OUTT / "TA_gene_fdr.csv.gz", index=False)
np.savez_compressed(OUTT / "TA_permutations.npz",
                    R2_null_time=R2_null_time, r1_null_time=r1_null_time,
                    R2_null_label=R2_null_label, r1_null_label=r1_null_label,
                    boot_r2=boot_r2, boot_r1=boot_r1)
summary = {
    "n_fish": int(len(ages)), "median_age": float(np.median(ages)),
    "n_short": int((groups == 0).sum()), "n_long": int((groups == 1).sum()),
    "audit_r2": float(R2), "audit_r1": float(r1), "audit_r2_comp2": float(r2),
    "p_time_perm": float(p_time), "p_label_perm": float(p_label),
    "boot_r2_ci": [float(ci_r2[0]), float(ci_r2[1])],
    "boot_r1_ci": [float(ci_r1[0]), float(ci_r1[1])],
    "loso": cv_table, "fdr_n_sig_q05": n_sig,
    "complexI": {"n_members": len(ci_present), "in_top100": ci_in_top,
                 "fisher_or": float(odds), "fisher_p": float(p_fisher),
                 "n_q05": ci_sig},
}
(OUTT / "TA_summary.json").write_text(json.dumps(summary, indent=2))
np.savez_compressed(DATA / "TA_model_state.npz",
                    T=T, ages=ages, groups=groups, rho_g=rho_g, q_g=q_g,
                    load_gene=load_gene, sel1000=sel1000)
print("\nTask A analysis done.")
