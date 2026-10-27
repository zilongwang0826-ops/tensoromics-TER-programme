# 12_inventory_Li2017.py
# Inventory the MMRN supplemental source data: subject x day coverage for each assay.
# Decides the multi-block (WP3) go/no-go and saves tidy copies of usable blocks.
import os
import re
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(os.environ.get("TER_ROOT", getwd()))
SRC = BASE / "data/Li2017_suppl/mmc2_extracted/SupplementalPackage-MMRN/sourcedata"
OUT = BASE / "data/Li2017_suppl/tidy"
OUT.mkdir(parents=True, exist_ok=True)
REPORT = BASE / "results/tables/wp0b_inventory.txt"

# transcriptome subjects (72)
import pandas as pd
samples = pd.read_csv(BASE / "data/GSE79396/samples_all.csv")
tx_subjects = sorted(samples["subject_id"].unique())
print(f"transcriptome subjects: {len(tx_subjects)}")

def norm_subject(s):
    # 'ZV1-01' -> 'ZV101'; 'ZV2-35' -> 'ZV235'
    m = re.match(r"^ZV(\d)-(\d+)$", str(s))
    return f"ZV{m.group(1)}{m.group(2)}" if m else str(s)

lines = []
def log(s):
    print(s)
    lines.append(s)

# ---------------- metabolomics ----------------
for mode in ["pos", "neg"]:
    path = SRC / f"ave_log2_QC_zerofil_vax05_{mode}_aplcms6_{'0904' if mode=='pos' else '0828'}.txt"
    df = pd.read_csv(path, sep="\t")
    meta_cols = ["mz", "time", "mz.min", "mz.max", "ave_CV"]
    smp_cols = [c for c in df.columns if c not in meta_cols]
    info = []
    for c in smp_cols:
        m = re.match(r"^(ZV\d-\d+)-d(\d+)$", c)
        if m:
            info.append((c, norm_subject(m.group(1)), int(m.group(2))))
    meta = pd.DataFrame(info, columns=["col", "subject", "day"])
    n_subj = meta["subject"].nunique()
    days = sorted(meta["day"].unique())
    log(f"[metab-{mode}] features={len(df)} samples={len(smp_cols)} subjects={n_subj} days={days}")
    log(f"  subjects also in transcriptome: {len(set(meta['subject']) & set(tx_subjects))}")
    for d in days:
        sub = meta[meta["day"] == d]
        log(f"  day {d}: {sub['subject'].nunique()} subjects, {len(sub)} samples")
    # subject x day table for the 72 transcriptome subjects
    key = meta[meta["subject"].isin(tx_subjects)]
    n_complete = (key.groupby("subject")["day"].nunique() == len(days)).sum()
    log(f"  transcriptome-subjects with all {len(days)} days present: {n_complete}/72")
    # save a tidy matrix for D0/1/3/7 subjects (features x samples) for later use
    keep_cols = key[key["day"].isin([0, 1, 3, 7])]["col"].tolist()
    mat = df[meta_cols + keep_cols].copy()
    mat.to_csv(OUT / f"metab_{mode}_d0_d7.csv.gz", index=False, compression="gzip")
    log(f"  saved {len(keep_cols)} D0-D7 columns -> metab_{mode}_d0_d7.csv.gz")

# ---------------- flow ----------------
flow = pd.read_csv(SRC / "reformatted4_Flow.txt", sep="\t")
subj_cols = [c for c in flow.columns if c != "Sample"]
flow_map = pd.DataFrame({"col": subj_cols, "subject": [norm_subject(c) for c in subj_cols]})
days_flow = flow["Sample"].str.extract(r"^D(\d+)_")[0].astype(float)
log(f"[flow] rows={len(flow)} (cell types x days) subject cols={len(subj_cols)} subjects={flow_map['subject'].nunique()}")
log(f"  days: {sorted(days_flow.dropna().unique())}")
log(f"  subjects in transcriptome: {len(set(flow_map['subject']) & set(tx_subjects))}")
flow.to_csv(OUT / "flow_raw.csv", index=False)

# ---------------- cytokines ----------------
cyto = pd.read_csv(SRC / "SL20160305_log2_formatted_cytokines3.txt", sep="\t")
cyto_cols = [c for c in cyto.columns if c != "Cytokines"]
cyto_days = cyto["Cytokines"].str.extract(r"^D(\d+)_")[0].astype(float)
log(f"[cytokines] rows={len(cyto)} subject cols={len(cyto_cols)} subjects={len(cyto_cols)}")
log(f"  days: {sorted(cyto_days.dropna().unique())}; cytokines: {cyto['Cytokines'].str.replace(r'^D\d+_','',regex=True).nunique()}")
log(f"  subjects in transcriptome: {len(set(cyto_cols) & set(tx_subjects))}")
cyto.to_csv(OUT / "cytokines_raw.csv", index=False)

# ---------------- BTMs ----------------
btm = pd.read_csv(SRC / "formatted_BTMs.txt", sep="\t")
btm_cols = [c for c in btm.columns if c != "BTM"]
btm_days = btm["BTM"].str.extract(r"^D(\d+)[/v]")[0].astype(float)
log(f"[BTM] rows={len(btm)} subjects={len(btm_cols)} days(prefix)={sorted(btm_days.dropna().unique())}")
log(f"  subjects in transcriptome: {len(set(btm_cols) & set(tx_subjects))}")
btm.to_csv(OUT / "btm_scores.csv", index=False)

# ---------------- responses ----------------
igg = pd.read_csv(SRC / "formated2_IgG_NLSrenalyzed.txt", sep="\t")
igg_cols = [c for c in igg.columns if c != "subj"]
log(f"[IgG] rows={igg['subj'].tolist()[:6]}... subjects={len(igg_cols)}")
log(f"  subjects in transcriptome: {len(set(igg_cols) & set(tx_subjects))}")
igg.to_csv(OUT / "igg_raw.csv", index=False)

tfh = pd.read_csv(SRC / "formatted2_TFH.txt", sep="\t")
tfh_cols = [c for c in tfh.columns if c != "Subj"]
log(f"[Tfh] rows={tfh['Subj'].tolist()[:6]}... subjects={len(tfh_cols)}")
log(f"  subjects in transcriptome: {len(set(tfh_cols) & set(tx_subjects))}")
tfh.to_csv(OUT / "tfh_raw.csv", index=False)

ifng = pd.read_csv(SRC / "Tcell_IFNg3.txt", sep="\t")
ifng_cols = [c for c in ifng.columns if c != "Subj"]
log(f"[Tcell IFNg] rows={ifng['Subj'].tolist()[:6]}... subjects={len(ifng_cols)}")
ifng.to_csv(OUT / "tcell_ifng_raw.csv", index=False)

cd4 = pd.read_csv(SRC / "formatted2_CD4.txt", sep="\t")
cd4_cols = [c for c in cd4.columns if c != "Subj"]
log(f"[CD4] rows={cd4['Subj'].tolist()[:6]}... subjects={len(cd4_cols)}")
cd4.to_csv(OUT / "cd4_raw.csv", index=False)

# demographics validation
demo = pd.read_csv(SRC / "VZV_age_gender2.txt", sep="\t")
demo["subject"] = demo["ZV"].map(norm_subject)
log(f"[demographics] n={len(demo)} subjects={demo['subject'].nunique()}")
m = samples.drop_duplicates("subject_id")[["subject_id", "age", "sex", "age_group"]].merge(
    demo, left_on="subject_id", right_on="subject", how="inner")
log(f"  demographics joined to transcriptome: {len(m)}"
    f" | age match: {(m['age'] == m['AGE']).all()}"
    f" | sex match after F/M: {(m['sex'].str[0].str.upper() == m['Gender']).all()}")

log("")
log("=== GO/NO-GO summary for multi-block (WP3) ===")
log("- metabolome pos+neg: subject-level, days 0/1/3/7 present for 72-subject set -> GO (primary block 2)")
log("- flow, cytokines: check coverage above; usable as smaller/additional blocks if days>=2 and n>=40")
log("- BTMs are transcriptome-derived: use for gene-set testing (WP2E), not as independent block")

REPORT.write_text("\n".join(lines), encoding="utf-8")
print("\nreport saved:", REPORT)
