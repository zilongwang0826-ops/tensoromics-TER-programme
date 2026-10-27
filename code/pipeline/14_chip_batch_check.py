import gzip
import os
import re
from pathlib import Path

import pandas as pd

supp = {}
gsms = None
ROOT = Path(os.environ.get("TER_ROOT", os.getcwd()))
SERIES = Path(os.environ.get("GSE79396_SERIES_MATRIX",
                             str(ROOT / "data" / "GSE79396" / "GSE79396_series_matrix.txt.gz")))
if not SERIES.exists():
    raise SystemExit(
        f"[14] series matrix not found: {SERIES}\n"
        "Set GSE79396_SERIES_MATRIX or place the file under data/GSE79396/."
    )
print(f"[14] series matrix: {SERIES}")

with gzip.open(SERIES, "rt",
               encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"):
            gsms = [p.strip('"') for p in line.rstrip().split("\t")[1:]]
        elif line.startswith("!Sample_supplementary_file"):
            vals = [p.strip('"') for p in line.rstrip().split("\t")[1:]]
            for g, v in zip(gsms, vals):
                m = re.search(r"_(\d{10})_", v)
                supp[g] = m.group(1) if m else "NA"
        elif line.startswith("!series_matrix_table_begin"):
            break

samples = pd.read_csv(ROOT / "data" / "GSE79396" / "samples_all.csv")
samples["chip"] = samples["gsm"].map(supp)
print("chip x cohort x age_group (samples):")
print(pd.crosstab([samples["cohort"], samples["chip"]], samples["age_group"]).to_string())
print()
print("chip x cohort (samples):")
print(pd.crosstab(samples["chip"], samples["cohort"]).to_string())
subj = samples.drop_duplicates("subject_id")
print()
print("chip x cohort (subjects):")
print(pd.crosstab(subj["chip"], subj["cohort"]).to_string())
print()
print("chip x age_group (subjects):")
print(pd.crosstab(subj["chip"], subj["age_group"]).to_string())
