#!/usr/bin/env python3
"""Build the Nfu_g_1_XXXXXX → annotation mapping from the original FLI gene build.

Source: Nfu_20150522.genes_20150922.gff3.gz (NFINgb, Leibniz-FLI) — the gene build that the
GSE150318 counts matrix is based on.
Output: gene_annotation.tsv with columns:
  gene_id, description, biotype, short_name, orthologs
"""
import os
import gzip
import csv
from urllib.parse import unquote

SRC = "Nfu_20150522.genes_20150922.gff3.gz"
OUT = "gene_annotation.tsv"

rows = []
with gzip.open(SRC, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        attrs = {}
        for kv in parts[8].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                attrs[k] = unquote(v)
        rows.append({
            "gene_id": attrs.get("ID", ""),
            "description": attrs.get("annotation", ""),
            "biotype": attrs.get("biotype", ""),
            "short_name": attrs.get("gene", ""),
            "orthologs": attrs.get("db_xref", ""),
        })

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene_id", "description", "biotype", "short_name", "orthologs"],
                       delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"wrote {OUT}: {len(rows)} genes")
biotypes = {}
for r in rows:
    biotypes[r["biotype"]] = biotypes.get(r["biotype"], 0) + 1
print("biotypes:", dict(sorted(biotypes.items(), key=lambda x: -x[1])[:8]))
