#!/usr/bin/env Rscript
# 25_fdr_corrections.R — Benjamini-Hochberg corrections within each analysis's natural
# family (P2B-audited). Reads the committed result tables and writes fdr_* candidate
# files + a summary. Raw p values are never altered.
suppressPackageStartupMessages(library(data.table))
PROJ <- Sys.getenv("TER_ROOT", getwd())
TAB  <- file.path(PROJ, "results/tables")
rd <- function(f) fread(file.path(TAB, f))
bh <- function(p) p.adjust(p, method = "BH")

log_rows <- list()
record <- function(family, variant, question, classification, tests, raw, note = "") {
  stopifnot(length(tests) == length(raw))
  padj <- bh(raw)
  out <- data.frame(test = tests, raw_p = raw, BH = padj, stringsAsFactors = FALSE)
  fwrite(out, file.path(TAB, paste0("fdr_", family, "__", variant, ".csv")))
  log_rows[[length(log_rows) + 1]] <<- data.frame(
    family = family, variant = variant, question = question, classification = classification,
    n_tests = length(raw), n_raw_sig = sum(raw < 0.05), n_BH_sig = sum(padj < 0.05),
    min_raw = min(raw), min_BH = min(padj), correction = "BH within variant",
    note = note, stringsAsFactors = FALSE)
  invisible(NULL)
}

## FAM21_POOL
bm <- rd("vacc_batch_metrics.csv")
for (i in seq_len(nrow(bm))) {
  record("FAM21_POOL", paste0("param_", bm$version[i]),
         "Cohort/age structure in tPCA/tPLS-DA per correction level", "confirmatory",
         c("cohort_tPCA", "age_tPCA", "age_tPLSDA"),
         c(bm$cohort_p_tPCA[i], bm$age_p_tPCA[i], bm$age_p_tPLSDA[i]),
         "parametric Pillai p")
  record("FAM21_POOL", paste0("perm_", bm$version[i]),
         "Same structure under subject-level permutation (N=299)", "confirmatory",
         c("perm_cohort_tPCA", "perm_age_tPCA", "perm_age_tPLSDA"),
         c(bm$p_perm_cohort_tPCA[i], bm$p_perm_age_tPCA[i], bm$p_perm_age_tPLSDA[i]), "")
}
lab <- character(); raw <- numeric()
for (i in seq_len(nrow(bm))) {
  lab <- c(lab, paste0(bm$version[i], "_", c("cohort_tPCA", "age_tPCA", "age_tPLSDA")))
  raw <- c(raw, bm$cohort_p_tPCA[i], bm$age_p_tPCA[i], bm$age_p_tPLSDA[i])
}
record("FAM21_POOL", "param_all9", "Strict single family across 9 parametric p values",
       "exploratory", lab, raw, "reference variant; mixes correction levels")

## FAM21B_COMP (full family if 21b2 output present)
rec_path <- file.path(TAB, "vacc_heterogeneity_comps_full.csv")
if (file.exists(rec_path)) {
  hc <- fread(rec_path)
  for (i in seq_len(nrow(hc))) {
    record("FAM21B_COMP", paste0(hc$cohort[i], "_", hc$version[i]),
           "Per-component age t-tests on within-cohort tPCA (5 components)", "confirmatory",
           paste0("comp", 1:5, "_p"),
           c(hc$comp1_p[i], hc$comp2_p[i], hc$comp3_p[i], hc$comp4_p[i], hc$comp5_p[i]),
           "best_comp_p (a min) excluded from BH")
  }
  record("FAM21B_COMP", "p5_perm", "Within-cohort 5-component permutation p (199 perms)",
         "confirmatory", paste0(hc$cohort, "_", hc$version), hc$p5_perm, "floor 0.005")
}

## FAM22_TPCA5
tp5 <- readLines(file.path(TAB, "vacc_chip_tpca5.txt"))
pparam <- as.numeric(sub(".*p_param[^0-9]*([0-9.]+).*", "\\1", grep("p_param", tp5, value = TRUE)))
pperm  <- as.numeric(sub(".*p_perm[^0-9]*([0-9.]+).*", "\\1", grep("p_perm", tp5, value = TRUE)))
record("FAM22_TPCA5", "param_perm", "Chip-corrected 5-component tPCA age association",
       "descriptive", c("p_param", "p_perm"), c(pparam, pperm), "n=2; permutation floor 0.002")

## FAM21C_CHIP
cc <- rd("vacc_chip_correction.csv")
for (i in seq_len(nrow(cc))) {
  record("FAM21C_CHIP", paste0("age_", cc$version[i]),
         "Age association before/after chip-level ComBat", "confirmatory",
         c("age_p_tPCA", "age_p_tPCA_perm", "age_p_tPLSDA_perm"),
         c(cc$age_p_tPCA[i], cc$age_p_tPCA_perm[i], cc$age_p_tPLSDA_perm[i]), "")
}
c1652 <- rd("vacc_chip1652_age.csv")
combat_row <- which(cc$version == "ComBat by chip")
record("FAM21C_CHIP", "age_after_chip_with1652",
       "Age null after chip correction incl. the mixed plate (n=24)", "confirmatory",
       c("age_ComBat_tPCA", "age_ComBat_tPCA_perm", "age_ComBat_tPLSDA_perm", "within_plate1652_age"),
       c(cc$age_p_tPCA[combat_row], cc$age_p_tPCA_perm[combat_row], cc$age_p_tPLSDA_perm[combat_row],
         c1652$p[1]), "primary family")
record("FAM21C_CHIP", "chip_effect", "Chip effect on tPCA before/after correction",
       "descriptive", paste0(cc$version, "_chip_p_tPCA"), cc$chip_p_tPCA, "separate question")

## FAM21D_COV
ca <- rd("vacc_chip_adjustment.csv")
record("FAM21D_COV", "all6", "Age term across six corrected spaces/models", "exploratory",
       paste0(ca$space, "_", ca$model), ca$age_stat, "chip_stat all NA in stored table")

## FAM21_WITHIN
wc <- rd("vacc_within_cohort.csv")
for (i in seq_len(nrow(wc))) {
  record("FAM21_WITHIN", paste0("cohort", i), "Within-cohort tPCA age association",
         "exploratory", c("tPCA_p_parametric", "tPCA_p_perm"),
         c(wc$tPCA_p_parametric[i], wc$tPCA_p_perm[i]), "perm floor 0.005")
}

## FAM20_PERM
pi_ <- rd("vacc_perm_inflation.csv")
record("FAM20_PERM", "ref_all6", "Permutation-scheme comparison panel", "descriptive",
       paste(pi_$analysis, pi_$scheme), pi_$p_perm, "BH reference only")

## FAM22_V1G
vg <- rd("vacc_v1genes_chip_adjust.csv")
record("FAM22_V1G", "age_chip_sex", "v1 genes adjusted for plate+sex", "confirmatory",
       vg$symbol, vg$p_age_chip_sex, "the reported adjusted null")
record("FAM22_V1G", "age_only", "v1 genes unadjusted", "exploratory", vg$symbol, vg$p_age_only)
lab <- character(); raw <- numeric()
for (i in seq_len(nrow(vg))) {
  lab <- c(lab, paste0(vg$symbol[i], "_age_only"), paste0(vg$symbol[i], "_age_chip_sex"))
  raw <- c(raw, vg$p_age_only[i], vg$p_age_chip_sex[i])
}
record("FAM22_V1G", "joint18", "Joint reference family", "exploratory", lab, raw)

## FAM22_TOP8
t8 <- rd("vacc_top8_mixed.csv")
record("FAM22_TOP8", "p_age", "Age term for the 8 top-loading genes", "exploratory",
       t8$symbol, t8$p_age, "genes selected by full-data loadings (selection bias)")
record("FAM22_TOP8", "p_interaction", "Age x day interaction for the same genes",
       "exploratory", t8$symbol, t8$p_interaction)
lab <- character(); raw <- numeric()
for (i in seq_len(nrow(t8))) {
  lab <- c(lab, paste0(t8$symbol[i], "_age"), paste0(t8$symbol[i], "_int"))
  raw <- c(raw, t8$p_age[i], t8$p_interaction[i])
}
record("FAM22_TOP8", "joint16", "Joint reference family", "exploratory", lab, raw)

## FAM22B_ANOVA / FAM22B_COMP
an <- rd("vacc_axis_dynamics_anova.csv")
record("FAM22B_ANOVA", "all5", "Day-resolved axis score ANOVA terms", "descriptive",
       an$effect, an$p, "model terms are not independent hypotheses")
comp_txt <- readLines(file.path(TAB, "vacc_tpls_comp_age_p.txt"))
vals <- as.numeric(regmatches(comp_txt, gregexpr("[0-9.]+e[-+][0-9]+", comp_txt))[[1]])
record("FAM22B_COMP", "comp1_comp2", "tPLS-DA component age LM (uncorrected)", "descriptive",
       c("comp1_p", "comp2_p"), vals[1:2])

## FAM30_BLOCK / CANON / ALIGN
bt <- rd("vacc_multiblock_blocktplsda.csv")
perm_txt <- readLines(file.path(TAB, "vacc_multiblock_blocktplsda_perm.txt"))
pm <- suppressWarnings(as.numeric(regmatches(perm_txt, regexpr("[0-9]*\\.?[0-9]+", perm_txt))))
record("FAM30_BLOCK", "blocks_and_perm", "Per-block Pillai + expr-block permutation",
       "exploratory", c(paste0(bt$block, "_pillai"), "expr_perm"), c(bt$p, pm[1]),
       "perm floor 0.010")
can <- rd("vacc_multiblock_canonical.csv")
record("FAM30_CANON", "comp1", "Canonical comp1 pairing permutation (n=1)", "exploratory",
       "pairing_perm_p", can$pairing_perm_p[!is.na(can$pairing_perm_p)], "floor 0.01")
al <- rd("vacc_multiblock_align.csv")
record("FAM30_ALIGN", "all9", "Alignment table (9 block x variable tests)", "exploratory",
       paste(al$block, al$variable), al$p, "Pillai p")
for (v in unique(al$variable)) {
  sub <- al[variable == v]
  record("FAM30_ALIGN", paste0("byvar_", v), paste0("Alignment for variable ", v),
         "exploratory", paste(sub$block, sub$variable), sub$p)
}

## FAM30_MANMETAB (metabolome audit; persisted by script 30 in the release pipeline)
ma_path <- file.path(TAB, "vacc_metab_audit.csv")
if (file.exists(ma_path)) {
  ma <- fread(ma_path)
  record("FAM30_MANMETAB", "pillai_all",
         "Metabolome tensor age/cohort/chip Pillai p (pre-integration audit)", "descriptive",
         paste0(ma$variable, "_pillai"), ma$p, "persisted from script 30; BH reference")
}

summary <- rbindlist(log_rows)
fwrite(summary, file.path(TAB, "fdr_corrections_summary.csv"))
cat("FDR corrections written;", nrow(summary), "variants across",
    length(unique(summary$family)), "families\n")
cat("DONE\n")

