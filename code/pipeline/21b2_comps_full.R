#!/usr/bin/env Rscript
# 21b2_comps_full.R — store the full per-component age t-tests for script 21b
# (comp1..comp5 + p5_perm), computed with the SAME code path and stored objects,
# validated against the committed summary columns before being trusted.
suppressPackageStartupMessages({ library(data.table); library(tensorOmics) })
PROJ <- Sys.getenv("TER_ROOT", getwd())
source(file.path(PROJ, "scripts/lib_tensor_helpers.R"))
prep <- readRDS(file.path(PROJ, "data/GSE79396/preprocessed_A.rds"))
corr <- readRDS(file.path(PROJ, "data/GSE79396/corrected_versions.rds"))
samples <- fread(file.path(PROJ, "data/GSE79396/samples_all.csv"))
tab_dir <- file.path(PROJ, "results/tables")
stored <- fread(file.path(tab_dir, "vacc_heterogeneity_comps.csv"))

subs <- prep$subs
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
cohort_subj <- factor(samples$cohort[match(subs, samples$subject_id)])
set.seed(11)

build_tensor_scaled <- function(Xrow) {
  Z <- scale(Xrow)
  Xt <- array(NA_real_, dim = c(length(subs), ncol(Z), 4),
              dimnames = list(subs, colnames(Z), paste0("D", 1:4)))
  for (s in subs) for (k in 1:4) {
    d <- c(0, 1, 3, 7)[k]
    Xt[s, , paste0("D", k)] <- Z[samples$subject_id == s & samples$day == d, ]
  }
  Xt
}
perm_p_unsup <- function(S, g, n = 199) {
  obs <- pillai(S, g)["trace"]
  (1 + sum(sapply(seq_len(n), function(b) pillai(S, sample(g))["trace"]) >= obs)) / (1 + n)
}
Xt_before <- prep$X_tensor
Xt_combat <- build_tensor_scaled(corr$combat_rows)

rows <- list()
for (ver in c("before", "combat")) {
  Xt <- if (ver == "before") Xt_before else Xt_combat
  for (co in levels(cohort_subj)) {
    idx <- which(cohort_subj == co)
    g <- droplevels(grp_subj[idx])
    tp <- tpca(Xt[idx, , , drop = FALSE], ncomp = 5)$variates
    pv <- sapply(1:5, function(k) t.test(tp[, k] ~ g)$p.value)
    rows[[length(rows) + 1]] <- data.frame(
      version = ver, cohort = co, n = length(idx),
      comp1_p = pv[1], comp2_p = pv[2], comp3_p = pv[3], comp4_p = pv[4], comp5_p = pv[5],
      best_comp_p = min(pv), best_comp = which.min(pv),
      p5_perm = perm_p_unsup(tp[, 1:5], g))
  }
}
rec <- rbindlist(rows)
m <- merge(rec, stored, by = c("version", "cohort"), suffixes = c("_new", "_stored"))
checks <- c(
  comp2 = max(abs(m$comp2_p_new - m$comp2_p_stored)),
  comp5 = max(abs(m$comp5_p_new - m$comp5_p_stored)),
  best  = max(abs(m$best_comp_p_new - m$best_comp_p_stored)),
  p5    = max(abs(m$p5_perm_new - m$p5_perm_stored))
)
fwrite(rec, file.path(tab_dir, "vacc_heterogeneity_comps_full.csv"))
cat("validation vs committed summary (max abs diff):\n"); print(checks)
cat("validation passed:", all(checks < 1e-10), "\nDONE\n")
