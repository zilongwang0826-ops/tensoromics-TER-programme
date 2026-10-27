# 21_vaccine_batch.R
# WP2D: cohort (site) detection and correction workflow, mirroring tensorOmics paper Appendix B.
#   - tPCA detection of cohort structure (before)
#   - ComBat (primary; preserves age/sex) and PLSDA_batch (paper-mirror)
#   - tPCA/tPLS-DA after correction, age association before vs after
#   - within-cohort tPLS-DA robustness (Atlanta / Denver separately)
suppressPackageStartupMessages({
  library(data.table)
  library(mixOmics)
  library(tensorOmics)
  library(sva)
  library(PLSDAbatch)
})
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_tensor_helpers.R"))

base_dir <- Sys.getenv("TER_ROOT", getwd())
prep <- readRDS(file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
dat <- readRDS(file.path(base_dir, "data/GSE79396/vaccine_data.rds"))
samples <- fread(file.path(base_dir, "data/GSE79396/samples_all.csv"))
tab_dir <- file.path(base_dir, "results/tables")
set.seed(2026)
N_PERM <- 299

sel <- prep$sel
days <- prep$days; subs <- prep$subs
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
grp_row <- factor(samples$age_group, levels = c("young", "elderly"))
cohort_row <- factor(samples$cohort)
cohort_subj <- factor(samples$cohort[match(subs, samples$subject_id)])
sex_row <- factor(samples$sex)
dct <- make_dct(length(days))

raw2000 <- t(dat$expr[sel, ])   # 288 x 2000 raw log2
stopifnot(rownames(raw2000) == samples$gsm)

build_tensor_scaled <- function(Xrow) {
  Z <- scale(Xrow)
  Xt <- array(NA_real_, dim = c(length(subs), ncol(Z), length(days)),
              dimnames = list(subs, colnames(Z), paste0("D", days)))
  for (s in subs) for (d in days) {
    Xt[s, , paste0("D", d)] <- Z[samples$subject_id == s & samples$day == d, ]
  }
  Xt
}

## ---- correction methods -----------------------------------------------------
mod <- model.matrix(~ grp_row + sex_row)
combat <- ComBat(dat = t(raw2000), batch = cohort_row, mod = mod)
combat_rows <- t(combat)
pb <- PLSDA_batch(X = combat_rows, Y.trt = grp_row, Y.bat = cohort_row,
                  ncomp.trt = 2, ncomp.bat = 2, mode = "regression")
plsdb_rows <- pb$X.nobatch
stopifnot(nrow(plsdb_rows) == 288)

Xt_list <- list(
  before    = prep$X_tensor,
  combat    = build_tensor_scaled(combat_rows),
  plsdabatch = build_tensor_scaled(plsdb_rows)
)

## ---- evaluate each version --------------------------------------------------
evaluate_tensor <- function(Xt) {
  tp <- tpca(Xt, ncomp = 2)
  ft <- tplsda(Xt, as.numeric(grp_subj == "elderly"), ncomp = 2)
  list(tp = tp$variates, tpls = ft$x_projected, fit = ft)
}
perm_p_unsup <- function(S, g, n = N_PERM) {
  obs <- pillai(S, g)["trace"]
  (1 + sum(sapply(seq_len(n), function(b) pillai(S, sample(g))["trace"]) >= obs)) / (1 + n)
}
perm_p_sup <- function(Xt, g, n = N_PERM) {
  obs <- pillai(tplsda(Xt, as.numeric(g == "elderly"), ncomp = 2)$x_projected, g)["trace"]
  (1 + sum(sapply(seq_len(n), function(b) {
    gg <- sample(g)
    pillai(tplsda(Xt, as.numeric(gg == "elderly"), ncomp = 2)$x_projected, gg)["trace"]
  }) >= obs)) / (1 + n)
}

set.seed(7)
eva <- lapply(Xt_list, evaluate_tensor)
res <- data.frame()
for (nm in names(eva)) {
  e <- eva[[nm]]
  cohort_tp <- pillai(e$tp, cohort_subj); age_tp <- pillai(e$tp, grp_subj)
  age_tpls <- pillai(e$tpls, grp_subj)
  res <- rbind(res, data.frame(
    version = nm,
    cohort_pillai_tPCA = cohort_tp["trace"], cohort_p_tPCA = cohort_tp["p"],
    age_pillai_tPCA = age_tp["trace"], age_p_tPCA = age_tp["p"],
    age_pillai_tPLSDA = age_tpls["trace"], age_p_tPLSDA = age_tpls["p"],
    eta2_cohort_tPCA1 = summary(lm(e$tp[, 1] ~ cohort_subj))$r.squared,
    eta2_cohort_tPCA2 = summary(lm(e$tp[, 2] ~ cohort_subj))$r.squared,
    eta2_age_tPLS1 = summary(lm(e$tpls[, 1] ~ grp_subj))$r.squared
  ))
}
# permutation p for the key contrasts (before + combat + plsdabatch)
res$p_perm_cohort_tPCA <- NA
res$p_perm_age_tPCA <- NA
res$p_perm_age_tPLSDA <- NA
for (nm in names(eva)) {
  res$p_perm_cohort_tPCA[res$version == nm] <- perm_p_unsup(eva[[nm]]$tp, cohort_subj)
  res$p_perm_age_tPCA[res$version == nm] <- perm_p_unsup(eva[[nm]]$tp, grp_subj)
  res$p_perm_age_tPLSDA[res$version == nm] <- perm_p_sup(Xt_list[[nm]], grp_subj)
}

## ---- within-cohort robustness (before correction) ---------------------------
pred_matrix_grouped <- function(Xm_sub, subj_sub, grp) {
  pred <- character(length(subj_sub))
  for (s in seq_along(subj_sub)) {
    hold <- which(subj_sub == subj_sub[s])
    f <- mixOmics::plsda(Xm_sub[-hold, ], grp[-hold], ncomp = 2)
    pr <- predict(f, newdata = Xm_sub[hold, , drop = FALSE])
    tb <- table(pr$class$max.dist[, 2])
    if (length(tb) == 1 || max(tb) > sum(tb) / 2) pred[s] <- names(which.max(tb)) else {
      mv <- colMeans(pr$predict[, , dim(pr$predict)[3]])
      pred[s] <- names(which.max(mv))
    }
  }
  pred
}
within <- data.frame()
Xm_full <- prep$Xm_scaled
for (co in levels(cohort_subj)) {
  idx <- which(cohort_subj == co)
  g <- droplevels(grp_subj[idx])
  Xt_c <- Xt_list$before[idx, , , drop = FALSE]
  man <- pillai(tpca(Xt_c, ncomp = 2)$variates, g)
  perm <- perm_p_unsup(tpca(Xt_c, ncomp = 2)$variates, g, n = 199)
  pred_t <- tplsda_loo_predict(Xt_c, g, dct, ncomp = 2)
  ba_t <- balanced_acc(pred_t, g)
  rows_c <- which(samples$subject_id %in% subs[idx])
  pred_m <- pred_matrix_grouped(Xm_full[rows_c, ], samples$subject_id[rows_c],
                                factor(samples$age_group[rows_c], levels = c("young", "elderly")))
  ba_m <- balanced_acc(pred_m, factor(samples$age_group[rows_c], levels = c("young", "elderly")))
  within <- rbind(within, data.frame(
    cohort = co, n = length(idx), n_young = sum(g == "young"), n_elderly = sum(g == "elderly"),
    tPCA_pillai = man["trace"], tPCA_p_parametric = man["p"], tPCA_p_perm = perm,
    tPLSDA_loo_balanced_acc = ba_t, matrixPLSDA_loo_balanced_acc = ba_m))
}

## ---- save -------------------------------------------------------------------
fwrite(res, file.path(tab_dir, "vacc_batch_metrics.csv"))
fwrite(within, file.path(tab_dir, "vacc_within_cohort.csv"))

scores <- data.frame(subject = subs, cohort = cohort_subj, age_group = grp_subj,
                     before_tPCA1 = eva$before$tp[, 1], before_tPCA2 = eva$before$tp[, 2],
                     before_tPLS1 = eva$before$tpls[, 1], before_tPLS2 = eva$before$tpls[, 2],
                     combat_tPCA1 = eva$combat$tp[, 1], combat_tPCA2 = eva$combat$tp[, 2],
                     combat_tPLS1 = eva$combat$tpls[, 1], combat_tPLS2 = eva$combat$tpls[, 2],
                     plsdb_tPCA1 = eva$plsdabatch$tp[, 1], plsdb_tPCA2 = eva$plsdabatch$tp[, 2],
                     plsdb_tPLS1 = eva$plsdabatch$tpls[, 1], plsdb_tPLS2 = eva$plsdabatch$tpls[, 2])
fwrite(scores, file.path(tab_dir, "vacc_batch_scores.csv"))

saveRDS(list(combat_rows = combat_rows, plsdb_rows = plsdb_rows),
        file.path(base_dir, "data/GSE79396/corrected_versions.rds"))

cat("\n=== batch correction metrics ===\n"); print(res)
cat("\n=== within-cohort ===\n"); print(within)
cat("\nDONE\n")
