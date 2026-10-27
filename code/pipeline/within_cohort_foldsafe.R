#!/usr/bin/env Rscript
# within_cohort_foldsafe.R — full fold-safe re-derivation of the within-cohort LOSO metrics
# (Denver / Atlanta) that are cited in report Section 5.3, plus in-process leaky parity runs.
# For each cohort subject fold: top-2000 HVGs selected and z-scored on the cohort's OTHER
# subjects only (fold-safe), versus the committed definition (global selection + global z tensor).
suppressPackageStartupMessages({ library(data.table); library(tensorOmics); library(MASS) })
PROJ <- Sys.getenv("TER_ROOT", getwd())
source(file.path(PROJ, "scripts/lib_tensor_helpers.R"))
dat     <- readRDS(file.path(PROJ, "data/GSE79396/vaccine_data.rds"))
prep    <- readRDS(file.path(PROJ, "data/GSE79396/preprocessed_A.rds"))
samples <- fread(file.path(PROJ, "data/GSE79396/samples_all.csv"))
tab_dir <- file.path(PROJ, "results/tables")
subs <- prep$subs; days <- prep$days
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
cohort_subj <- samples$cohort[match(subs, samples$subject_id)]
expr <- dat$expr
M_all <- t(expr[, samples$gsm, drop = FALSE])
TOP_N <- 2000L; dct <- make_dct(length(days))

build_tensor <- function(Zrows, subj_v, day_v, subs_v, days_v) {
  Xt <- array(NA_real_, dim = c(length(subs_v), ncol(Zrows), length(days_v)),
              dimnames = list(subs_v, colnames(Zrows), paste0("D", days_v)))
  key <- paste(subj_v, day_v, sep = "|")
  for (s in seq_along(subs_v)) for (d in seq_along(days_v)) {
    i <- match(paste(subs_v[s], days_v[d], sep = "|"), key)
    if (is.na(i)) stop("build_tensor: missing sample")
    Xt[s, , d] <- Zrows[i, ]
  }
  Xt
}
mat_agg <- function(pr) {
  cls <- as.character(pr$class$max.dist[, ncol(pr$class$max.dist)])
  tb <- table(cls)
  if (length(tb) == 1 || max(tb) > sum(tb) / 2) return(names(which.max(tb)))
  prob <- pr$predict[, , dim(pr$predict)[3], drop = FALSE]
  names(which.max(colMeans(prob[, , 1])))
}
bacc <- function(pred, grp) mean(sapply(levels(grp), function(cl) mean(pred[grp == cl] == cl)))

rows <- list(); pred_rows <- list()
for (co in levels(factor(cohort_subj))) {
  idx <- which(cohort_subj == co)
  g <- droplevels(grp_subj[idx])
  Xm_c_raw <- as.matrix(prep$Xm_scaled)[samples$subject_id %in% subs[idx], ]
  Xt_glob <- prep$X_tensor[idx, , , drop = FALSE]
  pr_lk_m <- character(length(idx)); pr_lk_t <- character(length(idx))
  pr_fs_m <- character(length(idx)); pr_fs_t <- character(length(idx))
  for (s in seq_along(idx)) {
    subj_s <- subs[idx[s]]
    hold   <- which(samples$subject_id == subj_s)
    tr_rows <- setdiff(seq_len(nrow(M_all)), hold)
    tr_subs <- idx[-s]
    ytr <- as.numeric(grp_subj[tr_subs] == "elderly")
    ## leaky (committed): global selection matrix + global-z tensor
    rows_c <- which(samples$subject_id %in% subs[idx])
    hloc   <- which(subs[idx] == subj_s)
    f_lk <- mixOmics::plsda(Xm_c_raw[samples$subject_id[rows_c] != subj_s, , drop = FALSE],
                            factor(samples$age_group[rows_c][samples$subject_id[rows_c] != subj_s],
                                   levels = c("young", "elderly")), ncomp = 2)
    pr_lk_m[s] <- mat_agg(predict(f_lk, newdata = Xm_c_raw[samples$subject_id[rows_c] == subj_s, , drop = FALSE]))
    prt_lk <- tplsda_project(Xt_glob[-hloc, , , drop = FALSE], ytr, Xt_glob[hloc, , , drop = FALSE], dct, ncomp = 2)
    ld_lk <- MASS::lda(prt_lk$b, grouping = factor(grp_subj[tr_subs]))
    pr_lk_t[s] <- as.character(predict(ld_lk, newdata = matrix(prt_lk$bnew, 1, 2))$class)
    ## fold-safe: selection + z from cohort training rows only
    v <- apply(M_all[tr_rows, , drop = FALSE], 2, var)
    sel <- order(v, decreasing = TRUE)[seq_len(min(TOP_N, ncol(M_all)))]
    mu <- colMeans(M_all[tr_rows, sel, drop = FALSE]); sdv <- apply(M_all[tr_rows, sel, drop = FALSE], 2, sd)
    Ztr <- sweep(sweep(M_all[tr_rows, sel, drop = FALSE], 2, mu, "-"), 2, sdv, "/")
    Zte <- sweep(sweep(M_all[hold, sel, drop = FALSE], 2, mu, "-"), 2, sdv, "/")
    colnames(Ztr) <- colnames(Zte) <- colnames(M_all)[sel]
    f_fs <- mixOmics::plsda(Ztr, factor(samples$age_group[tr_rows], levels = c("young", "elderly")), ncomp = 2)
    pr_fs_m[s] <- mat_agg(predict(f_fs, newdata = Zte))
    Xtr <- build_tensor(Ztr, samples$subject_id[tr_rows], samples$day[tr_rows], subs[tr_subs], days)
    Xte <- build_tensor(Zte, samples$subject_id[hold], samples$day[hold], subs[idx[s]], days)
    prt_fs <- tplsda_project(Xtr, ytr, Xte, dct, ncomp = 2)
    ld_fs <- MASS::lda(prt_fs$b, grouping = factor(grp_subj[tr_subs]))
    pr_fs_t[s] <- as.character(predict(ld_fs, newdata = matrix(prt_fs$bnew, 1, 2))$class)
  }
  rows[[length(rows) + 1]] <- data.frame(
    cohort = co, n = length(idx),
    leaky_matrix_acc = bacc(pr_lk_m, factor(grp_subj[idx])),
    foldsafe_matrix_acc = bacc(pr_fs_m, factor(grp_subj[idx])),
    leaky_tensor_acc = bacc(pr_lk_t, factor(grp_subj[idx])),
    foldsafe_tensor_acc = bacc(pr_fs_t, factor(grp_subj[idx])))
  pred_rows[[length(pred_rows) + 1]] <- data.frame(
    cohort = co, subject = subs[idx], grp = as.character(grp_subj[idx]),
    leaky_matrix = pr_lk_m, fs_matrix = pr_fs_m, leaky_tensor = pr_lk_t, fs_tensor = pr_fs_t)
  cat(sprintf("[within-fs] %s done (n=%d)\n", co, length(idx)))
}
res <- rbindlist(rows); preds <- rbindlist(pred_rows)
res[, `:=`(delta_matrix = foldsafe_matrix_acc - leaky_matrix_acc,
           delta_tensor = foldsafe_tensor_acc - leaky_tensor_acc)]
fwrite(res, file.path(tab_dir, "vacc_within_cohort_foldsafe.csv"))
fwrite(preds, file.path(tab_dir, "vacc_within_cohort_foldsafe_predictions.csv"))
print(res)
cat("stored (committed): Atlanta 0.379464 / 0.410714 ; Denver 0.866667 / 0.791667\n")
cat("DONE\n")

