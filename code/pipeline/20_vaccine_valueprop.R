# 20_vaccine_valueprop.R
# Vaccine core value-proposition analyses (WP2 A/B/C):
#   - v1 reproduction checks (targets from results_GSE79396_tPCA)
#   - matrix vs tensor: PCA / subject-averaged PCA / tPCA / PLS-DA / tPLS-DA
#   - permutation-scheme inflation: naive row-level vs subject-level
#   - LOSO CV: naive vs subject-grouped; faithful tPLS-DA projection for the tensor model
suppressPackageStartupMessages({
  library(data.table)
  library(mixOmics)
  library(tensorOmics)
  library(MASS)
})

base_dir <- Sys.getenv("TER_ROOT", getwd())
prep <- readRDS(file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
samples <- fread(file.path(base_dir, "data/GSE79396/samples_all.csv"))
tab_dir <- file.path(base_dir, "results/tables")

Xm <- prep$Xm_scaled                       # 288 x 2000
X_tensor <- prep$X_tensor                  # 72 x 2000 x 4
subs <- prep$subs; days <- prep$days
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)],
                   levels = c("young", "elderly"))
grp_row <- factor(samples$age_group, levels = c("young", "elderly"))
cohort_row <- samples$cohort
cohort_subj <- samples$cohort[match(subs, samples$subject_id)]
sex_subj <- samples$sex[match(subs, samples$subject_id)]

set.seed(2026)
N_PERM <- 499
dct <- dctii_m_transforms(length(days))

pillai <- function(scores, g) {
  m <- manova(as.matrix(scores) ~ g)
  s <- summary(m, test = "Pillai")$stats
  c(trace = unname(s[1, "Pillai"]), p = unname(s[1, "Pr(>F)"]))
}

## ================= 0. v1 reproduction checks ================================
S_tpca <- prep$scores_tensor               # 72 x 5 (exact v1)
S_mpca <- prep$scores_matrix               # 288 x 5 (exact v1)
r0 <- list(
  tpca_manova = pillai(S_tpca[, 1:2], grp_subj),
  tpca_c1 = coef(summary(lm(S_tpca[, 1] ~ grp_subj)))["grp_subjelderly", "Pr(>|t|)"],
  tpca_c2 = coef(summary(lm(S_tpca[, 2] ~ grp_subj)))["grp_subjelderly", "Pr(>|t|)"]
)
Xsub <- do.call(rbind, lapply(subs, function(s) colMeans(Xm[samples$subject_id == s, , drop = FALSE])))
pca_sub <- mixOmics::pca(Xsub, ncomp = 2, center = TRUE, scale = FALSE)
r0$mpca_sub_manova <- pillai(pca_sub$variates$X, grp_subj)

cat("=== v1 reproduction ===\n")
cat("tPCA MANOVA Pillai:", r0$tpca_manova["trace"], " p:", r0$tpca_manova["p"],
    " | targets 0.186 / 8.09e-4\n")
cat("tPCA comp1/comp2 age LM p:", r0$tpca_c1, "/", r0$tpca_c2, " | targets 0.932 / 1.53e-4\n")
cat("subject-averaged matrix PCA Pillai:", r0$mpca_sub_manova["trace"],
    " p:", r0$mpca_sub_manova["p"], " | target 0.182 / 9.72e-4\n")

## ================= 1. model fits ============================================
fit_mpls <- mixOmics::plsda(Xm, grp_row, ncomp = 2)
S_mpls <- fit_mpls$variates$X
tp <- tpca(X_tensor, ncomp = 2)
S_tpca2 <- tp$variates
fit_tpls <- tplsda(X_tensor, as.numeric(grp_subj == "elderly"), ncomp = 2)
S_tpls <- fit_tpls$x_projected

# sanity: our manual projection reproduces the package exactly
M_all <- apply(X_tensor, c(2, 3), mean)
xhat_all <- dct$m(sweep(X_tensor, c(2, 3), M_all, "-"))
b1_check <- xhat_all[, , fit_tpls$faces[1]] %*% fit_tpls$x_loadings[, 1]
stopifnot(max(abs(b1_check - S_tpls[, 1])) < 1e-8)

## ================= 2. MANOVA table ==========================================
man_mpca_row <- pillai(S_mpca[, 1:2], grp_row)
man_mpls <- pillai(S_mpls, grp_row)
man_tpls <- pillai(S_tpls, grp_subj)
man <- data.frame(
  analysis = c("matrix PCA (288 sample-time rows)", "subject-averaged matrix PCA (72)",
               "tPCA (72 subject trajectories)", "matrix PLS-DA (288 rows)",
               "tPLS-DA (72 subject trajectories)"),
  statistic = c(man_mpca_row["trace"], r0$mpca_sub_manova["trace"], r0$tpca_manova["trace"],
                man_mpls["trace"], man_tpls["trace"]),
  p_parametric = c(man_mpca_row["p"], r0$mpca_sub_manova["p"], r0$tpca_manova["p"],
                   man_mpls["p"], man_tpls["p"]))

## ================= 3. permutation inflation =================================
row_labels_of <- function(gs) gs[match(samples$subject_id, subs)]
p_mpca_row <- (1 + sum(sapply(1:N_PERM, function(b) pillai(S_mpca[, 1:2], sample(grp_row))["trace"]) >= man_mpca_row["trace"])) / (1 + N_PERM)
p_mpca_sub <- (1 + sum(sapply(1:N_PERM, function(b) pillai(S_mpca[, 1:2], row_labels_of(sample(grp_subj)))["trace"]) >= man_mpca_row["trace"])) / (1 + N_PERM)
p_sub_row <- (1 + sum(sapply(1:N_PERM, function(b) pillai(pca_sub$variates$X, sample(grp_subj))["trace"]) >= r0$mpca_sub_manova["trace"])) / (1 + N_PERM)
p_tpca_sub <- (1 + sum(sapply(1:N_PERM, function(b) pillai(S_tpca2, sample(grp_subj))["trace"]) >= r0$tpca_manova["trace"])) / (1 + N_PERM)

set.seed(101)
pm_row <- numeric(N_PERM); pm_sub <- numeric(N_PERM)
for (b in 1:N_PERM) {
  g1 <- sample(grp_row)
  f1 <- mixOmics::plsda(Xm, g1, ncomp = 2)
  pm_row[b] <- pillai(f1$variates$X, g1)["trace"]
  gs <- sample(grp_subj); g2 <- row_labels_of(gs)
  f2 <- mixOmics::plsda(Xm, g2, ncomp = 2)
  pm_sub[b] <- pillai(f2$variates$X, g2)["trace"]
}
p_mpls_row <- (1 + sum(pm_row >= man_mpls["trace"])) / (1 + N_PERM)
p_mpls_sub <- (1 + sum(pm_sub >= man_mpls["trace"])) / (1 + N_PERM)
set.seed(102)
pt <- numeric(N_PERM)
for (b in 1:N_PERM) {
  g <- sample(grp_subj)
  f <- tplsda(X_tensor, as.numeric(g == "elderly"), ncomp = 2)
  pt[b] <- pillai(f$x_projected, g)["trace"]
}
p_tpls_sub <- (1 + sum(pt >= man_tpls["trace"])) / (1 + N_PERM)

man$p_perm_subject <- c(p_mpca_sub, p_sub_row, p_tpca_sub, p_mpls_sub, p_tpls_sub)
perm_tab <- data.frame(
  analysis = c("matrix PCA", "matrix PCA", "subject-averaged matrix PCA",
               "matrix PLS-DA", "matrix PLS-DA", "tPLS-DA"),
  scheme = c("row-level (naive)", "subject-level", "subject-level",
             "row-level (naive)", "subject-level", "subject-level"),
  n_perm = N_PERM,
  p_perm = c(p_mpca_row, p_mpca_sub, p_sub_row, p_mpls_row, p_mpls_sub, p_tpls_sub))

## ==== 3.5 fold-safe helpers (P2B-verified) ======================================
dat_vac <- readRDS(file.path(base_dir, "data/GSE79396/vaccine_data.rds"))
M_all <- t(dat_vac$expr[, samples$gsm, drop = FALSE])   # 288 x 54715 raw log2
TOP_N <- 2000L
fold_select_z <- function(M, tr, te) {
  v <- apply(M[tr, , drop = FALSE], 2, var)
  sel <- order(v, decreasing = TRUE)[seq_len(min(TOP_N, ncol(M)))]
  mu <- colMeans(M[tr, sel, drop = FALSE]); sdv <- apply(M[tr, sel, drop = FALSE], 2, sd)
  list(sel = sel,
       Ztr = sweep(sweep(M[tr, sel, drop = FALSE], 2, mu, "-"), 2, sdv, "/"),
       Zte = sweep(sweep(M[te, sel, drop = FALSE], 2, mu, "-"), 2, sdv, "/"))
}
build_tensor_fs <- function(Zrows, subj_v, day_v, subs_v, days_v) {
  Xt <- array(NA_real_, dim = c(length(subs_v), ncol(Zrows), length(days_v)),
              dimnames = list(subs_v, colnames(Zrows), paste0("D", days_v)))
  key <- paste(subj_v, day_v, sep = "|")
  for (s in seq_along(subs_v)) for (d in seq_along(days_v)) {
    i <- match(paste(subs_v[s], days_v[d], sep = "|"), key)
    if (is.na(i)) stop("build_tensor_fs: missing sample")
    Xt[s, , d] <- Zrows[i, ]
  }
  Xt
}

## ================= 4. LOSO CV ================================================
## 4a. matrix PLS-DA: naive row-level vs subject-grouped
pred_class_m <- function(fit, Xnew) {
  pr <- predict(fit, newdata = Xnew)
  cls <- pr$class$max.dist
  pr$clsvec <- as.character(cls[, ncol(cls)])
  pr
}
pred_row <- character(nrow(Xm))
for (i in 1:nrow(Xm)) {
  f <- mixOmics::plsda(Xm[-i, ], grp_row[-i], ncomp = 2)
  pred_row[i] <- pred_class_m(f, Xm[i, , drop = FALSE])$clsvec
}
ba_row <- mean(sapply(levels(grp_subj), function(cl) mean(pred_row[grp_row == cl] == cl)))

pred_subj_m <- character(length(subs))
for (s in seq_along(subs)) {
  hold <- which(samples$subject_id == subs[s])
  fs <- fold_select_z(M_all, setdiff(seq_len(nrow(M_all)), hold), hold)
  f <- mixOmics::plsda(fs$Ztr, grp_row[-hold], ncomp = 2)
  pr <- pred_class_m(f, fs$Zte)
  tb <- table(pr$clsvec)
  if (length(tb) == 1 || max(tb) > sum(tb) / 2) {
    pred_subj_m[s] <- names(which.max(tb))
  } else {
    prob <- pr$predict[, , dim(pr$predict)[3], drop = FALSE]
    mv <- colMeans(prob[, , 1])
    pred_subj_m[s] <- names(which.max(mv))
  }
}
ba_group_m <- mean(sapply(levels(grp_subj), function(cl) mean(pred_subj_m[grp_subj == cl] == cl)))

## 4b. tPLS-DA: faithful subject-grouped LOSO via manual projection
tplsda_fold <- function(Xtr, ytr, Xnew) {
  fit <- tplsda(Xtr, ytr, ncomp = 2)
  Mt <- apply(Xtr, c(2, 3), mean)
  xh <- dct$m(sweep(Xtr, c(2, 3), Mt, "-"))
  a <- fit$x_loadings; f <- fit$faces
  b <- matrix(NA_real_, dim(Xtr)[1], 2); regs <- vector("list", 2)
  xst <- xh
  for (i in 1:2) {
    bi <- xst[, , f[i]] %*% a[, i]
    b[, i] <- bi
    ri <- crossprod(xst[, , f[i]], bi) / as.numeric(crossprod(bi))
    regs[[i]] <- as.numeric(ri)
    xst[, , f[i]] <- xst[, , f[i]] - tcrossprod(bi, ri)
  }
  xh_new <- dct$m(sweep(Xnew, c(2, 3), Mt, "-"))
  bnew <- numeric(2)
  for (i in 1:2) {
    bnew[i] <- as.numeric(xh_new[, , f[i]] %*% a[, i])
    xh_new[, , f[i]] <- xh_new[, , f[i]] - tcrossprod(bnew[i], regs[[i]])
  }
  list(b = b, bnew = bnew, fit = fit)
}
pred_subj_t <- character(length(subs))
for (s in seq_along(subs)) {
  tr <- setdiff(seq_along(subs), s)
  hold_rows <- which(samples$subject_id == subs[s])
  fs <- fold_select_z(M_all, setdiff(seq_len(nrow(M_all)), hold_rows), hold_rows)
  colnames(fs$Ztr) <- colnames(fs$Zte) <- colnames(M_all)[fs$sel]
  Xtr_fs <- build_tensor_fs(fs$Ztr, samples$subject_id[-hold_rows], samples$day[-hold_rows], subs[tr], days)
  Xte_fs <- build_tensor_fs(fs$Zte, samples$subject_id[hold_rows], samples$day[hold_rows], subs[s], days)
  res <- tplsda_fold(Xtr_fs,
                     as.numeric(grp_subj[tr] == "elderly"),
                     Xte_fs)
  ld <- MASS::lda(res$b, grouping = grp_subj[tr])
  pred_subj_t[s] <- as.character(predict(ld, newdata = matrix(res$bnew, 1, 2))$class)
}
ba_group_t <- mean(sapply(levels(grp_subj), function(cl) mean(pred_subj_t[grp_subj == cl] == cl)))

loso_tab <- data.frame(
  analysis = c("matrix PLS-DA, row-level LOSO (naive)", "matrix PLS-DA, subject-grouped LOSO",
               "tPLS-DA, subject-grouped LOSO"),
  balanced_accuracy = c(ba_row, ba_group_m, ba_group_t),
  n_predictions = c(nrow(Xm), length(subs), length(subs)))

## ================= 5. save ===================================================
fwrite(man, file.path(tab_dir, "vacc_manova.csv"))
fwrite(perm_tab, file.path(tab_dir, "vacc_perm_inflation.csv"))
fwrite(loso_tab, file.path(tab_dir, "vacc_loso.csv"))

scores_subj <- data.frame(
  subject = subs, age_group = grp_subj, cohort = cohort_subj, sex = sex_subj,
  mPCA1 = S_mpca[match(subs, samples$subject_id), 1],
  mPCA2 = S_mpca[match(subs, samples$subject_id), 2],
  sPCA1 = pca_sub$variates$X[, 1], sPCA2 = pca_sub$variates$X[, 2],
  tPCA1 = S_tpca2[, 1], tPCA2 = S_tpca2[, 2],
  mPLS1 = S_mpls[match(subs, samples$subject_id), 1],
  mPLS2 = S_mpls[match(subs, samples$subject_id), 2],
  tPLS1 = S_tpls[, 1], tPLS2 = S_tpls[, 2],
  pred_mpls_grouped = pred_subj_m, pred_tpls_grouped = pred_subj_t)
fwrite(scores_subj, file.path(tab_dir, "vacc_scores_subject.csv"))

scores_row <- data.frame(
  gsm = samples$gsm, subject = samples$subject_id, day = samples$day,
  age_group = samples$age_group, cohort = samples$cohort,
  mPCA1 = S_mpca[, 1], mPCA2 = S_mpca[, 2],
  mPLS1 = S_mpls[, 1], mPLS2 = S_mpls[, 2],
  pred_mpls_naive = pred_row)
fwrite(scores_row, file.path(tab_dir, "vacc_scores_row.csv"))

# loadings for downstream biological interpretation
ann <- fread(file.path(base_dir, "data/GPL13158/annot_clean.csv"))
sel <- prep$sel
loadings_tpls <- data.frame(probe = sel, symbol = ann$symbol[match(sel, ann$probe)],
                            tPLS1 = fit_tpls$x_loadings[, 1], tPLS2 = fit_tpls$x_loadings[, 2],
                            tPCA1 = tp$loadings[, 1], tPCA2 = tp$loadings[, 2])
fwrite(loadings_tpls, file.path(tab_dir, "vacc_loadings_tensor.csv"))

cat("\n=== MANOVA ===\n"); print(man)
cat("\n=== permutation inflation ===\n"); print(perm_tab)
cat("\n=== LOSO ===\n"); print(loso_tab)
cat("\nDONE\n")
