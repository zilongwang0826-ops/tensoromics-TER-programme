# 07_fish_comparison.R
# Fish comparison analysis for the revised report:
#   - matrix PCA / matrix PLS-DA / tPCA / tPLS-DA on the same 113-fish data
#   - permutation-scheme inflation (naive row-level vs subject-level)
#   - LOSO CV: naive row-level vs subject-grouped (matrix PLS-DA)
# Uses the official tensorOmics package for the tensor models.
suppressPackageStartupMessages({
  library(data.table)
  library(mixOmics)
  library(tensorOmics)
})

base_dir <- Sys.getenv("TER_ROOT", getwd())
dat_dir <- file.path(base_dir, "data/GSE150318")
tab_dir <- file.path(base_dir, "results/tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1)
N_PERM <- 499
N_BOOT <- NULL

X <- fread(file.path(dat_dir, "fish_X1000_long.csv"))
meta <- fread(file.path(dat_dir, "fish_meta.csv"))
genes <- setdiff(names(X), c("subj", "time"))
subj <- X$subj; time <- X$time
n_subj <- nrow(meta)
cat("rows:", nrow(X), "| genes:", length(genes), "| subjects:", n_subj, "\n")

Xm <- as.matrix(X[, ..genes])          # 226 x 1000 (subject-time rows)
rownames(Xm) <- paste0("s", subj, "_t", time)
grp_row <- factor(meta$long_lived[match(subj, meta$subj)],
                  levels = c(0, 1), labels = c("short_lived", "long_lived"))
grp_subj <- factor(meta$long_lived, levels = c(0, 1), labels = c("short_lived", "long_lived"))

# tensor: 113 x 1000 x 2 (order [10w, 20w] as exported)
X_tensor <- array(NA_real_, dim = c(n_subj, length(genes), 2),
                  dimnames = list(paste0("s", meta$subj), genes, c("10w", "20w")))
for (i in seq_len(n_subj)) for (k in 1:2) {
  X_tensor[i, , k] <- Xm[paste0("s", meta$subj[i], "_t", c(10, 20)[k]), ]
}
stopifnot(!anyNA(X_tensor))

pillai <- function(scores, g) {
  # MANOVA Pillai's trace + approx p; scores: n x q matrix/data.frame
  m <- manova(as.matrix(scores) ~ g)
  s <- summary(m, test = "Pillai")$stats
  c(trace = unname(s[1, "Pillai"]), p = unname(s[1, "Pr(>F)"]))
}

## ---------------- unsupervised: matrix PCA vs tPCA --------------------------
pca_m <- mixOmics::pca(Xm, ncomp = 2, center = TRUE, scale = FALSE)
S_pca <- pca_m$variates$X
tp <- tpca(X_tensor, ncomp = 2)
S_tpca <- tp$variates

res_man <- list()
res_man$m_pca <- pillai(S_pca, grp_row)
res_man$t_pca <- pillai(S_tpca, grp_subj)

## permutation for unsupervised scores: labels permuted, scores fixed
# subject-level permutation for tensor scores
perm_tpca <- sapply(1:N_PERM, function(b) pillai(S_tpca, sample(grp_subj))["trace"])
# matrix PCA: row-level labels (naive) and subject-level labels
perm_mpca_row <- sapply(1:N_PERM, function(b) pillai(S_pca, sample(grp_row))["trace"])
perm_mpca_sub <- sapply(1:N_PERM, function(b) {
  gs <- sample(grp_subj)
  pillai(S_pca, gs[match(subj, meta$subj)])["trace"]
})
p_mpca_row <- (1 + sum(perm_mpca_row >= res_man$m_pca["trace"])) / (1 + N_PERM)
p_mpca_sub <- (1 + sum(perm_mpca_sub >= res_man$m_pca["trace"])) / (1 + N_PERM)
p_tpca_sub <- (1 + sum(perm_tpca >= res_man$t_pca["trace"])) / (1 + N_PERM)

## ---------------- supervised: matrix PLS-DA vs tPLS-DA ----------------------
fit_m <- mixOmics::plsda(Xm, grp_row, ncomp = 2)
S_mpls <- fit_m$variates$X
res_man$m_plsda <- pillai(S_mpls, grp_row)

fit_t <- tplsda(X_tensor, as.numeric(grp_subj == "long_lived"), ncomp = 2)
S_tpls <- fit_t$x_projected
res_man$t_plsda <- pillai(S_tpls, grp_subj)

# supervised permutation requires refitting under the permuted labels
perm_sup <- function(kind, scheme, n = N_PERM) {
  tr <- numeric(n)
  if (kind == "m") {
    for (b in 1:n) {
      g <- if (scheme == "row") sample(grp_row) else {
        gs <- sample(grp_subj); gs[match(subj, meta$subj)]
      }
      f <- mixOmics::plsda(Xm, g, ncomp = 2)
      tr[b] <- pillai(f$variates$X, g)["trace"]
    }
  } else {
    for (b in 1:n) {
      g <- sample(grp_subj)
      f <- tplsda(X_tensor, as.numeric(g == "long_lived"), ncomp = 2)
      tr[b] <- pillai(f$x_projected, g)["trace"]
    }
  }
  tr
}
set.seed(2)
pt <- perm_sup("t", "subj")
set.seed(3)
pm_row <- perm_sup("m", "row")
set.seed(4)
pm_sub <- perm_sup("m", "subj")
p_mpls_row <- (1 + sum(pm_row >= res_man$m_plsda["trace"])) / (1 + N_PERM)
p_mpls_sub <- (1 + sum(pm_sub >= res_man$m_plsda["trace"])) / (1 + N_PERM)
p_tpls_sub <- (1 + sum(pt >= res_man$t_plsda["trace"])) / (1 + N_PERM)

## ---------------- LOSO CV: naive vs grouped (matrix PLS-DA) -----------------
latent_centroid_acc <- function(fit, Xnew, gnew) {
  pred <- predict(fit, newdata = Xnew)
  cls <- pred$class$max.dist
  as.character(cls[, ncol(cls)])
}

set.seed(5)
# naive row-level LOSO
pred_row <- character(nrow(Xm))
for (i in 1:nrow(Xm)) {
  f <- mixOmics::plsda(Xm[-i, ], grp_row[-i], ncomp = 2)
  pred_row[i] <- latent_centroid_acc(f, Xm[i, , drop = FALSE], grp_row)
}
# confusion and balanced accuracy
ba_row <- mean(sapply(levels(grp_row), function(cl) {
  mean(pred_row[grp_row == cl] == cl)
}))

# subject-grouped LOSO
pred_subj <- character(n_subj)
for (s in 1:n_subj) {
  hold <- which(subj == meta$subj[s])
  f <- mixOmics::plsda(Xm[-hold, ], grp_row[-hold], ncomp = 2)
  pr <- latent_centroid_acc(f, Xm[hold, , drop = FALSE], grp_subj)
  pred_subj[s] <- names(sort(table(pr), decreasing = TRUE))[1]  # majority of 2 rows
}
ba_group <- mean(sapply(levels(grp_subj), function(cl) {
  mean(pred_subj[grp_subj == cl] == cl)
}))

## ---------------- save --------------------------------------------------------
man <- rbind(
  data.frame(analysis = c("matrix PCA (row labels)", "tPCA (subject labels)",
                          "matrix PLS-DA", "tPLS-DA"),
             statistic = c(res_man$m_pca["trace"], res_man$t_pca["trace"],
                           res_man$m_plsda["trace"], res_man$t_plsda["trace"]),
             p_parametric = c(res_man$m_pca["p"], res_man$t_pca["p"],
                              res_man$m_plsda["p"], res_man$t_plsda["p"]))
)
fwrite(man, file.path(tab_dir, "fish_manova.csv"))

perm_tab <- data.frame(
  analysis = c("matrix PCA", "matrix PCA", "tPCA",
               "matrix PLS-DA", "matrix PLS-DA", "tPLS-DA"),
  scheme = c("row-level (naive)", "subject-level", "subject-level",
             "row-level (naive)", "subject-level", "subject-level"),
  n_perm = N_PERM,
  p_perm = c(p_mpca_row, p_mpca_sub, p_tpca_sub,
             p_mpls_row, p_mpls_sub, p_tpls_sub)
)
fwrite(perm_tab, file.path(tab_dir, "fish_perm_inflation.csv"))

loso_tab <- data.frame(
  analysis = c("matrix PLS-DA, row-level LOSO (naive)",
               "matrix PLS-DA, subject-grouped LOSO"),
  balanced_accuracy = c(ba_row, ba_group),
  n_predictions = c(nrow(Xm), n_subj)
)
fwrite(loso_tab, file.path(tab_dir, "fish_loso_naive_vs_grouped.csv"))

scores <- data.frame(
  subj = meta$subj, age_death = meta$age_death, long_lived = meta$long_lived,
  mPCA1 = S_pca[, 1], mPCA2 = S_pca[, 2],
  tPCA1 = S_tpca[, 1], tPCA2 = S_tpca[, 2],
  mPLS1 = S_mpls[, 1], mPLS2 = S_mpls[, 2],
  tPLS1 = S_tpls[, 1], tPLS2 = S_tpls[, 2]
)
fwrite(scores, file.path(dat_dir, "fish_scores_4panel.csv"))

cat("\n=== fish MANOVA ===\n"); print(man)
cat("\n=== permutation inflation ===\n"); print(perm_tab)
cat("\n=== LOSO naive vs grouped ===\n"); print(loso_tab)
cat("\nDONE\n")
