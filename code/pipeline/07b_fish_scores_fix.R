# 07b_fish_scores_fix.R : rebuild fish score tables with correct row alignment
suppressPackageStartupMessages({library(data.table); library(mixOmics); library(tensorOmics)})
base_dir <- Sys.getenv("TER_ROOT", getwd())
dat_dir <- file.path(base_dir, "data/GSE150318")

X <- fread(file.path(dat_dir, "fish_X1000_long.csv"))
meta <- fread(file.path(dat_dir, "fish_meta.csv"))
genes <- setdiff(names(X), c("subj", "time"))
subj <- X$subj; time <- X$time
Xm <- as.matrix(X[, ..genes])
rownames(Xm) <- paste0("s", subj, "_t", time)
grp_row <- factor(meta$long_lived[match(subj, meta$subj)], levels = c(0, 1),
                  labels = c("short_lived", "long_lived"))
grp_subj <- factor(meta$long_lived, levels = c(0, 1), labels = c("short_lived", "long_lived"))

pca_m <- mixOmics::pca(Xm, ncomp = 2, center = TRUE, scale = FALSE)
S_pca <- pca_m$variates$X
mpls <- mixOmics::plsda(Xm, grp_row, ncomp = 2)
S_mpls <- mpls$variates$X

n <- nrow(meta)
X_tensor <- array(NA_real_, dim = c(n, length(genes), 2), dimnames = list(paste0("s", meta$subj), genes, c("10w", "20w")))
for (i in 1:n) for (k in 1:2) X_tensor[i, , k] <- Xm[paste0("s", meta$subj[i], "_t", c(10, 20)[k]), ]
tp <- tpca(X_tensor, ncomp = 2)
ft <- tplsda(X_tensor, as.numeric(grp_subj == "long_lived"), ncomp = 2)

row_scores <- data.frame(
  subj = subj, time = time, group = grp_row,
  mPCA1 = S_pca[, 1], mPCA2 = S_pca[, 2],
  mPLS1 = S_mpls[, 1], mPLS2 = S_mpls[, 2])
subj_scores <- data.frame(
  subj = meta$subj, age_death = meta$age_death, group = grp_subj,
  tPCA1 = tp$variates[, 1], tPCA2 = tp$variates[, 2],
  tPLS1 = ft$x_projected[, 1], tPLS2 = ft$x_projected[, 2])
fwrite(row_scores, file.path(dat_dir, "fish_scores_row.csv"))
fwrite(subj_scores, file.path(dat_dir, "fish_scores_subject.csv"))
cat("saved:", nrow(row_scores), "rows /", nrow(subj_scores), "subjects\n")
cat("spot check mPCA cor with time:", cor(row_scores$mPCA1, row_scores$time), "\n")
