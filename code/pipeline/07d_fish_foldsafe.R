#!/usr/bin/env Rscript
# 07d_fish_foldsafe.R — fold-safe re-derivation of the fish grouped LOSO metrics.
# For each of 113 LOSO folds, the top-1000 HVGs are selected on training fish only.
# Covers: matrix PLS-DA grouped LOSO (07 target), tensor tPLS-DA grouped LOSO (07c target).
# Also recomputes the committed (leaky) versions in-process as a parity check.
suppressPackageStartupMessages({ library(data.table); library(mixOmics); library(tensorOmics); library(MASS) })
PROJ <- Sys.getenv("TER_ROOT", getwd())
source(file.path(PROJ, "scripts/lib_tensor_helpers.R"))
dat_dir <- file.path(PROJ, "data/GSE150318")
tab_dir <- file.path(PROJ, "results/tables")

full <- fread(file.path(dat_dir, "fish_full_log2.csv"))
x1k  <- fread(file.path(dat_dir, "fish_X1000_long.csv"))
meta <- fread(file.path(dat_dir, "fish_full_meta.csv"))
stopifnot(all(full$subj == x1k$subj), all(full$time == x1k$time))

genes_full <- setdiff(names(full), c("subj", "time"))
Xf <- as.matrix(full[, ..genes_full])
subj <- full$subj; time <- full$time
n <- nrow(meta)
grp_subj <- factor(meta$group, levels = c(0, 1), labels = c("short_lived", "long_lived"))
grp_row  <- factor(meta$group[match(subj, meta$subj)], levels = c(0, 1),
                   labels = c("short_lived", "long_lived"))
genes1k <- setdiff(names(x1k), c("subj", "time"))
Xm1k <- as.matrix(x1k[, ..genes1k])

build_tensor2 <- function(Mrows, subj_v, time_v, subs_v) {
  Xt <- array(NA_real_, dim = c(length(subs_v), ncol(Mrows), 2),
              dimnames = list(as.character(subs_v), colnames(Mrows), c("10w", "20w")))
  key <- paste(subj_v, time_v, sep = "|")
  for (i in seq_along(subs_v)) for (k in 1:2) {
    j <- match(paste(subs_v[i], c(10, 20)[k], sep = "|"), key)
    if (is.na(j)) stop("build_tensor2: missing row")
    Xt[i, , k] <- Mrows[j, ]
  }
  Xt
}
dct <- make_dct(2)
maj2 <- function(v) names(sort(table(v), decreasing = TRUE))[1]
ba <- function(pred, grp) mean(sapply(levels(grp), function(cl) mean(pred[grp == cl] == cl)))

## leaky parity run (global top-1000, as committed)
Xts <- build_tensor2(Xm1k, subj, time, meta$subj)
pred_lk_m <- character(n); pred_lk_t <- character(n)
for (s in seq_len(n)) {
  hold <- which(subj == meta$subj[s])
  f <- mixOmics::plsda(Xm1k[-hold, , drop = FALSE], grp_row[-hold], ncomp = 2)
  pr <- predict(f, newdata = Xm1k[hold, , drop = FALSE])$class$max.dist
  pred_lk_m[s] <- maj2(as.character(pr[, ncol(pr)]))
  prt <- tplsda_project(Xts[-s, , , drop = FALSE], as.numeric(grp_subj[-s] == "long_lived"),
                        Xts[s, , , drop = FALSE], dct, ncomp = 2)
  ld <- MASS::lda(prt$b, grouping = grp_subj[-s])
  pred_lk_t[s] <- as.character(predict(ld, newdata = matrix(prt$bnew, 1, 2))$class)
}

## fold-safe run (per-fold selection on training fish only)
pred_fs_m <- character(n); pred_fs_t <- character(n); overlaps <- integer(n)
for (s in seq_len(n)) {
  hold <- which(subj == meta$subj[s])
  tr <- setdiff(seq_len(nrow(Xf)), hold)
  v <- apply(Xf[tr, , drop = FALSE], 2, var)
  sel <- order(v, decreasing = TRUE)[seq_len(1000)]
  overlaps[s] <- length(intersect(colnames(Xf)[sel], genes1k))
  f <- mixOmics::plsda(Xf[tr, sel, drop = FALSE], grp_row[tr], ncomp = 2)
  pr <- predict(f, newdata = Xf[hold, sel, drop = FALSE])$class$max.dist
  pred_fs_m[s] <- maj2(as.character(pr[, ncol(pr)]))
  Xt_fs <- build_tensor2(Xf[tr, sel, drop = FALSE], subj[tr], time[tr], as.character(meta$subj[-s]))
  Xt_ne <- build_tensor2(Xf[hold, sel, drop = FALSE], subj[hold], time[hold], as.character(meta$subj[s]))
  prt <- tplsda_project(Xt_fs, as.numeric(grp_subj[-s] == "long_lived"), Xt_ne, dct, ncomp = 2)
  ld <- MASS::lda(prt$b, grouping = grp_subj[-s])
  pred_fs_t[s] <- as.character(predict(ld, newdata = matrix(prt$bnew, 1, 2))$class)
}

out <- data.frame(
  metric = c("matrix PLS-DA grouped LOSO", "tPLS-DA grouped LOSO"),
  leaky = c(ba(pred_lk_m, grp_subj), ba(pred_lk_t, grp_subj)),
  foldsafe = c(ba(pred_fs_m, grp_subj), ba(pred_fs_t, grp_subj)))
out$delta <- out$foldsafe - out$leaky
fwrite(out, file.path(tab_dir, "fish_foldsafe.csv"))
fwrite(data.frame(subj = meta$subj, grp = as.character(grp_subj),
                  leaky_matrix = pred_lk_m, fs_matrix = pred_fs_m,
                  leaky_tensor = pred_lk_t, fs_tensor = pred_fs_t,
                  sel_overlap = overlaps),
       file.path(tab_dir, "fish_foldsafe_predictions.csv"))
cat("leaky parity (expect 0.620056 / 0.681576):\n"); print(out)
cat(sprintf("mean overlap of per-fold selection with global top-1000: %.1f / 1000\n", mean(overlaps)))
cat("DONE\n")
