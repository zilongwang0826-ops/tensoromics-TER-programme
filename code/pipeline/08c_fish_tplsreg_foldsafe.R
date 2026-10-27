#!/usr/bin/env Rscript
# 08c_fish_tplsreg_foldsafe.R — fold-safe re-derivation of the tPLS-regression LOSO (script 08).
# Per fold: top-1000 HVGs selected on training fish only; tensor + manual regression projection;
# LM prediction on the held-out fish. Leaky parity run included for comparison.
suppressPackageStartupMessages({ library(data.table); library(tensorOmics) })
PROJ <- Sys.getenv("TER_ROOT", getwd())
source(file.path(PROJ, "scripts/lib_tensor_helpers.R"))
dat_dir <- file.path(PROJ, "data/GSE150318")
tab_dir <- file.path(PROJ, "results/tables")
full <- fread(file.path(dat_dir, "fish_full_log2.csv"))
meta <- fread(file.path(dat_dir, "fish_full_meta.csv"))
x1k  <- fread(file.path(dat_dir, "fish_X1000_long.csv"))
stopifnot(all(full$subj == x1k$subj), all(full$time == x1k$time))
genes_full <- setdiff(names(full), c("subj", "time"))
genes1k <- setdiff(names(x1k), c("subj", "time"))
Xf <- as.matrix(full[, ..genes_full]); subj <- full$subj; time <- full$time
n <- nrow(meta); ages <- meta$age_death
dct <- make_dct(2)

tpls_reg_project <- function(Xtr, ytr, Xnew, ncomp = 2) {
  Ytr <- array(rep(ytr, 2), dim = c(length(ytr), 2, 2))
  fit <- tpls(Xtr, Ytr, ncomp = ncomp, mode = "regression")
  Mt <- apply(Xtr, c(2, 3), mean)
  xh <- dct$m(sweep(Xtr, c(2, 3), Mt, "-"))
  a <- fit$x_loadings; f <- fit$faces
  b <- matrix(NA_real_, dim(Xtr)[1], ncomp); regs <- vector("list", ncomp)
  xst <- xh
  for (i in seq_len(ncomp)) {
    bi <- xst[, , f[i]] %*% a[, i]; b[, i] <- bi
    ri <- crossprod(xst[, , f[i]], bi) / as.numeric(crossprod(bi))
    regs[[i]] <- as.numeric(ri)
    xst[, , f[i]] <- xst[, , f[i]] - tcrossprod(bi, ri)
  }
  xh_new <- dct$m(sweep(Xnew, c(2, 3), Mt, "-"))
  bnew <- numeric(ncomp)
  for (i in seq_len(ncomp)) {
    bnew[i] <- as.numeric(xh_new[, , f[i]] %*% a[, i])
    xh_new[, , f[i]] <- xh_new[, , f[i]] - tcrossprod(bnew[i], regs[[i]])
  }
  list(b = b, bnew = bnew)
}
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

## leaky parity (global top-1000, as committed by 08)
Xm1k <- as.matrix(x1k[, ..genes1k])
Xts <- build_tensor2(Xm1k, subj, time, meta$subj)
pred_lk <- numeric(n)
for (s in 1:n) {
  tr <- setdiff(1:n, s)
  res <- tpls_reg_project(Xts[tr, , , drop = FALSE], ages[tr], Xts[s, , , drop = FALSE])
  md <- lm(ages[tr] ~ res$b); pred_lk[s] <- predict(md, newdata = data.frame(res$bnew))
}
r2_lk <- 1 - sum((ages - pred_lk)^2) / sum((ages - mean(ages))^2)

## fold-safe
pred_fs <- numeric(n); sel_ov <- integer(n)
for (s in 1:n) {
  hold <- which(subj == meta$subj[s]); tr <- setdiff(seq_len(nrow(Xf)), hold)
  v <- apply(Xf[tr, , drop = FALSE], 2, var)
  sel <- order(v, decreasing = TRUE)[seq_len(1000)]
  sel_ov[s] <- length(intersect(colnames(Xf)[sel], genes1k))
  Xt_fs <- build_tensor2(Xf[tr, sel, drop = FALSE], subj[tr], time[tr], as.character(meta$subj[-s]))
  Xt_ne <- build_tensor2(Xf[hold, sel, drop = FALSE], subj[hold], time[hold], as.character(meta$subj[s]))
  res <- tpls_reg_project(Xt_fs, ages[-s], Xt_ne)
  md <- lm(ages[-s] ~ res$b); pred_fs[s] <- predict(md, newdata = data.frame(res$bnew))
}
r2_fs <- 1 - sum((ages - pred_fs)^2) / sum((ages - mean(ages))^2)

out <- data.frame(metric = "tPLS regression LOSO R2",
                  leaky = r2_lk, foldsafe = r2_fs, delta = r2_fs - r2_lk,
                  mean_sel_overlap = mean(sel_ov))
fwrite(out, file.path(tab_dir, "fish_tplsreg_foldsafe.csv"))
fwrite(data.frame(subj = meta$subj, age = ages, leaky_pred = pred_lk, fs_pred = pred_fs,
                  sel_overlap = sel_ov),
       file.path(tab_dir, "fish_tplsreg_foldsafe_predictions.csv"))
print(out)
cat("DONE\n")

