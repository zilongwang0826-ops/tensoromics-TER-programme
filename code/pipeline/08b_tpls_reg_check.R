# 08b_tpls_reg_check.R : verify projection recipe for tpls regression + LOSO diagnostics
suppressPackageStartupMessages({library(data.table); library(tensorOmics)})
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_tensor_helpers.R"))
base_dir <- Sys.getenv("TER_ROOT", getwd())
dat_dir <- file.path(base_dir, "data/GSE150318")
X <- fread(file.path(dat_dir, "fish_X1000_long.csv"))
meta <- fread(file.path(dat_dir, "fish_meta.csv"))
genes <- setdiff(names(X), c("subj", "time"))
Xmat <- as.matrix(X[, ..genes]); rownames(Xmat) <- paste0("s", X$subj, "_t", X$time)
n <- nrow(meta)
X_tensor <- array(NA_real_, dim = c(n, length(genes), 2), dimnames = list(paste0("s", meta$subj), genes, c("10w", "20w")))
for (i in 1:n) for (k in 1:2) X_tensor[i, , k] <- Xmat[paste0("s", meta$subj[i], "_t", c(10, 20)[k]), ]
dct <- make_dct(2)
Y <- array(rep(meta$age_death, 2), dim = c(n, 2, 2))
fit <- tpls(X_tensor, Y, ncomp = 2, mode = "regression")

project_manual <- function(Xtr, fit, Xnew, ncomp = 2) {
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

chk <- project_manual(X_tensor, fit, X_tensor[1, , , drop = FALSE])
cat("recipe check (full-fit projections):\n")
for (i in 1:2) cat(sprintf("comp%d: cor=%.6f max|diff|=%.2e\n", i,
                           cor(chk$b[, i], fit$x_projected[, i]),
                           max(abs(chk$b[, i] - fit$x_projected[, i]))))

# LOSO with the same recipe
pred <- numeric(n); cor_fold <- numeric(n)
for (s in 1:n) {
  tr <- setdiff(1:n, s)
  Ytr <- array(rep(meta$age_death[tr], 2), dim = c(length(tr), 2, 2))
  ftr <- tpls(X_tensor[tr, , , drop = FALSE], Ytr, ncomp = 2, mode = "regression")
  res <- project_manual(X_tensor[tr, , , drop = FALSE], ftr, X_tensor[s, , , drop = FALSE])
  md <- lm(meta$age_death[tr] ~ res$b)
  pred[s] <- predict(md, newdata = data.frame(res$bnew))
  cor_fold[s] <- cor(res$b[, 1], meta$age_death[tr])
}
R2 <- 1 - sum((meta$age_death - pred)^2) / sum((meta$age_death - mean(meta$age_death))^2)
cat("LOSO R2:", round(R2, 3), "| cor(pred, age):", round(cor(pred, meta$age_death), 3), "\n")
cat("per-fold cor(b1, age):", paste(round(summary(cor_fold), 3), collapse = ", "), "\n")
bnew_all <- numeric(n)
for (s in 1:n) {
  tr <- setdiff(1:n, s)
  Ytr <- array(rep(meta$age_death[tr], 2), dim = c(length(tr), 2, 2))
  ftr <- tpls(X_tensor[tr, , , drop = FALSE], Ytr, ncomp = 2, mode = "regression")
  res <- project_manual(X_tensor[tr, , , drop = FALSE], ftr, X_tensor[s, , , drop = FALSE])
  bnew_all[s] <- res$bnew[1]
}
cat("cor(bnew_comp1 across folds, age):", round(cor(bnew_all, meta$age_death), 3), "\n")
