# lib_tensor_helpers.R : shared helpers for faithful tensorOmics projections and tests
suppressPackageStartupMessages({
  library(tensorOmics)
  library(MASS)
})

pillai <- function(scores, g) {
  m <- manova(as.matrix(scores) ~ g)
  s <- summary(m, test = "Pillai")$stats
  c(trace = unname(s[1, "Pillai"]), p = unname(s[1, "Pr(>F)"]))
}

make_dct <- function(t) dctii_m_transforms(t)

# faithful fit + projection of a new subject for tpls/tplsda (verified exact vs package)
tplsda_project <- function(Xtr, ytr, Xnew, dct, ncomp = 2) {
  fit <- tplsda(Xtr, ytr, ncomp = ncomp)
  Mt <- apply(Xtr, c(2, 3), mean)
  xh <- dct$m(sweep(Xtr, c(2, 3), Mt, "-"))
  a <- fit$x_loadings; f <- fit$faces
  b <- matrix(NA_real_, dim(Xtr)[1], ncomp); regs <- vector("list", ncomp)
  xst <- xh
  for (i in seq_len(ncomp)) {
    bi <- xst[, , f[i]] %*% a[, i]
    b[, i] <- bi
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
  list(b = b, bnew = bnew, fit = fit)
}

# subject-grouped LOO predictions for tPLS-DA (LDA on training component scores)
tplsda_loo_predict <- function(X_tensor, grp, dct, ncomp = 2) {
  n <- dim(X_tensor)[1]
  pred <- character(n)
  for (s in seq_len(n)) {
    tr <- setdiff(seq_len(n), s)
    res <- tplsda_project(X_tensor[tr, , , drop = FALSE],
                          as.numeric(grp[tr] == levels(grp)[2]),
                          X_tensor[s, , , drop = FALSE], dct, ncomp)
    ld <- MASS::lda(res$b, grouping = grp[tr])
    pred[s] <- as.character(predict(ld, newdata = matrix(res$bnew, 1, ncomp))$class)
  }
  pred
}

balanced_acc <- function(pred, grp) {
  mean(sapply(levels(grp), function(cl) mean(pred[grp == cl] == cl)))
}
