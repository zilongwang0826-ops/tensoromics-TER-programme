# 08_fish_tpls_regression.R : cover the tPLS regression mode (continuous age at death)
suppressPackageStartupMessages({library(data.table); library(tensorOmics)})
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_tensor_helpers.R"))
base_dir <- Sys.getenv("TER_ROOT", getwd())
dat_dir <- file.path(base_dir, "data/GSE150318")
tab_dir <- file.path(base_dir, "results/tables")

X <- fread(file.path(dat_dir, "fish_X1000_long.csv"))
meta <- fread(file.path(dat_dir, "fish_meta.csv"))
genes <- setdiff(names(X), c("subj", "time"))
Xmat <- as.matrix(X[, ..genes])
rownames(Xmat) <- paste0("s", X$subj, "_t", X$time)
n <- nrow(meta)
X_tensor <- array(NA_real_, dim = c(n, length(genes), 2),
                  dimnames = list(paste0("s", meta$subj), genes, c("10w", "20w")))
for (i in 1:n) for (k in 1:2) X_tensor[i, , k] <- Xmat[paste0("s", meta$subj[i], "_t", c(10, 20)[k]), ]
dct <- make_dct(2)

# faithful projection recipe generalised to tpls regression fit
tpls_reg_project <- function(Xtr, ytr, Xnew, ncomp = 2) {
  Ytr <- array(rep(ytr, 2), dim = c(length(ytr), 2, 2))
  fit <- tpls(Xtr, Ytr, ncomp = ncomp, mode = "regression")
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

# full fit
Y <- array(rep(meta$age_death, 2), dim = c(n, 2, 2))  # duplicated response column (implementation requires q>=2)
fit <- tpls(X_tensor, Y, ncomp = 2, mode = "regression")
S <- fit$x_projected
r1 <- cor(S[, 1], meta$age_death)
r2_model <- summary(lm(meta$age_death ~ S))$r.squared
cat("tPLS regression: r(comp1, age) =", round(r1, 3), "| in-sample R2 =", round(r2_model, 3), "\n")

set.seed(31)
nperm <- 199
p_perm <- (1 + sum(sapply(1:nperm, function(b) {
  g <- sample(meta$age_death)
  f <- tpls(X_tensor, array(rep(g, 2), dim = c(n, 2, 2)), ncomp = 2, mode = "regression")
  summary(lm(g ~ f$x_projected))$r.squared
}) >= r2_model)) / (1 + nperm)
cat("subject-level permutation p =", p_perm, "\n")

# LOSO CV
pred <- numeric(n)
for (s in 1:n) {
  tr <- setdiff(1:n, s)
  res <- tpls_reg_project(X_tensor[tr, , , drop = FALSE], meta$age_death[tr],
                          X_tensor[s, , , drop = FALSE], ncomp = 2)
  md <- lm(meta$age_death[tr] ~ res$b)
  pred[s] <- predict(md, newdata = data.frame(res$bnew))
}
cv_r2 <- 1 - sum((meta$age_death - pred)^2) / sum((meta$age_death - mean(meta$age_death))^2)
cat("LOSO CV R2 =", round(cv_r2, 3), "\n")

out <- data.frame(model = "tPLS regression (fish, age at death)",
                  r_comp1 = r1, in_sample_R2 = r2_model, perm_p = p_perm, loso_cv_R2 = cv_r2)
fwrite(out, file.path(tab_dir, "fish_tpls_regression.csv"))
fwrite(data.frame(subj = meta$subj, age_death = meta$age_death, comp1 = S[, 1], comp2 = S[, 2],
                  loso_pred = pred), file.path(dat_dir, "fish_tpls_reg_scores.csv"))
cat("DONE\n")
