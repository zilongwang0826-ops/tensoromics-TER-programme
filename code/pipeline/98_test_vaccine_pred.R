suppressPackageStartupMessages({library(tensorOmics)})
base_dir <- Sys.getenv("TER_ROOT", getwd())
prep <- readRDS(file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
samples <- read.csv(file.path(base_dir, "data/GSE79396/samples_all.csv"))
subs <- prep$subs
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
Xt <- prep$X_tensor

# 1) m_product vs direct m()
set.seed(11)  # fixed seed: the diagnostic input tensor is deterministic
Md <- dctii_m_transforms(4)
X <- array(rnorm(72 * 10 * 4), c(72, 10, 4))
A <- Md$m(X)
B <- m_product(X, m = Md$m, minv = Md$minv)
# corrected equivalence check: for a single input tensor, m_product(X) = minv(m(X))
cat("max|minv(m(X)) - m_product(X)|:", max(abs(Md$minv(A) - B)), "\n")

ft_model <- tplsda(Xt, as.numeric(grp_subj == "elderly"), ncomp = 2)
a1 <- ft_model$x_loadings[, 1]; a2 <- ft_model$x_loadings[, 2]
f1 <- ft_model$faces[1]; f2 <- ft_model$faces[2]

M <- apply(Xt, c(2, 3), mean)
Xc <- sweep(Xt, c(2, 3), M, "-")
xhat <- Md$m(Xc)

b1 <- xhat[, , f1] %*% a1
cat("direct m(): b1 vs proj1 cor =", cor(b1, ft_model$x_projected[, 1]),
    " max|diff| =", max(abs(b1 - ft_model$x_projected[, 1])), "\n")

r1 <- crossprod(xhat[, , f1], b1) / as.numeric(crossprod(b1))
xhat1 <- xhat
xhat1[, , f1] <- xhat[, , f1] - tcrossprod(b1, r1)
b2 <- xhat1[, , f2] %*% a2
cat("direct m(): b2 vs proj2 cor =", cor(b2, ft_model$x_projected[, 2]),
    " max|diff| =", max(abs(b2 - ft_model$x_projected[, 2])), "\n")
