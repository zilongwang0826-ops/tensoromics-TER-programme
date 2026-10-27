suppressPackageStartupMessages({library(PLSDAbatch); library(tensorOmics)})
cat("PLSDA_batch args:\n"); print(args(PLSDA_batch))
set.seed(1)
x1 <- array(rnorm(20 * 30 * 3), c(20, 30, 3))
x2 <- array(rnorm(20 * 10 * 3), c(20, 10, 3))
y <- rep(0:1, each = 10)
bt <- try(block_tplsda(x = list(A = x1, B = x2), y = y, ncomp = 2), silent = TRUE)
if (inherits(bt, "try-error")) {
  cat("block_tplsda error:", conditionMessage(attr(bt, "condition")), "\n")
} else {
  cat("block_tplsda names:", paste(names(bt), collapse = ", "), "\n")
  str(bt$projected, max.level = 1)  # correct field name
}
# also check block_tpls list input
Y3 <- array(rep(cbind(y, 1 - y), times = 3), dim = c(length(y), 2, 3))
bl <- try(block_tpls(a = list(A = x1, B = x2, Y = Y3), y_ind = 3, ncomp = 2), silent = TRUE)
if (inherits(bl, "try-error")) cat("block_tpls error:", conditionMessage(attr(bl, "condition")), "\n") else {
  cat("block_tpls names:", paste(names(bl), collapse = ", "), "\n")
}
