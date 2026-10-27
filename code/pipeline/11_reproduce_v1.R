# 11_reproduce_v1.R
# Reproduce the v1 vaccine tPCA pipeline (results_GSE79396_tPCA) exactly.
# Strategy: try candidate preprocessing pipelines, compare scores/loadings to v1 targets.
suppressPackageStartupMessages({
  library(data.table)
  library(mixOmics)
  library(tensorOmics)
})

base_dir <- Sys.getenv("TER_ROOT", getwd())
v1_dir <- Sys.getenv("TER_V1_DIR", file.path(base_dir, "data/v1_tables"))
if (!dir.exists(v1_dir)) {
  stop("v1 result tables not found: ", v1_dir,
       "\nSet TER_V1_DIR to the directory containing matrix_pca_scores.csv etc.")
}
dat <- readRDS(file.path(base_dir, "data/GSE79396/vaccine_data.rds"))
expr <- dat$expr          # probes x samples (log2 RMA)
samples <- dat$samples
stopifnot(ncol(expr) == nrow(samples))

v1_m <- fread(file.path(v1_dir, "matrix_pca_scores.csv"))
v1_t <- fread(file.path(v1_dir, "tensor_pca_scores.csv"))
v1_l1 <- fread(file.path(v1_dir, "top_loadings_pc1.csv"))
v1_l2 <- fread(file.path(v1_dir, "top_loadings_pc2.csv"))
cat("v1 matrix scores:", dim(v1_m), "| v1 tensor scores:", dim(v1_t), "\n")
cat("v1 top loadings pc1:", nrow(v1_l1), " pc2:", nrow(v1_l2), "\n")

## ---- candidate feature selection -------------------------------------------
# Candidate A: top 2000 by variance of raw log2 values (across all 288 samples)
v <- apply(expr, 1, var)
selA <- names(sort(v, decreasing = TRUE))[1:2000]
cat("selA probes:", length(selA), "\n")

## ---- candidate A: matrix PCA (center only, no scaling) ----------------------
Xm <- t(expr[selA, ])                      # 288 x 2000
rownames(Xm) <- samples$gsm

pca_mix <- mixOmics::pca(Xm, ncomp = 5, center = TRUE, scale = FALSE)
S_mix <- pca_mix$variates$X

pca_pr <- prcomp(Xm, center = TRUE, scale. = FALSE)
S_pr <- pca_pr$x[, 1:5]

pca_scale <- mixOmics::pca(Xm, ncomp = 5, center = TRUE, scale = TRUE)
S_sc <- pca_scale$variates$X

cmp <- merge(v1_m[, .(gsm = geo_accession, v1_1 = PC1, v1_2 = PC2, v1_3 = PC3, v1_4 = PC4, v1_5 = PC5)],
             data.table(gsm = rownames(Xm), mix1 = S_mix[,1], mix2 = S_mix[,2], mix3 = S_mix[,3], mix4 = S_mix[,4], mix5 = S_mix[,5],
                        pr1 = S_pr[,1], pr2 = S_pr[,2], pr3 = S_pr[,3], pr4 = S_pr[,4], pr5 = S_pr[,5],
                        sc1 = S_sc[,1], sc2 = S_sc[,2], sc3 = S_sc[,3], sc4 = S_sc[,4], sc5 = S_sc[,5]),
             by = "gsm")
stopifnot(nrow(cmp) == 288)

cors <- sapply(1:5, function(k) c(
  mixOmics = cor(cmp[[paste0("v1_", k)]], cmp[[paste0("mix", k)]]),
  prcomp   = cor(cmp[[paste0("v1_", k)]], cmp[[paste0("pr", k)]]),
  scaled   = cor(cmp[[paste0("v1_", k)]], cmp[[paste0("sc", k)]])))
colnames(cors) <- paste0("PC", 1:5)
cat("\n=== matrix PCA score correlations (v1 vs candidates) ===\n")
print(round(cors, 4))

## loading overlap check (top 50 by |loading|)
lam <- pca_mix$loadings$X[, 1:2]
top_mix <- lapply(1:2, function(k) {
  o <- order(abs(lam[, k]), decreasing = TRUE)[1:50]
  rownames(lam)[o]
})
ov <- sapply(1:2, function(k) {
  v1l <- if (k == 1) v1_l1$probe_id else v1_l2$probe_id
  sum(top_mix[[k]] %in% v1l)
})
cat("loading top50 overlap with v1 (mixOmics, center-only): PC1 =", ov[1], " PC2 =", ov[2], "\n")
# symbol check of v1 PC2 top probes
ann <- as.data.table(dat$annot)
cat("\nv1 PC2 top annotated probes:\n")
print(head(ann[probe %in% v1_l2$probe_id][order(match(probe, v1_l2$probe_id))], 10))

## ---- candidate A: tensor tPCA ------------------------------------------------
# build tensor: subjects x probes x days (30/10 ordering not important, join later)
subs <- unique(samples$subject_id)
days <- sort(unique(samples$day))
Z <- t(scale(t(expr[selA, ])))   # probe-wise z-score across all 288 samples
X_t <- array(NA_real_, dim = c(length(subs), length(selA), length(days)),
             dimnames = list(subs, selA, paste0("D", days)))
X_t_raw <- X_t
for (s in subs) for (d in days) {
  g <- samples$gsm[samples$subject_id == s & samples$day == d]
  X_t[s, , paste0("D", d)] <- Z[selA, g]
  X_t_raw[s, , paste0("D", d)] <- expr[selA, g]
}
stopifnot(!anyNA(X_t), !anyNA(X_t_raw))

tp  <- tpca(X_t, ncomp = 5)
tpR <- tpca(X_t_raw, ncomp = 5)
cat("\ntpca returned names:", paste(names(tp), collapse = ", "), "\n")
cat("tpca variates dims:", paste(dim(tp$variates), collapse = " x "), "\n")
S_t  <- tp$variates
S_tR <- tpR$variates
idx <- match(v1_t$subject_id, rownames(S_t))
cors_t  <- sapply(1:5, function(k) cor(v1_t[[paste0("Comp", k)]], S_t[idx, k]))
cors_tR <- sapply(1:5, function(k) cor(v1_t[[paste0("Comp", k)]], S_tR[idx, k]))
cat("=== tensor tPCA score correlations (v1 vs scaled-input) ===\n"); print(round(cors_t, 4))
cat("=== tensor tPCA score correlations (v1 vs raw-input) ===\n"); print(round(cors_tR, 4))

## ---- save reproduction report -------------------------------------------------
rep_lines <- c(
  "GSE79396 v1 reproduction check",
  paste("matrix PCA score correlations (rows = candidate, cols = PC1..PC5):"),
  capture.output(print(round(cors, 4))),
  paste("top50 loading overlap: PC1 =", ov[1], " PC2 =", ov[2]),
  paste("tensor tPCA score correlations (Comp1..Comp5):"),
  capture.output(print(round(cors_t, 4)))
)
writeLines(rep_lines, file.path(base_dir, "results/tables/repro_v1_check.txt"))
cat("\n", paste(rep_lines, collapse = "\n"), "\n")

## save candidate A preprocessed objects for later scripts
saveRDS(list(sel = selA, Z = Z, Xm_scaled = Xm, X_tensor = X_t, X_tensor_raw = X_t_raw,
             days = days, subs = subs,
             scores_matrix = S_mix, scores_tensor = S_t),
        file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
cat("\nDONE\n")
