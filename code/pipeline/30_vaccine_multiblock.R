# 30_vaccine_multiblock.R
# WP3: multi-block integration (transcriptome x metabolome), feature-reduced for tractability.
#   - top 500 transcriptome features; top 250 pos + 250 neg metabolites
#   - tPLS canonical + subject-pairing permutation
#   - block.tPLS-DA (3 blocks) on age; per-block age tests + permutation
#   - cross-layer correlation data (Fig 7C analog)
#   - chip-corrected sensitivity
suppressPackageStartupMessages({
  library(data.table); library(tensorOmics); library(MASS)
})
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_tensor_helpers.R"))

base_dir <- Sys.getenv("TER_ROOT", getwd())
prep <- readRDS(file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
corr <- readRDS(file.path(base_dir, "data/GSE79396/chip_correction.rds"))
mt <- readRDS(file.path(base_dir, "data/GSE79396/metab_tensors.rds"))
samples <- fread(file.path(base_dir, "data/GSE79396/samples_all.csv"))
tab_dir <- file.path(base_dir, "results/tables")
set.seed(2029)
subs <- prep$subs; days <- prep$days
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
cohort_subj <- factor(samples$cohort[match(subs, samples$subject_id)])
chip_subj <- factor(samples$chip[match(subs, samples$subject_id)])
n_subj <- length(subs)

# ---- feature reduction -------------------------------------------------------
X1 <- prep$X_tensor[, 1:500, ]                    # top-500 by variance (order in sel)
metab_pos <- mt$pos[, 1:250, ]
metab_neg <- mt$neg[, 1:250, ]
metab_all <- array(NA_real_, dim = c(n_subj, 500, length(days)),
                   dimnames = list(subs, c(dimnames(metab_pos)[[2]], dimnames(metab_neg)[[2]]), paste0("D", days)))
metab_all[, 1:250, ] <- metab_pos
metab_all[, 251:500, ] <- metab_neg
cat("dims: expr", dim(X1), " metab", dim(metab_all), "\n")

## ---------- metabolome audit --------------------------------------------------
tp_m <- tpca(metab_all, ncomp = 2)$variates
man_m <- rbind(age = pillai(tp_m, grp_subj), cohort = pillai(tp_m, cohort_subj), chip = pillai(tp_m, chip_subj))

## ---------- tPLS canonical ----------------------------------------------------
t0 <- Sys.time()
can <- tpls(X1, metab_all, ncomp = 2, mode = "canonical")
cat("tpls canonical:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
sx <- can$x_projected; sy <- can$y_projected
r_comp <- sapply(1:2, function(k) cor(sx[, k], sy[, k]))

nperm <- 99
null_r <- numeric(nperm)
for (b in 1:nperm) {
  perm <- sample(n_subj)
  cb <- tpls(X1, metab_all[perm, , , drop = FALSE], ncomp = 2, mode = "canonical")
  null_r[b] <- cor(cb$x_projected[, 1], cb$y_projected[, 1])
}
p_pair <- (1 + sum(abs(null_r) >= abs(r_comp[1]))) / (1 + nperm)
cat("comp1 r =", round(r_comp[1], 3), " pairing perm p =", p_pair, "\n")

align <- data.frame()
for (blk in c("transcriptome", "metabolome")) {
  S <- if (blk == "transcriptome") sx else sy
  for (grp_nm in c("age", "cohort", "chip")) {
    g <- switch(grp_nm, age = grp_subj, cohort = cohort_subj, chip = chip_subj)
    mv <- pillai(S, g)
    align <- rbind(align, data.frame(block = blk, variable = grp_nm,
                                     pillai = mv["trace"], p = mv["p"]))
  }
}
# chip-corrected sensitivity
Zc <- scale(corr$combat_chip_rows)[, prep$sel[1:500]]
X1c <- array(NA_real_, dim = dim(X1), dimnames = dimnames(X1))
for (s in subs) for (d in days) X1c[s, , paste0("D", d)] <- Zc[samples$subject_id == s & samples$day == d, ]
can_c <- tpls(X1c, metab_all, ncomp = 2, mode = "canonical")
r_comp_c <- cor(can_c$x_projected[, 1], can_c$y_projected[, 1])

## ---------- block.tPLS-DA ------------------------------------------------------
ybin <- as.numeric(grp_subj == "elderly")
t0 <- Sys.time()
bt <- block_tplsda(x = list(expr = X1, pos = metab_pos, neg = metab_neg), y = ybin, ncomp = 2)
cat("block_tplsda:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
bp <- bt$projected
block_man <- data.frame()
for (blk in c("expr", "pos", "neg")) {
  mv <- pillai(bp[[blk]], grp_subj)
  block_man <- rbind(block_man, data.frame(block = blk, pillai = mv["trace"], p = mv["p"]))
}
# permutation for the block model (age labels, refit)
nperm2 <- 99
null_bt <- numeric(nperm2)
for (b in 1:nperm2) {
  g <- sample(grp_subj)
  fb <- try(block_tplsda(x = list(expr = X1, pos = metab_pos, neg = metab_neg),
                         y = as.numeric(g == "elderly"), ncomp = 2), silent = TRUE)
  null_bt[b] <- if (inherits(fb, "try-error")) NA else pillai(fb$projected$expr, g)["trace"]
}
p_bt <- (1 + sum(null_bt >= block_man$pillai[block_man$block == "expr"], na.rm = TRUE)) / (1 + sum(!is.na(null_bt)))

## ---------- cross-layer correlation data ---------------------------------------
lx <- can$x_loadings[, 1]; ly <- can$y_loadings[, 1]
top_expr <- order(abs(lx), decreasing = TRUE)[1:30]
top_metab <- order(abs(ly), decreasing = TRUE)[1:30]
rows <- list()
for (d in 1:length(days)) {
  E <- X1[, top_expr, d]; M <- metab_all[, top_metab, d]
  R <- cor(E, M)   # 30 x 30
  rows[[d]] <- data.table(day = days[d],
                          expr_feat = rep(dimnames(X1)[[2]][top_expr], times = 30),
                          metab_feat = rep(dimnames(metab_all)[[2]][top_metab], each = 30),
                          r = as.vector(R))
}
cross <- rbindlist(rows)

## ---------- save ---------------------------------------------------------------
fwrite(data.table(comparison = c("canonical comp1", "canonical comp2"),
                  block_score_r = r_comp, pairing_perm_p = c(p_pair, NA),
                  chip_corrected_comp1_r = c(r_comp_c, NA)),
       file.path(tab_dir, "vacc_multiblock_canonical.csv"))
fwrite(align, file.path(tab_dir, "vacc_multiblock_align.csv"))
fwrite(block_man, file.path(tab_dir, "vacc_multiblock_blocktplsda.csv"))
writeLines(sprintf("block.tPLS-DA expr-block perm p = %.3f (n_perm=%d)", p_bt, sum(!is.na(null_bt))),
           file.path(tab_dir, "vacc_multiblock_blocktplsda_perm.txt"))
fwrite(cross, file.path(tab_dir, "vacc_multiblock_crosslayer.csv"))
fwrite(data.table(subject = subs, age_group = grp_subj, cohort = cohort_subj, chip = chip_subj,
                  can_tx1 = sx[, 1], can_tx2 = sx[, 2], can_mb1 = sy[, 1], can_mb2 = sy[, 2],
                  bt_expr1 = bp$expr[, 1], bt_pos1 = bp$pos[, 1], bt_neg1 = bp$neg[, 1]),
       file.path(tab_dir, "vacc_multiblock_scores.csv"))

cat("=== metabolome audit ===\n"); print(man_m)
fwrite(data.table(variable = rownames(man_m), as.data.frame(man_m)), file.path(tab_dir, "vacc_metab_audit.csv"))
cat("=== alignment ===\n"); print(align)
cat("=== block.tPLS-DA ===\n"); print(block_man); cat("expr block perm p:", p_bt, "\n")
cat("chip-corrected canonical comp1 r:", r_comp_c, "\n")
cat("DONE\n")
