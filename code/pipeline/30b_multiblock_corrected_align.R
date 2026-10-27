# 30b_multiblock_corrected_align.R : alignment of the chip-corrected canonical axis
suppressPackageStartupMessages({library(data.table); library(tensorOmics)})
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_tensor_helpers.R"))
base_dir <- Sys.getenv("TER_ROOT", getwd())
prep <- readRDS(file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
corr <- readRDS(file.path(base_dir, "data/GSE79396/chip_correction.rds"))
mt <- readRDS(file.path(base_dir, "data/GSE79396/metab_tensors.rds"))
samples <- fread(file.path(base_dir, "data/GSE79396/samples_all.csv"))
subs <- prep$subs; days <- prep$days
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
cohort_subj <- factor(samples$cohort[match(subs, samples$subject_id)])
chip_subj <- factor(samples$chip[match(subs, samples$subject_id)])

X1 <- prep$X_tensor[, 1:500, ]
metab_pos <- mt$pos[, 1:250, ]; metab_neg <- mt$neg[, 1:250, ]
metab_all <- array(NA_real_, dim = c(length(subs), 500, length(days)),
                   dimnames = list(subs, c(dimnames(metab_pos)[[2]], dimnames(metab_neg)[[2]]), paste0("D", days)))
metab_all[, 1:250, ] <- metab_pos; metab_all[, 251:500, ] <- metab_neg

Zc <- scale(corr$combat_chip_rows)[, prep$sel[1:500]]
X1c <- array(NA_real_, dim = dim(X1), dimnames = dimnames(X1))
for (s in subs) for (d in days) X1c[s, , paste0("D", d)] <- Zc[samples$subject_id == s & samples$day == d, ]
can_c <- tpls(X1c, metab_all, ncomp = 2, mode = "canonical")

rows <- data.frame()
for (blk in c("chip-corrected transcriptome")) {
  S <- can_c$x_projected
  for (grp_nm in c("age", "cohort", "chip")) {
    g <- switch(grp_nm, age = grp_subj, cohort = cohort_subj, chip = chip_subj)
    mv <- pillai(S, g)
    rows <- rbind(rows, data.frame(block = blk, variable = grp_nm, pillai = mv["trace"], p = mv["p"]))
  }
}
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_p2a_helpers.R"))
align_path <- file.path(base_dir, "results/tables/vacc_multiblock_align.csv")
old <- fread(align_path)
new_tab <- upsert_rows(old, rows, keys = c("block", "variable"))
atomic_write_csv(new_tab, align_path)
print(rows)
cat("DONE\n")
