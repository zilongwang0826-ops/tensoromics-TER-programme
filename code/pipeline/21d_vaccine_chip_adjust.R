# 21d_vaccine_chip_adjust.R
# Complementary check: adjust for chip as a covariate in the model (non-destructive)
# instead of correcting the data. Is there age signal beyond chip?
suppressPackageStartupMessages({library(data.table); library(tensorOmics)})
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_tensor_helpers.R"))
base_dir <- Sys.getenv("TER_ROOT", getwd())
prep <- readRDS(file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
samples <- fread(file.path(base_dir, "data/GSE79396/samples_all.csv"))
tab_dir <- file.path(base_dir, "results/tables")
subs <- prep$subs
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
chip_subj <- factor(samples$chip[match(subs, samples$subject_id)])
sex_subj <- factor(samples$sex[match(subs, samples$subject_id)])

S_tpca <- prep$scores_tensor[, 1:2]
ft <- tplsda(prep$X_tensor, as.numeric(grp_subj == "elderly"), ncomp = 2)
S_tpls <- ft$x_projected

out <- data.frame()
for (nm in c("tPCA", "tPLS-DA")) {
  S <- if (nm == "tPCA") S_tpca else S_tpls
  full <- lm(cbind(c1 = S[, 1], c2 = S[, 2]) ~ grp_subj + chip_subj + sex_subj)
  noage <- lm(cbind(c1 = S[, 1], c2 = S[, 2]) ~ chip_subj + sex_subj)
  cmp <- anova(noage, full, test = "Pillai")
  age_p_mv <- cmp[2, "Pr(>F)"]
  out <- rbind(out, data.frame(space = nm, model = "age + chip + sex; MANOVA age term",
                               age_stat = age_p_mv, chip_stat = NA))
  for (k in 1:2) {
    m <- lm(S[, k] ~ grp_subj + chip_subj + sex_subj)
    out <- rbind(out, data.frame(space = paste0(nm, " comp", k), model = "LM age coefficient",
                                 age_stat = coef(summary(m))["grp_subjelderly", 4], chip_stat = NA))
  }
}
print(out)
fwrite(out, file.path(tab_dir, "vacc_chip_adjustment.csv"))
cat("\nNote: chip x age MANOVA on 2 comps; age_p_mv comes from model comparison.\n")
cat("DONE\n")
