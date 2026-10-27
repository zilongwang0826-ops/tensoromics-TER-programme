# 22b_baseline_dynamics.R : proper baseline-vs-dynamics check with day-resolved axis score
suppressPackageStartupMessages({library(data.table); library(tensorOmics); library(lme4); library(lmerTest)})
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_tensor_helpers.R"))
base_dir <- Sys.getenv("TER_ROOT", getwd())
prep <- readRDS(file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
dat <- readRDS(file.path(base_dir, "data/GSE79396/vaccine_data.rds"))
samples <- fread(file.path(base_dir, "data/GSE79396/samples_all.csv"))
tab_dir <- file.path(base_dir, "results/tables")
subs <- prep$subs
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
grp_row <- factor(samples$age_group, levels = c("young", "elderly"))
chip_row <- factor(samples$chip); sex_row <- factor(samples$sex)

# tPLS-DA loadings -> top 20 features -> day-resolved axis score = mean z-expression
ft <- tplsda(prep$X_tensor, as.numeric(grp_subj == "elderly"), ncomp = 2)
top20 <- prep$sel[order(abs(ft$x_loadings[, 2]), decreasing = TRUE)[1:20]]
raw <- t(dat$expr[top20, ])                 # 288 x 20
Z <- scale(raw, center = TRUE, scale = TRUE) # z per feature across 288
score_day <- rowMeans(Z)

d <- data.frame(score = score_day, age = grp_row, day = factor(samples$day),
                chip = chip_row, sex = sex_row, subject = factor(samples$subject_id))
m <- lmer(score ~ age * day + chip + sex + (1 | subject), data = d)
a <- anova(m)
vc <- as.data.frame(VarCorr(m))
share <- vc$vcov[vc$grp == "subject"] / (vc$vcov[vc$grp == "subject"] + vc$vcov[vc$grp == "Residual"])
out <- data.frame(effect = rownames(a), p = a$`Pr(>F)`)
fwrite(out, file.path(tab_dir, "vacc_axis_dynamics_anova.csv"))
writeLines(sprintf("subject-level variance share of the day-resolved axis score: %.3f", share),
           file.path(tab_dir, "vacc_axis_dynamics_varshare.txt"))

# plain (uncorrected) component-level age LMs for the report
p_c1 <- coef(summary(lm(ft$x_projected[, 1] ~ grp_subj)))["grp_subjelderly", 4]
p_c2 <- coef(summary(lm(ft$x_projected[, 2] ~ grp_subj)))["grp_subjelderly", 4]
writeLines(sprintf("tPLS-DA comp1/comp2 age LM p (uncorrected): %.3g / %.3g", p_c1, p_c2),
           file.path(tab_dir, "vacc_tpls_comp_age_p.txt"))
print(out); cat("variance share:", share, "\n"); cat("comp p:", p_c1, p_c2, "\n")
cat("DONE\n")
