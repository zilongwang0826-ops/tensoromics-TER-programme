# 22_vaccine_deep.R
# WP2E/F: deep analyses after the chip-batch finding
#   - pooled 5-component check on chip-corrected tPCA
#   - v1 top-8 genes: age effect before vs after chip adjustment (mixed models)
#   - BTM module scores: age effects with chip adjustment (formal tests, FDR)
#   - baseline vs dynamics decomposition (mixed model on tPLS-DA comp1)
#   - boundary: correlations with IgG / Tfh + power analysis
suppressPackageStartupMessages({
  library(data.table); library(mixOmics); library(tensorOmics)
  library(lme4); library(lmerTest)
})
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_tensor_helpers.R"))

base_dir <- Sys.getenv("TER_ROOT", getwd())
prep <- readRDS(file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
corr <- readRDS(file.path(base_dir, "data/GSE79396/chip_correction.rds"))
dat <- readRDS(file.path(base_dir, "data/GSE79396/vaccine_data.rds"))
samples <- fread(file.path(base_dir, "data/GSE79396/samples_all.csv"))
tab_dir <- file.path(base_dir, "results/tables")
set.seed(2028)
subs <- prep$subs; days <- prep$days
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
grp_row <- factor(samples$age_group, levels = c("young", "elderly"))
chip_row <- factor(samples$chip); chip_subj <- factor(samples$chip[match(subs, samples$subject_id)])
sex_subj <- factor(samples$sex[match(subs, samples$subject_id)])
sex_row <- factor(samples$sex)

build_tensor_scaled <- function(Xrow) {
  Z <- scale(Xrow)
  Xt <- array(NA_real_, dim = c(length(subs), ncol(Z), length(days)),
              dimnames = list(subs, colnames(Z), paste0("D", days)))
  for (s in subs) for (d in days) Xt[s, , paste0("D", d)] <- Z[samples$subject_id == s & samples$day == d, ]
  Xt
}
Xt_chip <- build_tensor_scaled(corr$combat_chip_rows)

## ---- 0. pooled 5-component check after chip correction ----------------------
tp5 <- tpca(Xt_chip, ncomp = 5)$variates
man5 <- pillai(tp5, grp_subj)
p5 <- (1 + sum(sapply(1:499, function(b) pillai(tp5, sample(grp_subj))["trace"]) >= man5["trace"])) / 500
writeLines(sprintf("chip-corrected tPCA, 5 comps: Pillai=%.4f, p_param=%.4g, p_perm=%.3f",
                   man5["trace"], man5["p"], p5),
           file.path(tab_dir, "vacc_chip_tpca5.txt"))

## ---- 1. v1 top genes: before vs after chip adjustment ------------------------
raw2000 <- t(dat$expr[prep$sel, ])     # 288 x 2000
v1_genes <- c("201531_PM_at", "216841_PM_s_at", "202638_PM_s_at", "218810_PM_at",
              "202924_PM_s_at", "241889_PM_at", "37028_PM_at", "205067_PM_at", "220091_PM_at")
ann <- fread(file.path(base_dir, "data/GPL13158/annot_clean.csv"))
gene_tab <- data.frame()
for (probe in v1_genes) {
  df <- data.frame(y = raw2000[, probe], age = grp_row, chip = chip_row, sex = sex_row,
                   subject = factor(samples$subject_id))
  m0 <- try(lmer(y ~ age + (1 | subject), data = df), silent = TRUE)
  m1 <- try(lmer(y ~ age + chip + sex + (1 | subject), data = df), silent = TRUE)
  p0 <- if (inherits(m0, "try-error")) NA else summary(m0)$coefficients["ageelderly", "Pr(>|t|)"]
  p1 <- if (inherits(m1, "try-error")) NA else summary(m1)$coefficients["ageelderly", "Pr(>|t|)"]
  b0 <- if (inherits(m0, "try-error")) NA else summary(m0)$coefficients["ageelderly", "Estimate"]
  b1 <- if (inherits(m1, "try-error")) NA else summary(m1)$coefficients["ageelderly", "Estimate"]
  gene_tab <- rbind(gene_tab, data.frame(probe = probe, symbol = ann$symbol[match(probe, ann$probe)],
                                         beta_age_only = b0, p_age_only = p0,
                                         beta_age_chip_sex = b1, p_age_chip_sex = p1))
}
fwrite(gene_tab, file.path(tab_dir, "vacc_v1genes_chip_adjust.csv"))

## ---- 2. BTM modules: age effect with chip adjustment -------------------------
btm <- fread(file.path(base_dir, "data/Li2017_suppl/tidy/btm_scores.csv"))
btm_res <- data.frame()
for (i in seq_len(nrow(btm))) {
  row <- as.numeric(btm[i, -1]); names(row) <- colnames(btm)[-1]
  mod_name <- btm$BTM[i]
  ok <- intersect(names(row), subs)
  if (length(ok) < 60) next
  d <- data.frame(score = row[ok], age = grp_subj[subs %in% ok],
                  chip = chip_subj[subs %in% ok], sex = sex_subj[subs %in% ok])
  m <- lm(score ~ age + chip + sex, data = d)
  co <- summary(m)$coefficients
  if (!"ageelderly" %in% rownames(co)) next
  btm_res <- rbind(btm_res, data.frame(module = mod_name, beta = co["ageelderly", 1],
                                       p = co["ageelderly", 4]))
}
btm_res$FDR <- p.adjust(btm_res$p, "BH")
btm_res <- btm_res[order(btm_res$p), ]
fwrite(btm_res, file.path(tab_dir, "vacc_btm_age_chip.csv"))

## ---- 3. baseline vs dynamics (mixed model on tPLS-DA comp1) ------------------
ft <- tplsda(prep$X_tensor, as.numeric(grp_subj == "elderly"), ncomp = 2)
sc <- data.frame(score = ft$x_projected[, 1], subject = factor(samples$subject_id),
                 day = samples$day, age = grp_row, chip = chip_row, sex = sex_row)
mm <- lmer(score ~ age * factor(day) + chip + sex + (1 | subject), data = sc)
vc <- as.data.frame(VarCorr(mm))
var_sub <- vc$vcov[vc$grp == "subject"]; var_res <- vc$vcov[vc$grp == "Residual"]
mix_out <- list(
  interaction_p = anova(mm)["age:factor(day)", "Pr(>F)"],
  age_p = anova(mm)["age", "Pr(>F)"],
  subject_var_share = var_sub / (var_sub + var_res),
  table = as.data.frame(anova(mm)))
saveRDS(mix_out, file.path(tab_dir, "vacc_baseline_dynamics.rds"))
# per-gene (top 8 by |loading|) interaction test
top8 <- order(abs(ft$x_loadings[, 1]), decreasing = TRUE)[1:8]
t8 <- data.frame()
for (j in top8) {
  probe <- prep$sel[j]
  d <- data.frame(y = raw2000[, probe], age = grp_row, day = factor(samples$day),
                  chip = chip_row, sex = sex_row, subject = factor(samples$subject_id))
  m <- try(lmer(y ~ age * day + chip + sex + (1 | subject), data = d), silent = TRUE)
  if (inherits(m, "try-error")) next
  a <- anova(m)
  t8 <- rbind(t8, data.frame(probe = probe, symbol = ann$symbol[match(probe, ann$probe)],
                             loading = ft$x_loadings[j, 1],
                             p_age = a["age", "Pr(>F)"], p_interaction = a["age:day", "Pr(>F)"]))
}
fwrite(t8, file.path(tab_dir, "vacc_top8_mixed.csv"))

## ---- 4. boundary: IgG / Tfh correlations + power ------------------------------
igg <- fread(file.path(base_dir, "data/Li2017_suppl/tidy/igg_raw.csv"))
tfh <- fread(file.path(base_dir, "data/Li2017_suppl/tidy/tfh_raw.csv"))
S_tpca2 <- prep$scores_tensor[, 1:2]; S_tpls <- ft$x_projected
cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  r <- cor(x[ok], y[ok])
  n <- sum(ok)
  bt <- replicate(1000, cor(x[ok][sample(n, n, TRUE)], y[ok][sample(n, n, TRUE)]))
  data.frame(r = r, n = n, p = cor.test(x[ok], y[ok])$p.value,
             ci_lo = quantile(bt, 0.025), ci_hi = quantile(bt, 0.975))
}
bound <- data.frame()
for (resp in c("D30/0_IgG", "D180/0_IgG")) {
  rv <- as.numeric(igg[igg$subj == resp, -1]); names(rv) <- colnames(igg)[-1]
  x <- S_tpls[, 1]; names(x) <- subs
  ct <- cor_test(x[names(rv)[names(rv) %in% subs]], rv[names(rv) %in% subs])
  bound <- rbind(bound, cbind(space = "tPLS1", response = resp, ct))
}
for (resp in c("D7/0_TFH.CXCR3pos", "D7/0_TFH.CXCR3neg")) {
  rv <- as.numeric(tfh[tfh$Subj == resp, -1]); names(rv) <- colnames(tfh)[-1]
  ct <- cor_test(S_tpls[subs %in% names(rv), 1], rv[subs[subs %in% names(rv)]])
  bound <- rbind(bound, cbind(space = "tPLS1", response = resp, ct))
}
# MDE for n = 72 (80% power, two-sided alpha 0.05)
n <- 72
zcrit <- qnorm(0.975) + qnorm(0.8)
mde_r <- tanh(zcrit / sqrt(n - 3))
writeLines(sprintf("MDE |r| for n=72 at 80%% power: %.3f", mde_r),
           file.path(tab_dir, "vacc_mde.txt"))
fwrite(bound, file.path(tab_dir, "vacc_boundary.csv"))

cat("=== v1 top genes before/after chip adjustment ===\n"); print(gene_tab)
cat("\n=== BTM top modules (age, chip-adjusted) ===\n"); print(head(btm_res, 10))
cat("\n=== baseline vs dynamics (tPLS1) ===\n"); print(mix_out$table); cat("subject var share:", mix_out$subject_var_share, "\n")
cat("\n=== boundary ===\n"); print(bound); cat("MDE:", mde_r, "\n")
cat("\nDONE\n")
