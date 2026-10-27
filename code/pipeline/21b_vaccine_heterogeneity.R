# 21b_vaccine_heterogeneity.R
# Robustness of the cohort-heterogeneity finding:
#   - within-cohort tPCA with 5 components (before and after ComBat)
#   - age distributions by cohort and age group
suppressPackageStartupMessages({
  library(data.table); library(tensorOmics)
})
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_tensor_helpers.R"))
base_dir <- Sys.getenv("TER_ROOT", getwd())
prep <- readRDS(file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
corr <- readRDS(file.path(base_dir, "data/GSE79396/corrected_versions.rds"))
samples <- fread(file.path(base_dir, "data/GSE79396/samples_all.csv"))
tab_dir <- file.path(base_dir, "results/tables")
subs <- prep$subs
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
cohort_subj <- factor(samples$cohort[match(subs, samples$subject_id)])
set.seed(11)

build_tensor_scaled <- function(Xrow) {
  Z <- scale(Xrow)
  Xt <- array(NA_real_, dim = c(length(subs), ncol(Z), 4),
              dimnames = list(subs, colnames(Z), paste0("D", 1:4)))
  for (s in subs) for (k in 1:4) {
    d <- c(0, 1, 3, 7)[k]
    Xt[s, , paste0("D", k)] <- Z[samples$subject_id == s & samples$day == d, ]
  }
  Xt
}
Xt_before <- prep$X_tensor
Xt_combat <- build_tensor_scaled(corr$combat_rows)

perm_p_unsup <- function(S, g, n = 199) {
  obs <- pillai(S, g)["trace"]
  (1 + sum(sapply(seq_len(n), function(b) pillai(S, sample(g))["trace"]) >= obs)) / (1 + n)
}

res <- data.frame()
for (ver in c("before", "combat")) {
  Xt <- if (ver == "before") Xt_before else Xt_combat
  for (co in levels(cohort_subj)) {
    idx <- which(cohort_subj == co)
    g <- droplevels(grp_subj[idx])
    tp <- tpca(Xt[idx, , , drop = FALSE], ncomp = 5)$variates
    man2 <- pillai(tp[, 1:2], g); man5 <- pillai(tp[, 1:5], g)
    p5 <- perm_p_unsup(tp[, 1:5], g)
    # per-component age LM with FDR (5 comps)
    pv <- sapply(1:5, function(k) t.test(tp[, k] ~ g)$p.value)
    res <- rbind(res, data.frame(
      version = ver, cohort = co, n = length(idx),
      pillai2 = man2["trace"], p2 = man2["p"],
      pillai5 = man5["trace"], p5_parametric = man5["p"], p5_perm = p5,
      best_comp_p = min(pv), best_comp = which.min(pv),
      comp2_p = pv[2], comp5_p = pv[5]
    ))
  }
}

# age distributions
dis <- samples[!duplicated(subject_id), .(subject_id, cohort, age_group, age)]
dis[, age := as.numeric(age)]
dist_tab <- dis[, .(n = .N, mean_age = round(mean(age), 1), min_age = min(age),
                    median_age = as.numeric(median(age)), max_age = max(age)),
                by = .(cohort, age_group)]

print(res)
print(dist_tab)
fwrite(res, file.path(tab_dir, "vacc_heterogeneity_comps.csv"))
fwrite(dist_tab, file.path(tab_dir, "vacc_age_distribution.csv"))
cat("\nDONE\n")
