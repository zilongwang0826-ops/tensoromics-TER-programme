# 21c_vaccine_chip_batch.R
# Deeper batch structure: CEL plate/chip serials are nested in cohort and confounded with age.
# Correct at the finest level (chip) and re-test the age association.
suppressPackageStartupMessages({
  library(data.table); library(tensorOmics); library(sva)
})
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_tensor_helpers.R"))

base_dir <- Sys.getenv("TER_ROOT", getwd())
prep <- readRDS(file.path(base_dir, "data/GSE79396/preprocessed_A.rds"))
dat <- readRDS(file.path(base_dir, "data/GSE79396/vaccine_data.rds"))
samples <- fread(file.path(base_dir, "data/GSE79396/samples_all.csv"))
tab_dir <- file.path(base_dir, "results/tables")
set.seed(2027)

# chip from CEL filename
series_matrix <- Sys.getenv("GSE79396_SERIES_MATRIX",
                            file.path(base_dir, "data/GSE79396/GSE79396_series_matrix.txt.gz"))
if (!file.exists(series_matrix)) {
  stop("GSE79396 series matrix not found: ", series_matrix,
       "\nSet GSE79396_SERIES_MATRIX or place the file under data/GSE79396/.")
}
supp_lines <- readLines(gzfile(series_matrix), n = 200)
supp_files <- supp_lines[startsWith(supp_lines, "!Sample_supplementary_file")]
geo_files <- supp_lines[startsWith(supp_lines, "!Sample_geo_accession")]
gsms <- gsub('"', "", strsplit(geo_files, "\t")[[1]][-1])
cel <- strsplit(supp_files, "\t")[[1]][-1]
chip_by_gsm <- setNames(sub("^.*_([0-9]{10})_.*$", "\\1", cel), gsms)
missing_gsm <- setdiff(samples$gsm, names(chip_by_gsm))
if (length(missing_gsm)) stop("21c: GSM(s) missing from series matrix: ",
                            paste(head(missing_gsm, 5), collapse = ", "))
samples$chip <- unname(chip_by_gsm[samples$gsm])

cat("chip x cohort:\n"); print(table(samples$chip, samples$cohort))
cat("\nchip x age_group:\n"); print(table(samples$chip, samples$age_group))
cat("\nchip x cohort x age (subjects):\n")
print(table(samples[!duplicated(subject_id), .(chip, cohort, age_group)]))
source(file.path(Sys.getenv("TER_ROOT", getwd()), "code", "pipeline", "lib_p2a_helpers.R"))
changed <- write_chip_idempotent(samples, samples$chip,
                                 file.path(base_dir, "data/GSE79396/samples_all.csv"))
cat(if (changed) "samples_all.csv updated (chip column)\n" else
    "samples_all.csv unchanged (chip column identical) - write skipped\n")

subs <- prep$subs
days <- prep$days
grp_subj <- factor(samples$age_group[match(subs, samples$subject_id)], levels = c("young", "elderly"))
cohort_subj <- factor(samples$cohort[match(subs, samples$subject_id)])
chip_row <- factor(samples$chip); chip_subj <- factor(samples$chip[match(subs, samples$subject_id)])
grp_row <- factor(samples$age_group, levels = c("young", "elderly"))
sex_row <- factor(samples$sex)

sel <- prep$sel
raw2000 <- t(dat$expr[sel, ])

build_tensor_scaled <- function(Xrow) {
  Z <- scale(Xrow)
  Xt <- array(NA_real_, dim = c(length(subs), ncol(Z), length(days)),
              dimnames = list(subs, colnames(Z), paste0("D", days)))
  for (s in subs) for (d in days) Xt[s, , paste0("D", d)] <- Z[samples$subject_id == s & samples$day == d, ]
  Xt
}

mod <- model.matrix(~ grp_row + sex_row)
combat_chip <- ComBat(dat = t(raw2000), batch = chip_row, mod = mod)
Xt_chip <- build_tensor_scaled(t(combat_chip))

eval_age <- function(Xt, nperm = 299) {
  tp <- tpca(Xt, ncomp = 2)$variates
  ft <- tplsda(Xt, as.numeric(grp_subj == "elderly"), ncomp = 2)$x_projected
  obs_tp <- pillai(tp, grp_subj)["trace"]
  obs_tp5 <- NULL
  obs_ft <- pillai(ft, grp_subj)["trace"]
  p_tp <- (1 + sum(sapply(seq_len(nperm), function(b) pillai(tp, sample(grp_subj))["trace"]) >= obs_tp)) / (1 + nperm)
  p_ft <- (1 + sum(sapply(seq_len(nperm), function(b) {
    g <- sample(grp_subj); pillai(tplsda(Xt, as.numeric(g == "elderly"), ncomp = 2)$x_projected, g)["trace"]
  }) >= obs_ft)) / (1 + nperm)
  chip_tp <- pillai(tp, chip_subj)
  data.frame(age_pillai_tPCA = obs_tp, age_p_tPCA = pillai(tp, grp_subj)["p"], age_p_tPCA_perm = p_tp,
             age_pillai_tPLSDA = obs_ft, age_p_tPLSDA_perm = p_ft,
             chip_pillai_tPCA = chip_tp["trace"], chip_p_tPCA = chip_tp["p"])
}

res_before <- eval_age(prep$X_tensor)
res_chip   <- eval_age(Xt_chip)

# within-chip 1652 (mixed cohort, both age groups)
idx <- which(chip_subj == "8000001652")
g <- droplevels(grp_subj[idx])
tp1652 <- tpca(prep$X_tensor[idx, , , drop = FALSE], ncomp = 2)$variates
man1652 <- pillai(tp1652, g)
cat("\nchip 1652 (n=", length(idx), ", young=", sum(g=="young"), ", elderly=", sum(g=="elderly"), "): Pillai ", man1652["trace"], " p ", man1652["p"], "\n", sep="")

out <- rbind(cbind(version = "before (none)", res_before),
             cbind(version = "ComBat by chip", res_chip))
fwrite(out, file.path(tab_dir, "vacc_chip_correction.csv"))
fwrite(data.frame(chip1652_n = length(idx), young = sum(g == "young"), elderly = sum(g == "elderly"),
                  pillai = man1652["trace"], p = man1652["p"]),
       file.path(tab_dir, "vacc_chip1652_age.csv"))
print(out)
saveRDS(list(combat_chip_rows = t(combat_chip), chip = samples$chip),
        file.path(base_dir, "data/GSE79396/chip_correction.rds"))
cat("\nDONE\n")
