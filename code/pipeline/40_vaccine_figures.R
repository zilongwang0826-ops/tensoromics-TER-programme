# 40_vaccine_figures.R : report figures for the vaccine dataset (WP2 + WP3)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
base_dir <- Sys.getenv("TER_ROOT", getwd())
tab_dir <- file.path(base_dir, "results/tables")
fig_dir <- file.path(Sys.getenv("TER_ROOT", getwd()), "results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

col_age <- c(young = "#4477AA", elderly = "#CC6677")
col_co  <- c(Atlanta = "#228833", Denver = "#EE7733")
col_chip <- c("8000000153" = "#66CCEE", "8000001652" = "#AA3377", "8000013046" = "#BBBBBB")

row_s <- fread(file.path(tab_dir, "vacc_scores_row.csv"))
sub_s <- fread(file.path(tab_dir, "vacc_scores_subject.csv"))
perm  <- fread(file.path(tab_dir, "vacc_perm_inflation.csv"))
loso  <- fread(file.path(tab_dir, "vacc_loso.csv"))

# ---------------- Figure V1: subject-level representation ---------------------
pA <- ggplot(row_s, aes(mPCA1, mPCA2, group = subject)) +
  geom_line(color = "grey80", linewidth = 0.25, alpha = 0.5) +
  geom_point(aes(color = age_group), size = 1.1, alpha = 0.8) +
  scale_color_manual(values = col_age) +
  labs(title = "(a) matrix PCA — 288 sample-time points", subtitle = "same subject = 4 points (linked)",
       x = "PC1", y = "PC2", color = "Age group") + theme_bw(base_size = 10)

pB <- ggplot(sub_s, aes(tPCA1, tPCA2)) +
  geom_point(aes(color = age_group, shape = cohort), size = 2.4) +
  scale_color_manual(values = col_age) + scale_shape_manual(values = c(16, 17)) +
  labs(title = "(b) tPCA — 72 subject trajectories", subtitle = "one point per participant",
       x = "Comp 1", y = "Comp 2", color = "Age group", shape = "Cohort") + theme_bw(base_size = 10)

pC <- ggplot(row_s, aes(mPLS1, mPLS2, group = subject)) +
  geom_line(color = "grey80", linewidth = 0.25, alpha = 0.5) +
  geom_point(aes(color = age_group), size = 1.1, alpha = 0.8) +
  scale_color_manual(values = col_age) +
  labs(title = "(c) matrix PLS-DA — 288 rows", subtitle = "labels fitted on repeated rows",
       x = "Comp 1", y = "Comp 2", color = "Age group") + theme_bw(base_size = 10)

pD <- ggplot(sub_s, aes(tPLS1, tPLS2)) +
  geom_point(aes(color = age_group, shape = cohort), size = 2.4) +
  scale_color_manual(values = col_age) + scale_shape_manual(values = c(16, 17)) +
  labs(title = "(d) tPLS-DA — 72 subject trajectories", subtitle = "supervised tensor components",
       x = "Comp 1", y = "Comp 2", color = "Age group", shape = "Cohort") + theme_bw(base_size = 10)

figV1 <- (pA | pB) / (pC | pD)
ggsave(file.path(fig_dir, "FV1_fragmentation.png"), figV1, width = 10, height = 8, dpi = 200)

# ---------------- Figure V2: inference inflation + honest CV -------------------
perm[, p_lab := sprintf("%.3f", p_perm)]
perm[, analysis2 := factor(analysis, levels = c("matrix PCA", "subject-averaged matrix PCA", "matrix PLS-DA", "tPLS-DA"), labels = c("matrix PCA", "subject-avg PCA", "matrix PLS-DA", "tPLS-DA"))]
pE <- ggplot(perm, aes(analysis2, -log10(p_perm), fill = scheme)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red") +
  geom_text(aes(label = p_lab), position = position_dodge(0.7), vjust = -0.3, size = 3) +
  labs(title = "(a) Permutation scheme changes significance",
       subtitle = "matrix PLS-DA: p 0.002 (naive rows) vs 0.094 (subject-level)",
       x = NULL, y = expression(-log[10]~p), fill = "Permutation") +
  theme_bw(base_size = 10) + scale_fill_manual(values = c("row-level (naive)" = "#EE8866", "subject-level" = "#99DDFF")) + theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 8))

loso[, analysis2 := factor(analysis, levels = loso$analysis, labels = c("matrix PLS-DA — row-level (naive)", "matrix PLS-DA — subject-grouped", "tPLS-DA — subject-grouped"))]
pF <- ggplot(loso, aes(analysis2, balanced_accuracy, fill = analysis2)) +
  geom_col(width = 0.6) + geom_hline(yintercept = 0.5, linetype = 2) +
  geom_text(aes(label = sprintf("%.3f", balanced_accuracy)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("#EE8866", "#99DDFF", "#44AA99")) +
  coord_cartesian(ylim = c(0.4, 0.78)) +
  labs(title = "(b) Leave-one-out CV: leakage vs subject-grouped",
       subtitle = "row-level CV inflates accuracy; subject-grouped CV is lower",
       x = NULL, y = "Balanced accuracy") +
  theme_bw(base_size = 10) + theme(legend.position = "none",
                                   axis.text.x = element_text(angle = 15, hjust = 1, size = 7.5))

figV2 <- (pE | pF)
ggsave(file.path(fig_dir, "FV2_inflation_cv.png"), figV2, width = 11.5, height = 4.6, dpi = 200)

# ---------------- Figure V3: batch correction ladder + v1 genes ----------------
bs <- fread(file.path(tab_dir, "vacc_batch_scores.csv"))
pG <- ggplot(bs, aes(before_tPCA1, before_tPCA2)) +
  geom_point(aes(color = cohort, shape = age_group), size = 2.2) +
  scale_color_manual(values = col_co) + scale_shape_manual(values = c(16, 17)) +
  labs(title = "(a) tPCA before correction: cohort axis", x = "Comp 1", y = "Comp 2",
       color = "Cohort", shape = "Age group") + theme_bw(base_size = 10)

pH <- ggplot(bs, aes(combat_tPCA1, combat_tPCA2)) +
  geom_point(aes(color = age_group, shape = cohort), size = 2.2) +
  scale_color_manual(values = col_age) + scale_shape_manual(values = c(16, 17)) +
  labs(title = "(b) tPCA after plate-level ComBat", x = "Comp 1", y = "Comp 2",
       color = "Age group", shape = "Cohort") + theme_bw(base_size = 10)

ladder <- data.table(
  test = c("pooled, no correction", "pooled, cohort ComBat", "pooled, plate ComBat",
           "pooled, plate covariate", "Denver only", "Atlanta only", "single plate (n=24)"),
  p = c(0.0033, 0.0033, 0.1767, 0.2951, 0.005, 0.43, 0.158),
  group = c("pooled", "pooled", "pooled", "pooled", "within", "within", "within"))
ladder[, test := factor(test, levels = rev(test))]
pI <- ggplot(ladder, aes(test, -log10(p), fill = group)) +
  geom_col(width = 0.65) + geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red") +
  geom_text(aes(label = sprintf("%.3g", p)), hjust = -0.08, size = 2.8) +
  coord_flip() + scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) + scale_fill_manual(values = c(pooled = "#88CCEE", within = "#DDCC77")) +
  labs(title = "(c) The age association ladder", subtitle = "tPCA / audit age-group test (permutation p where available)",
       x = NULL, y = expression(-log[10]~p)) +
  theme_bw(base_size = 10) + theme(legend.position = "none")

gt <- fread(file.path(tab_dir, "vacc_v1genes_chip_adjust.csv"))
gt[, symbol := ifelse(is.na(symbol) | symbol == "", probe, symbol)]
g2 <- melt(gt[, .(symbol, `age only` = p_age_only, `age + plate + sex` = p_age_chip_sex)], id.vars = "symbol")
g2[, symbol := factor(symbol, levels = gt$symbol)]
pJ <- ggplot(g2, aes(symbol, -log10(value), color = variable, group = symbol)) +
  geom_line(color = "grey70") + geom_point(size = 2.2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red") +
  scale_color_manual(values = c("age only" = "#EE6677", "age + plate + sex" = "#4477AA")) +
  labs(title = "(d) v1 inflammaging genes lose age association", subtitle = "mixed models across 288 sample-days",
       x = NULL, y = expression(-log[10]~p), color = "Model") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

figV3 <- (pG | pH) / (pI | pJ)
ggsave(file.path(fig_dir, "FV3_batch_ladder.png"), figV3, width = 11, height = 9, dpi = 200)

# ---------------- Figure V4: multi-block integration ---------------------------
ms <- fread(file.path(tab_dir, "vacc_multiblock_scores.csv"))
ms[, chip := factor(chip)]
pK <- ggplot(ms, aes(can_tx1, can_mb1)) +
  geom_point(aes(color = chip, shape = cohort), size = 2.4) +
  scale_color_manual(values = col_chip) + scale_shape_manual(values = c(16, 17)) +
  labs(title = "(a) tPLS canonical: shared axis", subtitle = "r = 0.833; pairing permutation p = 0.01",
       x = "Transcriptome comp 1", y = "Metabolome comp 1", color = "Plate", shape = "Cohort") +
  theme_bw(base_size = 10)

pL <- ggplot(ms, aes(can_tx1, can_mb1)) +
  geom_point(aes(color = age_group), size = 2.4) +
  scale_color_manual(values = col_age) +
  labs(title = "(b) same axis coloured by age", subtitle = "alignment dominated by plate/cohort, not age",
       x = "Transcriptome comp 1", y = "Metabolome comp 1", color = "Age group") +
  theme_bw(base_size = 10)

cross <- fread(file.path(tab_dir, "vacc_multiblock_crosslayer.csv"))
pM <- ggplot(cross, aes(factor(day), abs(r))) +
  geom_violin(fill = "#CCBB44", alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.4) +
  labs(title = "(c) cross-layer correlations of top features", subtitle = "top 30 x 30 pairs, comp 1 loadings",
       x = "Day", y = "|Pearson r| across subjects") +
  theme_bw(base_size = 10)

align <- fread(file.path(tab_dir, "vacc_multiblock_align.csv"))
align[, block := sub("^chip-corrected", "plate-corrected", block)]
align[, variable := ifelse(variable == "chip", "plate", variable)]
align[, lab := paste(block, variable)]
align[, lab := factor(lab, levels = align[order(pillai), lab])]
pN <- ggplot(align, aes(lab, pillai, fill = variable)) +
  geom_col(width = 0.6) + coord_flip() + scale_fill_manual(values = c(age = "#4477AA", cohort = "#228833", plate = "#AA3377")) +
  labs(title = "(d) what does the shared axis track?", x = NULL, y = "Pillai trace",
       fill = "Variable") + theme_bw(base_size = 10)

figV4 <- (pK | pL) / (pM | pN)
ggsave(file.path(fig_dir, "FV4_multiblock.png"), figV4, width = 11, height = 9, dpi = 200)
cat("vaccine figures saved\n")
