# 41_fish_figures.R : fish comparison figures
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
base_dir <- Sys.getenv("TER_ROOT", getwd())
dat_dir <- file.path(base_dir, "data/GSE150318")
tab_dir <- file.path(base_dir, "results/tables")
fig_dir <- file.path(Sys.getenv("TER_ROOT", getwd()), "results", "figures")

col_g <- c(short_lived = "#CC6677", long_lived = "#4477AA")
row_s <- fread(file.path(dat_dir, "fish_scores_row.csv"))
sub_s <- fread(file.path(dat_dir, "fish_scores_subject.csv"))

pA <- ggplot(row_s, aes(mPCA1, mPCA2, group = subj)) +
  geom_line(color = "grey80", linewidth = 0.25, alpha = 0.5) +
  geom_point(aes(color = group, shape = factor(time)), size = 1.4) +
  scale_color_manual(values = col_g) + scale_shape_manual(values = c(16, 17), name = "Age (weeks)") +
  labs(title = "(a) matrix PCA — 226 subject-time rows", subtitle = "each fish appears twice (linked)",
       x = "PC1", y = "PC2", color = "Lifespan group") + theme_bw(base_size = 10)

pB <- ggplot(sub_s, aes(tPCA1, tPCA2)) +
  geom_point(aes(color = group), size = 2.4) +
  scale_color_manual(values = col_g) +
  labs(title = "(b) tPCA — 113 fish trajectories", subtitle = "one point per fish (10w+20w)",
       x = "Comp 1", y = "Comp 2", color = "Lifespan group") + theme_bw(base_size = 10)

pC <- ggplot(row_s, aes(mPLS1, mPLS2, group = subj)) +
  geom_line(color = "grey80", linewidth = 0.25, alpha = 0.5) +
  geom_point(aes(color = group, shape = factor(time)), size = 1.4) +
  scale_color_manual(values = col_g) + scale_shape_manual(values = c(16, 17), name = "Age (weeks)") +
  labs(title = "(c) matrix PLS-DA — 226 rows", subtitle = "supervised on lifespan label",
       x = "Comp 1", y = "Comp 2", color = "Lifespan group") + theme_bw(base_size = 10)

pD <- ggplot(sub_s, aes(tPLS1, tPLS2)) +
  geom_point(aes(color = group), size = 2.4) +
  scale_color_manual(values = col_g) +
  labs(title = "(d) tPLS-DA — 113 fish trajectories", subtitle = "Pillai = 0.446 (v1 report reproduction)",
       x = "Comp 1", y = "Comp 2", color = "Lifespan group") + theme_bw(base_size = 10)

figF1 <- (pA | pB) / (pC | pD)
ggsave(file.path(fig_dir, "FF1_fish_fragmentation.png"), figF1, width = 10, height = 8, dpi = 200)

perm <- fread(file.path(tab_dir, "fish_perm_inflation.csv"))
perm[, analysis2 := factor(analysis, levels = c("matrix PCA", "tPCA", "matrix PLS-DA", "tPLS-DA"))]
perm[, p_lab := sprintf("%.3f", p_perm)]
pE <- ggplot(perm, aes(analysis2, -log10(p_perm), fill = scheme)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red") +
  geom_text(aes(label = p_lab), position = position_dodge(0.7), vjust = -0.3, size = 3) +
  labs(title = "(a) fish: permutation-scheme inflation", subtitle = "matrix PLS-DA p 0.008 (naive) vs 0.106 (subject-level)",
       x = NULL, y = expression(-log[10]~p), fill = "Permutation") +
  theme_bw(base_size = 10) + scale_fill_manual(values = c("row-level (naive)" = "#EE8866", "subject-level" = "#99DDFF"))

loso <- fread(file.path(tab_dir, "fish_loso_naive_vs_grouped.csv"))
loso[, analysis2 := factor(analysis, levels = analysis, labels = c("matrix PLS-DA — row-level (naive)", "matrix PLS-DA — subject-grouped", "tPLS-DA — subject-grouped"))]
pF <- ggplot(loso, aes(analysis2, balanced_accuracy, fill = analysis2)) +
  geom_col(width = 0.55) + geom_hline(yintercept = 0.5, linetype = 2) +
  geom_text(aes(label = sprintf("%.3f", balanced_accuracy)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("#EE8866", "#99DDFF", "#44AA99")) +
  coord_cartesian(ylim = c(0.5, 0.76)) +
  labs(title = "(b) fish: LOSO accuracy", subtitle = "row-level CV inflates accuracy; tensor model performs best",
       x = NULL, y = "Balanced accuracy") +
  theme_bw(base_size = 10) + theme(legend.position = "none", axis.text.x = element_text(angle = 12, hjust = 1, size = 7.5))

figF2 <- (pE | pF)
ggsave(file.path(fig_dir, "FF2_fish_inflation_cv.png"), figF2, width = 11, height = 4.8, dpi = 200)
cat("fish figures saved\n")
