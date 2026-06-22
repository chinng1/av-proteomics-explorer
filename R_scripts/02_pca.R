# 02_pca.R
# PCA on log2-RFU expression matrix; ANOVA testing draw and surgery effects on PCs.
# Outputs: plots/02_pca_pc1_pc2.png, 02_pca_pc3_pc4.png, 02_pca_scores_by_group.png
#          tables/02_pca_variance.csv, 02_pca_loadings.csv, 02_pca_anova.csv

data_dir <- tryCatch(
  normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)
source(file.path(data_dir, "R_scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(lme4); library(car) })

# ── PCA ────────────────────────────────────────────────────────────────────────
pca_res <- prcomp(expr_all, scale. = TRUE)
pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as_tibble(pca_res$x[, 1:4]) |>
  mutate(
    SampleId      = soma_samples$SampleId,
    patient       = soma_samples$patient,
    draw          = factor(soma_samples$draw, levels = c("Venous", "Arterial")),
    surgery_group = factor(soma_samples$surgery_group),
    flagged       = soma_samples$flagged,
    male          = soma_samples$male,
    obese         = soma_samples$obese,
    hypothermia   = soma_samples$hypothermia,
    hyperglycemia = soma_samples$hyperglycemia
  )

# Variance explained table
var_df <- tibble(
  PC      = paste0("PC", seq_along(pca_res$sdev)),
  sdev    = pca_res$sdev,
  var_pct = round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2),
  cum_pct = cumsum(round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2))
)
write_csv(var_df |> head(20), file.path(table_dir, "02_pca_variance.csv"))

# ── Loadings ───────────────────────────────────────────────────────────────────
loadings_df <- as_tibble(pca_res$rotation[, 1:4], rownames = "AptName") |>
  left_join(analyte_info |> select(AptName, Target, EntrezGeneSymbol), by = "AptName") |>
  arrange(desc(abs(PC1)))
write_csv(loadings_df, file.path(table_dir, "02_pca_loadings.csv"))

# PCA sample scores — used by shiny_explorer
write_csv(pca_df, file.path(table_dir, "02_pca_scores.csv"))

# ── PC1 vs PC2 plot ────────────────────────────────────────────────────────────
p1 <- ggplot(pca_df, aes(PC1, PC2, color = draw, shape = surgery_group, label = patient)) +
  geom_line(aes(group = patient), color = "gray70", linewidth = 0.5) +
  geom_point(aes(size = flagged), alpha = 0.9) +
  geom_text_repel(size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = col_av) +
  scale_shape_manual(values = c("Head/Neck" = 17, "Laparoscopic" = 15, "Spine" = 16)) +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5), guide = "none") +
  labs(
    title    = "PCA: SomaScan 7,481 proteins (log2 RFU, scaled)",
    subtitle = "Lines connect arterial–venous pairs. Large points = RowCheck FLAG.",
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)"),
    color = "Draw", shape = "Surgery"
  ) +
  theme_bw(base_size = 13) + theme(legend.position = "right")

ggsave(file.path(plot_dir, "02_pca_pc1_pc2.png"), p1, width = 9, height = 7, dpi = 300)

# ── PC3 vs PC4 plot ────────────────────────────────────────────────────────────
p2 <- ggplot(pca_df, aes(PC3, PC4, color = draw, shape = surgery_group, label = patient)) +
  geom_line(aes(group = patient), color = "gray70", linewidth = 0.5) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = col_av) +
  scale_shape_manual(values = c("Head/Neck" = 17, "Laparoscopic" = 15, "Spine" = 16)) +
  labs(
    x = paste0("PC3 (", pct_var[3], "%)"),
    y = paste0("PC4 (", pct_var[4], "%)"),
    color = "Draw", shape = "Surgery"
  ) +
  theme_bw(base_size = 13) + theme(legend.position = "right")

ggsave(file.path(plot_dir, "02_pca_pc3_pc4.png"), p2, width = 9, height = 7, dpi = 300)

# ── ANOVA: draw and surgery effects on PC1–PC4 ────────────────────────────────
fit_pc <- function(pc_col) {
  form <- as.formula(paste0(pc_col, " ~ draw + surgery_group + (1 | patient)"))
  m    <- lmer(form, data = pca_df, REML = FALSE)
  anov <- car::Anova(m, type = "III")
  tibble(
    PC        = pc_col,
    term      = rownames(anov),
    Chisq     = round(anov$Chisq, 3),
    df        = anov$Df,
    p_value   = signif(anov$`Pr(>Chisq)`, 3)
  ) |> dplyr::filter(term != "(Intercept)")
}

anova_results <- map_dfr(c("PC1", "PC2", "PC3", "PC4"), fit_pc)
cat("Mixed-model ANOVA (Type III):\n")
print(anova_results)
write_csv(anova_results, file.path(table_dir, "02_pca_anova.csv"))

# ── PC scores by group dot plot ────────────────────────────────────────────────
p3 <- pca_df |>
  pivot_longer(c("PC1","PC2","PC3","PC4"), names_to = "PC", values_to = "score") |>
  mutate(PC = factor(PC, levels = c("PC1","PC2","PC3","PC4"))) |>
  ggplot(aes(x = surgery_group, y = score, color = draw, shape = draw)) +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_jitter(width = 0.15, size = 3, alpha = 0.85) +
  stat_summary(aes(group = draw), fun = mean, geom = "crossbar",
               width = 0.4, linewidth = 0.5, show.legend = FALSE) +
  facet_wrap(~ PC, scales = "free_y", ncol = 2) +
  scale_color_manual(values = col_av) +
  scale_shape_manual(values = c("Arterial" = 17, "Venous" = 16)) +
  labs(
    title    = "PC scores by surgery type and draw",
    subtitle = "Crossbars = group mean. Mixed-model ANOVA p-values in table 02_pca_anova.csv.",
    x = NULL, y = "PC score", color = NULL, shape = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(plot_dir, "02_pca_scores_by_group.png"), p3, width = 10, height = 7, dpi = 300)

cat("\nOutputs saved to outputs/plots/ and outputs/tables/\n")
