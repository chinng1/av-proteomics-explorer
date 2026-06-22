# 07_platform_concordance.R
# Cross-platform concordance: Luminex 65-plex vs. SomaScan.
#   (a) A-V log2FC correlation across matched proteins (effect-size concordance)
#   (b) Per-patient log2-MFI vs. log2-RFU correlation within arterial/venous
# Outputs: plots/07_logfc_scatter.png, 07_per_patient_scatter.png, 07_concordance_bar.png
#          tables/07_logfc_concordance.csv, 07_per_patient_cor.csv

data_dir <- tryCatch(
  normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)
source(file.path(data_dir, "R_scripts", "00_setup.R"))

# ── Load Luminex ───────────────────────────────────────────────────────────────
plate_map <- data.frame(
  unknown_id = paste0("Unknown", 1:24),
  patient    = c("0179","0193","0209","0215","0179","0193","0209","0215",
                 "0223","0245","0250","0254","0223","0245","0250","0254",
                 "0268","0275","0284","0291","0268","0275","0284","0291"),
  draw       = c(rep("Venous",4), rep("Arterial",4), rep("Venous",4), rep("Arterial",4),
                 rep("Venous",4), rep("Arterial",4))
)

lx_raw <- read_csv(file.path(data_dir, "AP_65plex_byJMR.csv"),
                    skip = 54, show_col_types = FALSE) |>
  dplyr::filter(str_starts(Sample, "Unknown"))
lx_cols <- lx_raw |> select(-Location, -Sample, -`Total Events`) |> names()

lx_avg <- lx_raw |>
  mutate(across(all_of(lx_cols), ~suppressWarnings(as.numeric(.)))) |>
  left_join(plate_map, by = c("Sample" = "unknown_id")) |>
  group_by(patient, draw) |>
  summarise(across(all_of(lx_cols), ~mean(., na.rm = TRUE)), .groups = "drop")

lx_art <- lx_avg |> dplyr::filter(draw == "Arterial") |> arrange(patient)
lx_ven <- lx_avg |> dplyr::filter(draw == "Venous")   |> arrange(patient)

lx_results <- map_dfr(lx_cols, function(cyt) {
  a <- lx_art[[cyt]]; v <- lx_ven[[cyt]]
  test <- tryCatch(wilcox.test(a, v, paired = TRUE, exact = FALSE),
                   error = function(e) list(p.value = NA))
  tibble(luminex_name = cyt,
         log2FC_lx    = log2((mean(a, na.rm=TRUE)+0.1)/(mean(v, na.rm=TRUE)+0.1)),
         p_value_lx   = test$p.value)
}) |> mutate(p_adj_lx = p.adjust(p_value_lx, method = "BH"))

# ── Re-run SomaScan A-V limma (minimal) ───────────────────────────────────────
patient_f <- factor(soma_samples$patient)
draw_f    <- factor(soma_samples$draw, levels = c("Venous", "Arterial"))
design_av <- model.matrix(~ patient_f + draw_f)
fit_av    <- lmFit(t(expr_all), design_av)
fit_av    <- eBayes(fit_av, trend = TRUE, robust = TRUE)
av_results <- topTable(fit_av, coef = "draw_fArterial", number = Inf, sort.by = "P") |>
  rownames_to_column("AptName") |> as_tibble() |>
  left_join(analyte_info |> select(AptName, Target, EntrezGeneSymbol), by = "AptName")

soma_best <- av_results |>
  group_by(EntrezGeneSymbol) |>
  slice_min(P.Value, n = 1, with_ties = FALSE) |>
  ungroup()

# ── (a) A-V log2FC concordance ────────────────────────────────────────────────
platform_df <- luminex_gene_map |>
  left_join(lx_results, by = "luminex_name") |>
  left_join(soma_best |> select(EntrezGeneSymbol,
                                 logFC_soma   = logFC,
                                 P.Value_soma = P.Value,
                                 AptName_soma = AptName),
            by = c("gene_symbol" = "EntrezGeneSymbol")) |>
  mutate(
    in_soma    = !is.na(logFC_soma),
    dir_lx     = if_else(log2FC_lx  > 0, "Art > Ven", "Ven > Art"),
    dir_soma   = if_else(logFC_soma > 0, "Art > Ven", "Ven > Art"),
    concordant = in_soma & (dir_lx == dir_soma),
    sig_lx     = p_value_lx < 0.05,
    sig_soma   = P.Value_soma < 0.05
  )

cat("Matched in SomaScan:", sum(platform_df$in_soma), "of", nrow(platform_df), "\n")
cat("Direction concordant:", sum(platform_df$concordant, na.rm = TRUE), "\n")
r_fc <- cor(platform_df$log2FC_lx, platform_df$logFC_soma,
            use = "complete.obs", method = "pearson")
r_sp <- cor(platform_df$log2FC_lx, platform_df$logFC_soma,
            use = "complete.obs", method = "spearman")
cat("Pearson r of fold-changes:", round(r_fc, 3), "\n")

write_csv(platform_df, file.path(table_dir, "07_logfc_concordance.csv"))

scatter_df <- platform_df |>
  dplyr::filter(in_soma) |>
  mutate(
    sig_category = case_when(
      sig_lx & sig_soma ~ "Both p < 0.05",
      sig_lx            ~ "Luminex only",
      sig_soma          ~ "SomaScan only",
      TRUE              ~ "Neither"
    ),
    label_pt = if_else(sig_lx | sig_soma, gene_symbol, NA_character_)
  )

p_logfc <- ggplot(scatter_df,
       aes(x = log2FC_lx, y = logFC_soma, color = sig_category, label = label_pt)) +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(size = 2.8, max.overlaps = 25, box.padding = 0.3,
                  na.rm = TRUE, show.legend = FALSE) +
  scale_color_manual(values = c("Both p < 0.05" = "#8e44ad", "Luminex only" = "#e67e22",
                                 "SomaScan only" = "#27ae60", "Neither" = "gray60")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3.5, color = "gray20",
           label = sprintf("Pearson r = %.2f\nSpearman ρ = %.2f", r_fc, r_sp)) +
  labs(
    title    = "Platform concordance: Luminex vs. SomaScan A-V log2FC",
    subtitle = sprintf("%d of 65 cytokines matched. Dashed = identity line.",
                       sum(scatter_df$in_soma)),
    x = "log2FC Luminex (Arterial / Venous)", y = "log2FC SomaScan", color = NULL
  ) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "07_logfc_scatter.png"), p_logfc, width = 9, height = 7, dpi = 300)

# ── (b) Per-patient concordance ───────────────────────────────────────────────
apt_mean_rfu <- soma_samples |> select(all_of(pass_analytes)) |>
  summarise(across(everything(), ~mean(., na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "AptName", values_to = "mean_rfu")

soma_gene_best <- analyte_info |> select(AptName, Target, EntrezGeneSymbol) |>
  dplyr::filter(AptName %in% pass_analytes,
                EntrezGeneSymbol %in% luminex_gene_map$gene_symbol) |>
  left_join(apt_mean_rfu, by = "AptName") |>
  group_by(EntrezGeneSymbol) |> slice_max(mean_rfu, n = 1, with_ties = FALSE) |> ungroup()

per_pt_crosswalk <- luminex_gene_map |>
  dplyr::filter(luminex_name %in% lx_cols) |>
  left_join(soma_gene_best |> select(EntrezGeneSymbol, AptName, Target),
            by = c("gene_symbol" = "EntrezGeneSymbol")) |>
  dplyr::filter(!is.na(AptName))

lumi_long <- lx_avg |>
  select(patient, draw, all_of(per_pt_crosswalk$luminex_name)) |>
  pivot_longer(all_of(per_pt_crosswalk$luminex_name),
               names_to = "luminex_name", values_to = "mfi") |>
  mutate(log2_mfi = log2(pmax(mfi, 1)))

soma_long <- soma_samples |>
  select(patient, draw, all_of(per_pt_crosswalk$AptName)) |>
  mutate(across(all_of(per_pt_crosswalk$AptName), log2)) |>
  pivot_longer(all_of(per_pt_crosswalk$AptName),
               names_to = "AptName", values_to = "log2_rfu")

per_pt_df <- per_pt_crosswalk |> select(luminex_name, AptName, Target) |>
  left_join(lumi_long, by = "luminex_name") |>
  left_join(soma_long, by = c("AptName", "patient", "draw")) |>
  mutate(draw = factor(draw, levels = c("Venous", "Arterial")))

per_pt_cor <- per_pt_df |>
  dplyr::filter(!is.na(draw)) |>
  group_by(Target, luminex_name, draw) |>
  summarise(
    n          = sum(!is.na(log2_mfi) & !is.na(log2_rfu)),
    r_spearman = if (n >= 3)
      tryCatch(cor(log2_mfi, log2_rfu, method = "spearman", use = "complete.obs"),
               error = function(e) NA_real_)
    else NA_real_,
    .groups = "drop"
  ) |>
  dplyr::filter(!is.na(r_spearman))

write_csv(per_pt_cor, file.path(table_dir, "07_per_patient_cor.csv"))

# Per-patient scatter grid
p_scatter <- per_pt_df |>
  mutate(panel_label = paste0(Target, "\n(", luminex_name, ")")) |>
  ggplot(aes(x = log2_mfi, y = log2_rfu, color = draw, shape = draw)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  facet_wrap(~ panel_label, scales = "free", ncol = 6) +
  scale_color_manual(values = col_av) +
  scale_shape_manual(values = c("Arterial" = 17, "Venous" = 16)) +
  labs(
    title    = "Per-patient concordance: Luminex log2-MFI vs. SomaScan log2-RFU",
    subtitle = "Each point = one patient. Lines = per-draw linear fit.",
    x = "Luminex log2-MFI (negatives floored at 1)",
    y = "SomaScan log2-RFU", color = NULL, shape = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom", strip.text = element_text(size = 6))

ggsave(file.path(plot_dir, "07_per_patient_scatter.png"), p_scatter,
       width = 18, height = 22, dpi = 300)

# Concordance bar chart
p_bar <- per_pt_cor |>
  mutate(Target = fct_reorder(Target, r_spearman, .fun = mean, .na_rm = TRUE)) |>
  ggplot(aes(x = Target, y = r_spearman, fill = draw)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = c(0, 0.7), linetype = c("solid", "dashed"),
             color = c("gray30", "gray60")) +
  scale_fill_manual(values = col_av) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  coord_flip() +
  labs(
    title    = "Per-patient Spearman r: Luminex log2-MFI vs. SomaScan log2-RFU",
    subtitle = "Dashed = r = 0.70. Sorted by mean r.",
    x = NULL, y = "Spearman r (N = 12 patients per draw)", fill = NULL
  ) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "07_concordance_bar.png"), p_bar, width = 10, height = 8, dpi = 300)

cat("\nOutputs saved to outputs/plots/ and outputs/tables/\n")
