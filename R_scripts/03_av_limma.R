# 03_av_limma.R
# Paired limma: arterial vs. venous across 7,481 SomaScan proteins.
# Outputs: plots/03_volcano.png, 03_ma_plot.png, 03_paired_dots.png, 03_heatmap_top50.pdf
#          tables/03_av_results_full.csv, 03_av_results_nominal.csv

data_dir <- tryCatch(
  normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)
source(file.path(data_dir, "R_scripts", "00_setup.R"))

# ── Paired limma ───────────────────────────────────────────────────────────────
patient_f <- factor(soma_samples$patient)
draw_f    <- factor(soma_samples$draw, levels = c("Venous", "Arterial"))
design_av <- model.matrix(~ patient_f + draw_f)

fit_av <- lmFit(t(expr_all), design_av)
fit_av <- eBayes(fit_av, trend = TRUE, robust = TRUE)

av_results <- topTable(fit_av, coef = "draw_fArterial", number = Inf, sort.by = "P") |>
  rownames_to_column("AptName") |>
  as_tibble() |>
  left_join(
    analyte_info |> select(AptName, Target, TargetFullName, EntrezGeneSymbol,
                            EntrezGeneID, UniProt, Dilution),
    by = "AptName"
  ) |>
  mutate(
    direction = if_else(logFC > 0, "Arterial > Venous", "Venous > Arterial"),
    sig_nom   = P.Value   < 0.05,
    sig_fdr   = adj.P.Val < 0.05
  ) |>
  arrange(P.Value)

cat("Nominal p<0.05:", sum(av_results$sig_nom), "/ FDR<0.05:", sum(av_results$sig_fdr), "\n")
cat("Direction (nominal hits):\n")
print(av_results |> dplyr::filter(sig_nom) |> count(direction))

write_csv(av_results, file.path(table_dir, "03_av_results_full.csv"))
write_csv(av_results |> dplyr::filter(sig_nom),
          file.path(table_dir, "03_av_results_nominal.csv"))

# ── Volcano plot ───────────────────────────────────────────────────────────────
label_df <- av_results |> slice_min(P.Value, n = 25)

p_volcano <- ggplot(av_results,
       aes(x = logFC, y = -log10(P.Value), color = direction, alpha = sig_nom)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_point(size = 1.2) +
  geom_text_repel(data = label_df, aes(label = Target),
                  size = 2.5, max.overlaps = 30, box.padding = 0.35,
                  show.legend = FALSE) +
  scale_color_manual(values = c("Arterial > Venous" = "#c0392b",
                                 "Venous > Arterial" = "#2980b9")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.25), guide = "none") +
  labs(
    title    = "SomaScan: Arterial vs. Venous (paired limma, N = 12, 7,481 proteins)",
    subtitle = "Dashed = nominal p < 0.05. Top 25 labeled.",
    x = "log2FC (Arterial / Venous)", y = "-log10(p-value)", color = NULL
  ) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "03_volcano.png"), p_volcano, width = 10, height = 8, dpi = 300)

# ── MA plot ────────────────────────────────────────────────────────────────────
p_ma <- ggplot(av_results, aes(x = AveExpr, y = logFC, color = sig_nom)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = 0, color = "gray40") +
  geom_text_repel(data = av_results |> slice_min(P.Value, n = 20),
                  aes(label = Target), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "#c0392b"),
                     labels = c("p ≥ 0.05", "p < 0.05")) +
  labs(
    title = "MA plot: SomaScan A-V",
    x = "Average log2 expression (AveExpr)",
    y = "log2FC (Arterial / Venous)", color = NULL
  ) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "03_ma_plot.png"), p_ma, width = 10, height = 7, dpi = 300)

# ── Paired trajectory plot — top 12 ───────────────────────────────────────────
top_apts <- av_results |> slice_min(P.Value, n = 12) |> pull(AptName)
top_tgts <- av_results |> slice_min(P.Value, n = 12) |>
  mutate(label = if_else(is.na(Target) | Target == "", AptName, Target)) |>
  pull(label)
names(top_tgts) <- top_apts

p_paired <- soma_samples |>
  select(patient, draw, surgery_group, all_of(top_apts)) |>
  mutate(across(all_of(top_apts), log2)) |>
  pivot_longer(all_of(top_apts), names_to = "AptName", values_to = "log2_rfu") |>
  mutate(
    draw    = factor(draw, levels = c("Venous", "Arterial")),
    protein = factor(top_tgts[AptName], levels = top_tgts)
  ) |>
  ggplot(aes(x = draw, y = log2_rfu, group = patient, color = surgery_group)) +
  geom_line(alpha = 0.5) + geom_point(size = 2) +
  facet_wrap(~ protein, scales = "free_y", ncol = 4) +
  scale_color_manual(values = col_surg) +
  labs(title = "Top 12 proteins: paired arterial vs. venous per patient",
       x = NULL, y = "log2 RFU", color = "Surgery") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "03_paired_dots_top12.png"), p_paired,
       width = 12, height = 8, dpi = 300)

# ── Spotlight: 4 named proteins (SPINK9, PDXK, KLK7, transgelin) ─────────────
spotlight_genes <- c("SPINK9", "PDXK", "KLK7", "TAGLN")
spotlight_df <- av_results |>
  dplyr::filter(EntrezGeneSymbol %in% spotlight_genes) |>
  group_by(EntrezGeneSymbol) |>
  slice_min(P.Value, n = 1) |>
  ungroup() |>
  arrange(P.Value)

if (nrow(spotlight_df) > 0) {
  sp_apts <- spotlight_df$AptName

  p_spotlight <- soma_samples |>
    select(patient, draw, surgery_group, all_of(sp_apts)) |>
    mutate(across(all_of(sp_apts), log2)) |>
    pivot_longer(all_of(sp_apts), names_to = "AptName", values_to = "log2_rfu") |>
    left_join(spotlight_df |> select(AptName, EntrezGeneSymbol, P.Value, logFC),
              by = "AptName") |>
    mutate(
      draw        = factor(draw, levels = c("Venous", "Arterial")),
      panel_label = sprintf("%s\n(logFC = %.2f, p = %.4f)",
                            EntrezGeneSymbol, logFC, P.Value),
      panel_label = factor(panel_label,
                           levels = unique(panel_label[order(spotlight_df$P.Value[
                             match(AptName, spotlight_df$AptName)])]))
    ) |>
    ggplot(aes(x = draw, y = log2_rfu, group = patient, color = surgery_group)) +
    geom_line(alpha = 0.6, linewidth = 0.8) +
    geom_point(size = 3) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = col_surg) +
    labs(
      title    = "Representative top hits: paired arterial vs. venous per patient",
      subtitle = "Lines connect the same patient. Downward slope = higher in venous blood.",
      x        = NULL, y = "log2 RFU", color = "Surgery"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", strip.text = element_text(size = 10))

  ggsave(file.path(plot_dir, "03_spotlight_top4.png"), p_spotlight,
         width = 8, height = 7, dpi = 300)
}

# ── Heatmap — top 50 ──────────────────────────────────────────────────────────
top50_apts <- av_results |> slice_min(P.Value, n = 50) |> pull(AptName)
hm_mat <- t(expr_all[, top50_apts])
colnames(hm_mat) <- soma_samples$SampleId
rownames(hm_mat) <- av_results |>
  dplyr::filter(AptName %in% top50_apts) |>
  arrange(match(AptName, top50_apts)) |>
  mutate(label = if_else(is.na(Target) | Target == "", AptName, Target)) |>
  pull(label)

col_anno <- data.frame(
  Draw    = soma_samples$draw,
  Surgery = soma_samples$surgery_group,
  Flagged = ifelse(soma_samples$flagged, "Yes", "No"),
  row.names = soma_samples$SampleId
)
anno_colors <- list(
  Draw    = col_av,
  Surgery = col_surg,
  Flagged = c("Yes" = "orange", "No" = "white")
)

pheatmap(
  hm_mat,
  scale             = "row",
  annotation_col    = col_anno,
  annotation_colors = anno_colors,
  show_colnames     = FALSE,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  fontsize_row      = 7,
  main              = "Top 50 A-V proteins (row-scaled log2 RFU)",
  filename          = file.path(plot_dir, "03_heatmap_top50.pdf"),
  width = 10, height = 14
)

cat("\nOutputs saved to outputs/plots/ and outputs/tables/\n")
