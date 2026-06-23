# 06_av_delta.R
# A-V gradient (log2 arterial – log2 venous) by patient covariate.
# Tests whether the magnitude or direction of the A-V delta differs by subgroup.
# Outputs: plots/06_delta_dotplots.png
#          tables/06_delta_results.csv, 06_delta_top_hits.csv

data_dir <- tryCatch(
  normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)
source(file.path(data_dir, "R_scripts", "00_setup.R"))

# ── Compute per-patient A-V delta ─────────────────────────────────────────────
art_s <- soma_samples |> dplyr::filter(draw == "Arterial") |> arrange(patient)
ven_s <- soma_samples |> dplyr::filter(draw == "Venous")   |> arrange(patient)
stopifnot(all(art_s$patient == ven_s$patient))

expr_art  <- expr_all[art_s$SampleId, ]
expr_venm <- expr_all[ven_s$SampleId, ]
delta_mat <- expr_art - expr_venm
rownames(delta_mat) <- art_s$patient

delta_meta <- art_s |> select(patient, all_of(binary_covariates))

# ── Limma on delta per covariate ──────────────────────────────────────────────
run_delta_cov <- function(dmat, dmeta, cov_name) {
  group <- dmeta[[cov_name]]
  keep  <- !is.na(group)
  if (sum(keep) < 4 || length(unique(group[keep])) < 2) return(NULL)
  design_d <- model.matrix(~ group[keep])
  fit_d    <- lmFit(t(dmat[keep, ]), design_d)
  fit_d    <- eBayes(fit_d, trend = TRUE, robust = TRUE)
  topTable(fit_d, coef = 2, number = Inf, sort.by = "P") |>
    rownames_to_column("AptName") |>
    as_tibble() |>
    left_join(analyte_info |> select(AptName, Target, EntrezGeneSymbol), by = "AptName") |>
    mutate(covariate = cov_name,
           n_pos = sum(group[keep] == TRUE,  na.rm = TRUE),
           n_neg = sum(group[keep] == FALSE, na.rm = TRUE))
}

delta_all <- map(binary_covariates, ~run_delta_cov(delta_mat, delta_meta, .x)) |>
  bind_rows()

cat("A-V delta hits (p<0.05) per covariate:\n")
print(delta_all |> dplyr::filter(P.Value < 0.05) |> count(covariate))

write_csv(delta_all, file.path(table_dir, "06_delta_results.csv"))
write_csv(
  delta_all |> dplyr::filter(P.Value < 0.05) |>
    group_by(covariate) |> slice_min(P.Value, n = 10) |> ungroup(),
  file.path(table_dir, "06_delta_top_hits.csv")
)

# ── Dot plots: top 3 per covariate ────────────────────────────────────────────
top_delta <- delta_all |>
  dplyr::filter(P.Value < 0.10) |>
  group_by(covariate) |> slice_min(P.Value, n = 3) |> ungroup()

if (nrow(top_delta) > 0) {
  delta_long <- as_tibble(delta_mat, rownames = "patient") |>
    pivot_longer(-patient, names_to = "AptName", values_to = "delta")
  delta_cov_long <- delta_meta |>
    pivot_longer(all_of(binary_covariates), names_to = "covariate", values_to = "group_val")

  plot_delta <- top_delta |>
    select(covariate, AptName) |>
    left_join(delta_long, by = "AptName") |>
    left_join(delta_cov_long, by = c("patient", "covariate")) |>
    left_join(analyte_info |> select(AptName, Target), by = "AptName") |>
    dplyr::filter(!is.na(group_val)) |>
    mutate(group_label = if_else(group_val, "Yes", "No"),
           protein     = if_else(is.na(Target) | Target == "", AptName, Target),
           panel_label = paste0(covariate, ": ", protein))

  if (nrow(plot_delta) > 0) {
    p_delta <- ggplot(plot_delta, aes(x = group_label, y = delta, color = group_label)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_jitter(width = 0.1, size = 2.5, alpha = 0.8) +
      stat_summary(fun = median, geom = "crossbar", width = 0.4,
                   color = "black", linewidth = 0.5) +
      facet_wrap(~ panel_label, scales = "free_y", ncol = 3) +
      scale_color_manual(values = col_yn) +
      labs(
        title = "A-V delta (log2 arterial − venous) by covariate — hits at p < 0.10",
        x = NULL, y = "log2 A-V delta", color = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none", strip.text = element_text(size = 8))

    ggsave(file.path(plot_dir, "06_delta_dotplots.png"), p_delta,
           width = 12, height = 10, dpi = 300)
  }
} else {
  cat("No hits at p < 0.10 — no dot plot generated.\n")
}

cat("\nOutputs saved to outputs/plots/ and outputs/tables/\n")

# ── A-V delta by surgery type (omnibus F + pairwise contrasts) ────────────────
delta_meta_surg <- art_s |>
  select(patient, surgery_group) |>
  mutate(surgery_group = factor(surgery_group))

surg_levels <- levels(delta_meta_surg$surgery_group)   # Head/Neck, Laparoscopic, Spine
safe_levels <- make.names(surg_levels)                  # Head.Neck, Laparoscopic, Spine

design_surg <- model.matrix(~ 0 + surgery_group, data = delta_meta_surg)
colnames(design_surg) <- safe_levels

contrasts_surg <- makeContrasts(
  Spine_vs_Laparoscopic    = Spine - Laparoscopic,
  Spine_vs_HeadNeck        = Spine - Head.Neck,
  Laparoscopic_vs_HeadNeck = Laparoscopic - Head.Neck,
  levels = safe_levels
)

fit_surg  <- lmFit(t(delta_mat), design_surg)
fit_surg2 <- contrasts.fit(fit_surg, contrasts_surg)
fit_surg2 <- eBayes(fit_surg2, trend = TRUE, robust = TRUE)

# Omnibus F-test: does A-V delta differ across any surgery group?
surg_omnibus <- topTable(fit_surg2, number = Inf, sort.by = "F") |>
  rownames_to_column("AptName") |>
  as_tibble() |>
  left_join(analyte_info |> select(AptName, Target, EntrezGeneSymbol), by = "AptName")

cat("\nA-V delta surgery omnibus hits (F-test p<0.05):", sum(surg_omnibus$P.Value < 0.05), "\n")

write_csv(surg_omnibus, file.path(table_dir, "06_delta_surgery_omnibus.csv"))
write_csv(
  surg_omnibus |> dplyr::filter(P.Value < 0.05) |> slice_min(P.Value, n = 20),
  file.path(table_dir, "06_delta_surgery_omnibus_top.csv")
)

# Pairwise contrasts
surg_pairwise <- map(colnames(contrasts_surg), function(ct) {
  topTable(fit_surg2, coef = ct, number = Inf, sort.by = "P") |>
    rownames_to_column("AptName") |>
    as_tibble() |>
    left_join(analyte_info |> select(AptName, Target, EntrezGeneSymbol), by = "AptName") |>
    mutate(contrast = ct)
}) |> bind_rows()

cat("A-V delta surgery pairwise hits (p<0.05) per contrast:\n")
print(surg_pairwise |> dplyr::filter(P.Value < 0.05) |> count(contrast))

write_csv(surg_pairwise, file.path(table_dir, "06_delta_surgery_pairwise.csv"))
write_csv(
  surg_pairwise |> dplyr::filter(P.Value < 0.05) |>
    group_by(contrast) |> slice_min(P.Value, n = 10) |> ungroup(),
  file.path(table_dir, "06_delta_surgery_pairwise_top.csv")
)

# ── PCA of A-V delta matrix, colored by surgery type ─────────────────────────
pca_delta <- prcomp(delta_mat, center = TRUE, scale. = TRUE)
pca_delta_df <- as_tibble(pca_delta$x[, 1:2], rownames = "patient") |>
  left_join(delta_meta_surg, by = "patient")

var_exp <- round(summary(pca_delta)$importance[2, 1:2] * 100, 1)

p_pca_surg <- ggplot(pca_delta_df, aes(x = PC1, y = PC2, color = surgery_group,
                                        label = patient)) +
  geom_point(size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_manual(values = col_surg) +
  labs(
    title = "PCA of A-V delta (log2 A − V) — whole panel",
    x = paste0("PC1 (", var_exp[1], "% var)"),
    y = paste0("PC2 (", var_exp[2], "% var)"),
    color = "Surgery type"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(plot_dir, "06_delta_surgery_pca.png"), p_pca_surg,
       width = 7, height = 6, dpi = 300)
cat("Surgery PCA saved.\n")

# ── Heatmap: all omnibus-significant proteins ──────────────────────────────────
sig_apts <- surg_omnibus |>
  dplyr::filter(P.Value < 0.05) |>
  pull(AptName)

if (length(sig_apts) > 0) {
  hm_mat <- t(delta_mat[, sig_apts])
  rownames(hm_mat) <- surg_omnibus$Target[match(sig_apts, surg_omnibus$AptName)] |>
    coalesce(sig_apts)

  annot_col <- data.frame(
    Surgery = delta_meta_surg$surgery_group,
    row.names = delta_meta_surg$patient
  )
  annot_colors <- list(Surgery = col_surg)

  pheatmap(
    hm_mat,
    annotation_col  = annot_col,
    annotation_colors = annot_colors,
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    scale           = "row",
    color           = colorRampPalette(c("#2980b9", "white", "#c0392b"))(50),
    fontsize_row    = 7,
    fontsize_col    = 8,
    main            = paste0("A-V delta by surgery type — ", length(sig_apts),
                             " omnibus-significant proteins (p<0.05), row-scaled"),
    border_color    = NA,
    angle_col       = 45,
    filename        = file.path(plot_dir, "06_delta_surgery_heatmap.pdf"),
    width = 10, height = max(6, length(sig_apts) * 0.15 + 3)
  )
  cat("Surgery heatmap saved (", length(sig_apts), "proteins).\n")
} else {
  cat("No omnibus hits at p<0.05 — no heatmap generated.\n")
}

# ── Dot plots: top omnibus hits, all 3 surgery groups ─────────────────────────
top_surg <- surg_omnibus |>
  dplyr::filter(P.Value < 0.10) |>
  slice_min(P.Value, n = 9)

if (nrow(top_surg) > 0) {
  delta_long_s <- as_tibble(delta_mat, rownames = "patient") |>
    pivot_longer(-patient, names_to = "AptName", values_to = "delta") |>
    inner_join(top_surg |> select(AptName, Target), by = "AptName") |>
    left_join(delta_meta_surg, by = "patient") |>
    mutate(protein     = if_else(is.na(Target) | Target == "", AptName, Target),
           panel_label = paste0("F p=", formatC(
             surg_omnibus$P.Value[match(AptName, surg_omnibus$AptName)], format = "g", digits = 2
           ), "  ", protein))

  p_surg <- ggplot(delta_long_s,
                   aes(x = surgery_group, y = delta, color = surgery_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4,
                 color = "black", linewidth = 0.5) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = 3) +
    scale_color_manual(values = col_surg) +
    labs(
      title = "A-V delta (log2 arterial − venous) by surgery type — top omnibus F hits",
      x = NULL, y = "log2 A-V delta", color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", strip.text = element_text(size = 7),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())

  ggsave(file.path(plot_dir, "06_delta_surgery_dotplots.png"), p_surg,
         width = 12, height = 10, dpi = 300)
  cat("Surgery dot plot saved.\n")
} else {
  cat("No omnibus hits at p < 0.10 — no surgery dot plot generated.\n")
}

cat("\nSurgery-type outputs saved to outputs/plots/ and outputs/tables/\n")
