# 05_covariate_analysis.R
# Limma covariate analysis on arterial and venous samples independently.
# Covariates: obesity, hypothermia, hyperglycemia, sex (binary) + surgery type (3-level).
# Outputs:
#   tables/05_{arterial,venous}_results.csv
#   tables/05_{arterial,venous}_top_hits.csv
#   tables/05_{arterial,venous}_surgery_omnibus.csv
#   tables/05_{arterial,venous}_surgery_pairwise.csv
#   plots/05_{arterial,venous}_covariate_heatmap.pdf
#   plots/05_{arterial,venous}_surgery_heatmap.pdf
#   plots/05_{arterial,venous}_pca_{covariate}.png
#   plots/05_{arterial,venous}_surgery_pca.png

data_dir <- tryCatch(
  normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)
source(file.path(data_dir, "R_scripts", "00_setup.R"))

# ── Split samples by draw ──────────────────────────────────────────────────────
art_idx  <- soma_samples$draw == "Arterial"
ven_idx  <- soma_samples$draw == "Venous"
soma_art <- soma_samples[art_idx, ]
soma_ven <- soma_samples[ven_idx, ]
expr_art <- expr_all[art_idx, ]
expr_ven <- expr_all[ven_idx, ]

# ── Binary covariate limma ─────────────────────────────────────────────────────
run_covariate_limma <- function(expr_mat, sample_df, cov_name, draw_label) {
  group <- sample_df[[cov_name]]
  keep  <- !is.na(group)
  if (sum(keep) < 4 || length(unique(group[keep])) < 2) {
    message("Skipping ", draw_label, " / ", cov_name, ": insufficient observations")
    return(NULL)
  }
  design_cov <- model.matrix(~ group[keep])
  fit_cov    <- lmFit(t(expr_mat[keep, ]), design_cov)
  fit_cov    <- eBayes(fit_cov, trend = TRUE, robust = TRUE)
  topTable(fit_cov, coef = 2, number = Inf, sort.by = "P") |>
    rownames_to_column("AptName") |>
    as_tibble() |>
    left_join(analyte_info |> select(AptName, Target, EntrezGeneSymbol, UniProt),
              by = "AptName") |>
    mutate(covariate = cov_name, draw = draw_label,
           n_pos = sum(group[keep] == TRUE,  na.rm = TRUE),
           n_neg = sum(group[keep] == FALSE, na.rm = TRUE))
}

# ── Surgery type limma (omnibus F + pairwise) ──────────────────────────────────
run_surgery_limma <- function(expr_mat, sample_df, draw_label) {
  surg <- factor(sample_df$surgery_group)
  safe <- make.names(levels(surg))

  design_s <- model.matrix(~ 0 + surg)
  colnames(design_s) <- safe

  contrasts_s <- makeContrasts(
    Spine_vs_Laparoscopic    = Spine - Laparoscopic,
    Spine_vs_HeadNeck        = Spine - Head.Neck,
    Laparoscopic_vs_HeadNeck = Laparoscopic - Head.Neck,
    levels = safe
  )

  fit_s  <- lmFit(t(expr_mat), design_s)
  fit_s2 <- contrasts.fit(fit_s, contrasts_s)
  fit_s2 <- eBayes(fit_s2, trend = TRUE, robust = TRUE)

  omnibus <- topTable(fit_s2, number = Inf, sort.by = "F") |>
    rownames_to_column("AptName") |>
    as_tibble() |>
    left_join(analyte_info |> select(AptName, Target, EntrezGeneSymbol), by = "AptName") |>
    mutate(draw = draw_label)

  pairwise <- map(colnames(contrasts_s), function(ct) {
    topTable(fit_s2, coef = ct, number = Inf, sort.by = "P") |>
      rownames_to_column("AptName") |>
      as_tibble() |>
      left_join(analyte_info |> select(AptName, Target, EntrezGeneSymbol), by = "AptName") |>
      mutate(contrast = ct, draw = draw_label)
  }) |> bind_rows()

  list(omnibus = omnibus, pairwise = pairwise)
}

# ── PCA helper ─────────────────────────────────────────────────────────────────
save_pca_plot <- function(expr_mat, sample_df, color_col, palette, draw_label,
                          cov_label, out_path) {
  keep <- !is.na(sample_df[[color_col]])
  pca  <- prcomp(expr_mat[keep, ], center = TRUE, scale. = TRUE)
  var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

  df <- as_tibble(pca$x[, 1:2], rownames = "patient") |>
    mutate(group = sample_df[[color_col]][keep],
           patient = sample_df$patient[keep])

  p <- ggplot(df, aes(x = PC1, y = PC2, color = as.factor(group), label = patient)) +
    geom_point(size = 3.5, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
    scale_color_manual(values = palette, name = cov_label) +
    labs(title = paste0(draw_label, " — PCA colored by ", cov_label),
         x = paste0("PC1 (", var_exp[1], "% var)"),
         y = paste0("PC2 (", var_exp[2], "% var)")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")

  ggsave(out_path, p, width = 7, height = 6, dpi = 300)
}

# ── Heatmap helper ─────────────────────────────────────────────────────────────
save_covariate_heatmap <- function(results_df, expr_mat, sample_df,
                                   draw_label, out_path) {
  hit_apts <- results_df |> dplyr::filter(P.Value < 0.05) |> pull(AptName) |> unique()
  if (length(hit_apts) == 0) {
    message("No hits at p<0.05 for ", draw_label, " covariate heatmap — skipping.")
    return(invisible(NULL))
  }

  hm_mat <- results_df |>
    dplyr::filter(AptName %in% hit_apts) |>
    mutate(label       = if_else(is.na(Target) | Target == "", AptName, Target),
           neg_log10_p = pmin(-log10(P.Value), 4)) |>
    group_by(covariate, label) |>
    summarise(neg_log10_p = max(neg_log10_p), .groups = "drop") |>
    pivot_wider(names_from = covariate, values_from = neg_log10_p, values_fill = 0) |>
    column_to_rownames("label") |>
    as.matrix()

  pheatmap(
    hm_mat,
    color        = colorRampPalette(c("white", "#f7dc6f", "#c0392b"))(50),
    breaks       = seq(0, 4, length.out = 51),
    cluster_rows = TRUE, cluster_cols = FALSE,
    fontsize_row = 6, fontsize_col = 9,
    main         = paste0(draw_label, " — –log10(p) by protein and covariate (capped at 4)"),
    border_color = NA, angle_col = 45,
    filename     = out_path,
    width = 9, height = 14
  )
}

# ── Surgery heatmap helper ─────────────────────────────────────────────────────
save_surgery_heatmap <- function(omnibus_df, expr_mat, sample_df,
                                  draw_label, out_path) {
  sig_apts <- omnibus_df |> dplyr::filter(P.Value < 0.05) |> pull(AptName)
  if (length(sig_apts) == 0) {
    message("No omnibus hits at p<0.05 for ", draw_label, " surgery heatmap — skipping.")
    return(invisible(NULL))
  }

  hm_mat <- t(expr_mat[, sig_apts])
  rownames(hm_mat) <- coalesce(
    omnibus_df$Target[match(sig_apts, omnibus_df$AptName)], sig_apts
  )
  colnames(hm_mat) <- sample_df$patient

  annot_col <- data.frame(Surgery = sample_df$surgery_group,
                           row.names = sample_df$patient)

  pheatmap(
    hm_mat,
    annotation_col    = annot_col,
    annotation_colors = list(Surgery = col_surg),
    cluster_rows      = TRUE, cluster_cols  = TRUE,
    scale             = "row",
    color             = colorRampPalette(c("#2980b9", "white", "#c0392b"))(50),
    fontsize_row      = 7, fontsize_col = 8,
    main              = paste0(draw_label, " — ", length(sig_apts),
                               " surgery-omnibus proteins (p<0.05), row-scaled"),
    border_color      = NA, angle_col = 45,
    filename          = out_path,
    width = 10, height = max(6, length(sig_apts) * 0.15 + 3)
  )
  cat(draw_label, "surgery heatmap:", length(sig_apts), "proteins.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# Run for each draw
# ══════════════════════════════════════════════════════════════════════════════
draws <- list(
  Arterial = list(expr = expr_art, meta = soma_art),
  Venous   = list(expr = expr_ven, meta = soma_ven)
)

for (draw_label in names(draws)) {
  expr_d <- draws[[draw_label]]$expr
  meta_d <- draws[[draw_label]]$meta
  cat("\n══ Processing:", draw_label, "══\n")

  # ── Binary covariates ────────────────────────────────────────────────────────
  cov_results <- map(binary_covariates,
                     ~run_covariate_limma(expr_d, meta_d, .x, draw_label)) |>
    bind_rows()

  cat("Nominal hits (p<0.05) per covariate:\n")
  print(cov_results |> dplyr::filter(P.Value < 0.05) |> count(covariate))

  slug <- tolower(draw_label)
  write_csv(cov_results,
            file.path(table_dir, paste0("05_", slug, "_results.csv")))
  write_csv(
    cov_results |> dplyr::filter(P.Value < 0.05) |>
      group_by(covariate) |> slice_min(P.Value, n = 20) |> ungroup(),
    file.path(table_dir, paste0("05_", slug, "_top_hits.csv"))
  )

  save_covariate_heatmap(
    cov_results, expr_d, meta_d, draw_label,
    file.path(plot_dir, paste0("05_", slug, "_covariate_heatmap.pdf"))
  )

  # PCA per binary covariate
  bool_palette <- c("TRUE" = "#c0392b", "FALSE" = "#2980b9")
  for (cov in binary_covariates) {
    save_pca_plot(expr_d, meta_d, cov, bool_palette, draw_label, cov,
                  file.path(plot_dir, paste0("05_", slug, "_pca_", cov, ".png")))
  }

  # PCA by surgery type
  save_pca_plot(expr_d, meta_d, "surgery_group", col_surg, draw_label,
                "surgery_group",
                file.path(plot_dir, paste0("05_", slug, "_surgery_pca.png")))

  # ── Surgery type ─────────────────────────────────────────────────────────────
  surg_res <- run_surgery_limma(expr_d, meta_d, draw_label)

  cat("Surgery omnibus hits (F p<0.05):", sum(surg_res$omnibus$P.Value < 0.05), "\n")
  cat("Surgery pairwise hits (p<0.05) per contrast:\n")
  print(surg_res$pairwise |> dplyr::filter(P.Value < 0.05) |> count(contrast))

  write_csv(surg_res$omnibus,
            file.path(table_dir, paste0("05_", slug, "_surgery_omnibus.csv")))
  write_csv(
    surg_res$omnibus |> dplyr::filter(P.Value < 0.05) |> slice_min(P.Value, n = 20),
    file.path(table_dir, paste0("05_", slug, "_surgery_omnibus_top.csv"))
  )
  write_csv(surg_res$pairwise,
            file.path(table_dir, paste0("05_", slug, "_surgery_pairwise.csv")))
  write_csv(
    surg_res$pairwise |> dplyr::filter(P.Value < 0.05) |>
      group_by(contrast) |> slice_min(P.Value, n = 10) |> ungroup(),
    file.path(table_dir, paste0("05_", slug, "_surgery_pairwise_top.csv"))
  )

  save_surgery_heatmap(
    surg_res$omnibus, expr_d, meta_d, draw_label,
    file.path(plot_dir, paste0("05_", slug, "_surgery_heatmap.pdf"))
  )
}

cat("\nAll outputs saved to outputs/plots/ and outputs/tables/\n")
