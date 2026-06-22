# 05_covariate_analysis.R
# Limma covariate analysis on venous samples: obesity, hypothermia, hyperglycemia, sex.
# Outputs: plots/05_covariate_heatmap.pdf, 05_covariate_dotplots.png
#          tables/05_covariate_results.csv, 05_covariate_top_hits.csv

data_dir <- tryCatch(
  normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)
source(file.path(data_dir, "R_scripts", "00_setup.R"))

# ── Venous samples only ────────────────────────────────────────────────────────
ven_idx  <- soma_samples$draw == "Venous"
soma_ven <- soma_samples[ven_idx, ]
expr_ven <- expr_all[ven_idx, ]

# ── Limma per covariate ────────────────────────────────────────────────────────
run_soma_covariate <- function(expr_mat, sample_df, cov_name) {
  group <- sample_df[[cov_name]]
  keep  <- !is.na(group)
  if (sum(keep) < 6 || length(unique(group[keep])) < 2) {
    message("Skipping ", cov_name, ": insufficient observations")
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
    mutate(covariate = cov_name,
           n_pos = sum(group[keep] == TRUE,  na.rm = TRUE),
           n_neg = sum(group[keep] == FALSE, na.rm = TRUE))
}

cov_results_all <- map(binary_covariates, ~run_soma_covariate(expr_ven, soma_ven, .x)) |>
  bind_rows()

cat("Nominal hits (p<0.05) per covariate:\n")
print(cov_results_all |> dplyr::filter(P.Value < 0.05) |> count(covariate))

write_csv(cov_results_all, file.path(table_dir, "05_covariate_results.csv"))
write_csv(
  cov_results_all |> dplyr::filter(P.Value < 0.05) |>
    group_by(covariate) |> slice_min(P.Value, n = 20) |> ungroup(),
  file.path(table_dir, "05_covariate_top_hits.csv")
)

# ── Heatmap: –log10(p) across covariates ──────────────────────────────────────
hit_apts <- cov_results_all |> dplyr::filter(P.Value < 0.05) |> pull(AptName) |> unique()

if (length(hit_apts) > 0) {
  hm_cov <- cov_results_all |>
    dplyr::filter(AptName %in% hit_apts) |>
    mutate(label       = if_else(is.na(Target) | Target == "", AptName, Target),
           neg_log10_p = pmin(-log10(P.Value), 4)) |>
    group_by(covariate, label) |>
    summarise(neg_log10_p = max(neg_log10_p), .groups = "drop") |>
    pivot_wider(names_from = covariate, values_from = neg_log10_p, values_fill = 0) |>
    column_to_rownames("label") |>
    as.matrix()

  pheatmap(
    hm_cov,
    color        = colorRampPalette(c("white", "#f7dc6f", "#c0392b"))(50),
    breaks       = seq(0, 4, length.out = 51),
    cluster_rows = TRUE, cluster_cols = FALSE,
    fontsize_row = 6, fontsize_col = 9,
    main         = "–log10(p) by protein and covariate (venous; capped at 4)",
    border_color = NA, angle_col = 45,
    filename     = file.path(plot_dir, "05_covariate_heatmap.pdf"),
    width = 9, height = 14
  )
}

# ── Dot plots: top 3 hits per covariate ───────────────────────────────────────
top3 <- cov_results_all |>
  dplyr::filter(P.Value < 0.05, covariate %in% binary_covariates) |>
  group_by(covariate) |> slice_min(P.Value, n = 3) |> ungroup()

if (nrow(top3) > 0) {
  cov_long <- soma_ven |>
    select(patient, all_of(binary_covariates)) |>
    pivot_longer(all_of(binary_covariates), names_to = "covariate", values_to = "group_val")
  apt_long <- as_tibble(expr_ven, rownames = "SampleId") |>
    mutate(patient = soma_ven$patient) |>
    select(patient, all_of(top3$AptName |> unique())) |>
    pivot_longer(-patient, names_to = "AptName", values_to = "log2_rfu")

  plot_data <- top3 |>
    select(covariate, AptName, Target) |>
    left_join(apt_long, by = "AptName") |>
    left_join(cov_long, by = c("patient", "covariate")) |>
    dplyr::filter(!is.na(group_val)) |>
    mutate(group_label = if_else(group_val, "Yes", "No"),
           protein     = if_else(is.na(Target) | Target == "", AptName, Target),
           panel_label = paste0(covariate, ": ", protein))

  p_dots <- ggplot(plot_data, aes(x = group_label, y = log2_rfu, color = group_label)) +
    geom_jitter(width = 0.1, size = 2.5, alpha = 0.8) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4,
                 color = "black", linewidth = 0.5) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = 3) +
    scale_color_manual(values = col_yn) +
    labs(title = "Top 3 proteins per covariate (venous log2 RFU)",
         x = NULL, y = "log2 RFU", color = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", strip.text = element_text(size = 8))

  ggsave(file.path(plot_dir, "05_covariate_dotplots.png"), p_dots,
         width = 12, height = 10, dpi = 300)
}

cat("\nOutputs saved to outputs/plots/ and outputs/tables/\n")
