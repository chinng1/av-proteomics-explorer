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
