# 09_opioid_analysis.R
# Focused analysis of endogenous opioid system (EOS) proteins — endorphins,
# enkephalins, dynorphins, nociceptin — in the perioperative A-V proteomics dataset.
#
# Reads: outputs/tables/03_*, 05_*, 06_* CSVs (from somascan pipeline)
#        SS-229807...anmlSMP.adat (for per-patient expression; requires SomaDataIO)
#
# Outputs:
#   outputs/tables/09_eos_*.csv
#   outputs/plots/09_eos_*.png / *.pdf

library(tidyverse)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

data_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path) |> dirname(),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)

out_tables <- file.path(data_dir, "outputs", "tables")
out_plots  <- file.path(data_dir, "outputs", "plots")

col_av   <- c("Arterial > Venous" = "#c0392b", "Venous > Arterial" = "#2980b9")
col_surg <- setNames(RColorBrewer::brewer.pal(3, "Set2"),
                     c("Head/Neck", "Laparoscopic", "Spine"))
col_yn   <- c("Yes" = "#c0392b", "No" = "#2980b9")

# ---------------------------------------------------------------------------
# 1. EOS gene universe
# ---------------------------------------------------------------------------

eos_genes <- tribble(
  ~gene_symbol, ~short_name,          ~eos_category,     ~eos_role,
  # Precursors
  "POMC",       "POMC",               "Precursor",        "Pro-opiomelanocortin → β-endorphin, ACTH, α-MSH",
  "PENK",       "Proenkephalin",      "Precursor",        "Proenkephalin → met-enkephalin, leu-enkephalin",
  "PDYN",       "Prodynorphin",       "Precursor",        "Prodynorphin → dynorphin A/B, α-neoendorphin",
  "PNOC",       "Pronociceptin",      "Precursor",        "Pronociceptin → nociceptin / orphanin FQ",
  # Receptors
  "OPRM1",      "MOR (μ)",            "Receptor",         "μ-opioid receptor: β-endorphin / enkephalin target",
  "OPRD1",      "DOR (δ)",            "Receptor",         "δ-opioid receptor: enkephalin-preferring",
  "OPRK1",      "KOR (κ)",            "Receptor",         "κ-opioid receptor: dynorphin target",
  "OPRL1",      "NOP",                "Receptor",         "Nociceptin/OFQ receptor",
  # Degrading enzymes
  "MME",        "Neprilysin",         "Degradation",      "Primary enkephalin endopeptidase; soluble form in plasma",
  "ANPEP",      "Aminopeptidase N",   "Degradation",      "N-terminal enkephalin cleavage (CD13)",
  "ACE",        "ACE",                "Degradation",      "Cleaves enkephalin C-terminus; also degrades bradykinin",
  "PREP",       "Prolyl endopeptidase","Degradation",     "Neuropeptide degradation including opioid fragments",
  # Processing enzymes (POMC → mature peptides)
  "PCSK1",      "PC1/3",              "Processing",       "POMC → ACTH / β-endorphin cleavage (pituitary)",
  "PCSK2",      "PC2",                "Processing",       "POMC → α-MSH / β-endorphin final processing",
  "CPE",        "Carboxypeptidase E", "Processing",       "Opioid peptide C-terminal maturation; circulates in plasma",
  # HPA axis / upstream drivers
  "CRH",        "CRH",                "HPA axis",         "Triggers ACTH/β-endorphin release from pituitary",
  "CGA",        "Glycoprotein α",     "HPA axis",         "Pituitary secretory marker (α-subunit common chain)",
  # ECS-EOS crosstalk
  "FAAH",       "FAAH",               "ECS-EOS crosstalk","Degrades AEA; ECS-opioid interaction node",
  "CNR1",       "CB1",                "ECS-EOS crosstalk","CB1/MOR heterodimer modulates opioid signaling",
  "BDNF",       "BDNF",               "ECS-EOS crosstalk","Upregulated by opioids and ECS; pain sensitization",
  "SIGMAR1",    "Sigma-1R",           "ECS-EOS crosstalk","Modulates both opioid and cannabinoid signaling"
)

cat("EOS gene universe:", nrow(eos_genes), "genes\n\n")

# ---------------------------------------------------------------------------
# 2. Load pre-computed SomaScan results
# ---------------------------------------------------------------------------

av_soma  <- read_csv(file.path(out_tables, "03_av_results_full.csv"),     show_col_types = FALSE)
cov_soma <- read_csv(file.path(out_tables, "05_covariate_results.csv"),   show_col_types = FALSE)
dlt_soma <- read_csv(file.path(out_tables, "06_delta_results.csv"),       show_col_types = FALSE)

# Join to EOS universe
eos_av  <- inner_join(av_soma,  eos_genes, by = c("EntrezGeneSymbol" = "gene_symbol"))
eos_cov <- inner_join(cov_soma, eos_genes, by = c("EntrezGeneSymbol" = "gene_symbol"))
eos_dlt <- inner_join(dlt_soma, eos_genes, by = c("EntrezGeneSymbol" = "gene_symbol"))

# ---------------------------------------------------------------------------
# 3. Tables
# ---------------------------------------------------------------------------

eos_av_out <- eos_av |>
  arrange(P.Value) |>
  select(AptName, Target, EntrezGeneSymbol, short_name, eos_category, eos_role,
         logFC, AveExpr, t, P.Value, adj.P.Val, direction)

eos_cov_out <- eos_cov |>
  arrange(P.Value) |>
  select(AptName, covariate, Target, EntrezGeneSymbol, short_name, eos_category, eos_role,
         logFC, t, P.Value, adj.P.Val)

eos_dlt_out <- eos_dlt |>
  arrange(P.Value) |>
  select(AptName, covariate, Target, EntrezGeneSymbol, short_name, eos_category, eos_role,
         logFC, t, P.Value, adj.P.Val)

write_csv(eos_av_out,  file.path(out_tables, "09_eos_av_full.csv"))
write_csv(eos_cov_out, file.path(out_tables, "09_eos_covariate_full.csv"))
write_csv(eos_dlt_out, file.path(out_tables, "09_eos_delta_full.csv"))

# Combined nominal summary
eos_summary <- bind_rows(
  eos_av_out  |> filter(P.Value < 0.05) |>
    transmute(analysis = "A-V gradient", covariate = NA_character_,
              EntrezGeneSymbol, short_name, eos_category, eos_role,
              logFC, p_value = P.Value, direction),
  eos_cov_out |> filter(P.Value < 0.05) |>
    transmute(analysis = "Covariate", covariate,
              EntrezGeneSymbol, short_name, eos_category, eos_role,
              logFC, p_value = P.Value,
              direction = if_else(logFC > 0, "Higher in Yes", "Higher in No")),
  eos_dlt_out |> filter(P.Value < 0.05) |>
    transmute(analysis = "A-V delta", covariate,
              EntrezGeneSymbol, short_name, eos_category, eos_role,
              logFC, p_value = P.Value,
              direction = if_else(logFC > 0, "Larger gradient in Yes", "Smaller gradient in Yes"))
) |>
  arrange(p_value)

write_csv(eos_summary, file.path(out_tables, "09_eos_summary_hits.csv"))

cat("Tables written.\n")
cat("  A-V full:        ", nrow(eos_av_out), "rows\n")
cat("  Covariate full:  ", nrow(eos_cov_out), "rows\n")
cat("  Delta full:      ", nrow(eos_dlt_out), "rows\n")
cat("  Summary nominal: ", nrow(eos_summary), "rows\n\n")

# ---------------------------------------------------------------------------
# 4. Plot A: Forest plot — A-V gradient for all detected EOS genes
# ---------------------------------------------------------------------------

forest_df <- eos_av_out |>
  # Where a gene has multiple aptamers, keep the one with smallest p-value
  group_by(EntrezGeneSymbol) |>
  slice_min(P.Value, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    label     = if_else(is.na(Target) | Target == "", EntrezGeneSymbol, Target),
    label     = fct_reorder(label, logFC),
    sig       = P.Value < 0.05,
    dot_color = if_else(logFC > 0, "Arterial > Venous", "Venous > Arterial"),
    se        = logFC / t,
    ci_lo     = logFC - 1.96 * se,
    ci_hi     = logFC + 1.96 * se
  )

p_forest <- ggplot(forest_df,
                   aes(x = logFC, y = label, color = dot_color, size = sig)) +
  geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3, linewidth = 0.5,
                 alpha = 0.6) +
  geom_point(alpha = 0.9) +
  geom_text(
    data = forest_df |> filter(sig),
    aes(label = sprintf("p = %.3f", P.Value)),
    hjust = -0.15, size = 2.6, show.legend = FALSE
  ) +
  scale_color_manual(values = col_av) +
  scale_size_manual(values = c("FALSE" = 2.5, "TRUE" = 4.5), guide = "none") +
  facet_grid(eos_category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title    = "Endogenous opioid system: arterial–venous gradient",
    subtitle = "log2FC (Arterial / Venous), paired limma, N = 12. Bars = ±1.96 SE. Larger points = p < 0.05.",
    x        = "log2FC (Arterial / Venous)",
    y        = NULL,
    color    = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text.y    = element_text(size = 8, angle = 0),
    strip.placement = "outside",
    panel.spacing   = unit(0.25, "lines")
  )

ggsave(file.path(out_plots, "09_eos_forest_av.png"),
       p_forest, width = 9, height = nrow(forest_df) * 0.42 + 3.5,
       limitsize = FALSE, dpi = 150)
cat("Saved: 09_eos_forest_av.png\n")

# ---------------------------------------------------------------------------
# 5. Plot B: Covariate heatmap — –log10(p) across EOS genes × covariates
# ---------------------------------------------------------------------------

cov_hm_df <- eos_cov_out |>
  group_by(EntrezGeneSymbol, covariate) |>
  slice_min(P.Value, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    label       = if_else(is.na(Target) | Target == "",
                          EntrezGeneSymbol,
                          paste0(Target, " (", EntrezGeneSymbol, ")")),
    neg_log10_p = pmin(-log10(P.Value), 4)
  ) |>
  select(label, covariate, neg_log10_p) |>
  pivot_wider(names_from = covariate, values_from = neg_log10_p, values_fill = 0) |>
  column_to_rownames("label") |>
  as.matrix()

# Only keep rows with at least one nominal hit
row_max <- apply(cov_hm_df, 1, max)
cov_hm_df <- cov_hm_df[row_max >= -log10(0.05), , drop = FALSE]

if (nrow(cov_hm_df) > 0) {
  colnames(cov_hm_df) <- recode(colnames(cov_hm_df),
    obese         = "Obese",
    hypothermia   = "Hypothermia",
    hyperglycemia = "Hyperglycemia",
    male          = "Male sex"
  )

  pheatmap(
    cov_hm_df,
    color        = colorRampPalette(c("white", "#f7dc6f", "#c0392b"))(50),
    breaks       = seq(0, 4, length.out = 51),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 9,
    fontsize_col = 10,
    main         = "EOS proteins × covariate: –log10(p)\n(venous samples; cap = 4; only genes with ≥1 nominal hit shown)",
    border_color = NA,
    angle_col    = 45,
    filename     = file.path(out_plots, "09_eos_covariate_heatmap.pdf"),
    width = 6, height = nrow(cov_hm_df) * 0.35 + 2.5
  )
  cat("Saved: 09_eos_covariate_heatmap.pdf\n")
}

# ---------------------------------------------------------------------------
# 6. Plot C: A-V delta heatmap — gradient magnitude × covariate
# ---------------------------------------------------------------------------

dlt_hm_df <- eos_dlt_out |>
  group_by(EntrezGeneSymbol, covariate) |>
  slice_min(P.Value, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    label       = if_else(is.na(Target) | Target == "",
                          EntrezGeneSymbol,
                          paste0(Target, " (", EntrezGeneSymbol, ")")),
    neg_log10_p = pmin(-log10(P.Value), 4)
  ) |>
  select(label, covariate, neg_log10_p) |>
  pivot_wider(names_from = covariate, values_from = neg_log10_p, values_fill = 0) |>
  column_to_rownames("label") |>
  as.matrix()

row_max_dlt <- apply(dlt_hm_df, 1, max)
dlt_hm_df <- dlt_hm_df[row_max_dlt >= -log10(0.05), , drop = FALSE]

if (nrow(dlt_hm_df) > 0) {
  colnames(dlt_hm_df) <- recode(colnames(dlt_hm_df),
    obese         = "Obese",
    hypothermia   = "Hypothermia",
    hyperglycemia = "Hyperglycemia",
    male          = "Male sex"
  )

  pheatmap(
    dlt_hm_df,
    color        = colorRampPalette(c("white", "#a8d8ea", "#2c3e50"))(50),
    breaks       = seq(0, 4, length.out = 51),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 9,
    fontsize_col = 10,
    main         = "EOS proteins × covariate: A-V delta –log10(p)\n(does gradient magnitude differ by subgroup?; cap = 4)",
    border_color = NA,
    angle_col    = 45,
    filename     = file.path(out_plots, "09_eos_delta_heatmap.pdf"),
    width = 6, height = nrow(dlt_hm_df) * 0.35 + 2.5
  )
  cat("Saved: 09_eos_delta_heatmap.pdf\n")
}

# ---------------------------------------------------------------------------
# 7. Per-patient dot plots (require SomaDataIO for raw expression)
# ---------------------------------------------------------------------------

has_somaio <- requireNamespace("SomaDataIO", quietly = TRUE)

if (!has_somaio) {
  cat("\nSomaDataIO not available — skipping per-patient paired dot plots.\n")
} else {
  library(SomaDataIO)

  adat_file <- file.path(
    data_dir,
    "SS-229807_v4.1_EDTA_Plasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat"
  )
  soma <- read_adat(adat_file)

  analyte_info <- getAnalyteInfo(soma)
  pass_analytes <- analyte_info |> filter(ColCheck == "PASS") |> pull(AptName)

  surgery_map <- tribble(
    ~patient, ~surgery_group,
    "0179", "Spine",
    "0193", "Laparoscopic",
    "0209", "Spine",
    "0215", "Spine",
    "0223", "Head/Neck",
    "0245", "Laparoscopic",
    "0250", "Laparoscopic",
    "0254", "Head/Neck",
    "0268", "Head/Neck",
    "0275", "Head/Neck",
    "0284", "Spine",
    "0291", "Spine"
  )

  soma_samples <- soma |>
    filter(SampleType == "Sample") |>
    mutate(
      patient = str_extract(SampleId, "\\d{4}"),
      draw    = if_else(str_ends(SampleId, "-2"), "Venous", "Arterial")
    ) |>
    left_join(surgery_map, by = "patient") |>
    arrange(patient, draw)

  covariate_map <- tribble(
    ~patient, ~obese, ~hypothermia, ~hyperglycemia, ~male,
    "0179", FALSE, TRUE,  FALSE, TRUE,
    "0193", TRUE,  FALSE, FALSE, FALSE,
    "0209", FALSE, TRUE,  NA,    FALSE,
    "0215", FALSE, TRUE,  FALSE, TRUE,
    "0223", FALSE, TRUE,  FALSE, TRUE,
    "0245", TRUE,  FALSE, TRUE,  FALSE,
    "0250", TRUE,  TRUE,  TRUE,  FALSE,
    "0254", FALSE, FALSE, NA,    TRUE,
    "0268", FALSE, FALSE, FALSE, FALSE,
    "0275", FALSE, TRUE,  FALSE, FALSE,
    "0284", FALSE, TRUE,  FALSE, FALSE,
    "0291", TRUE,  TRUE,  FALSE, FALSE
  )

  soma_samples <- soma_samples |> left_join(covariate_map, by = "patient")

  expr_log2 <- soma_samples |>
    select(all_of(pass_analytes)) |>
    as.matrix() |>
    log2()
  rownames(expr_log2) <- soma_samples$SampleId

  # ── Plot D: Paired A vs V for top EOS hits ──────────────────────────────────

  # Best aptamer per EOS gene by smallest A-V p-value
  top_av_apts <- eos_av_out |>
    filter(AptName %in% pass_analytes) |>
    group_by(EntrezGeneSymbol) |>
    slice_min(P.Value, n = 1, with_ties = FALSE) |>
    ungroup() |>
    slice_min(P.Value, n = 12) |>
    mutate(panel_label = paste0(
      if_else(is.na(Target) | Target == "", EntrezGeneSymbol, Target),
      "\n(p=", signif(P.Value, 2), ")"
    ))

  if (nrow(top_av_apts) > 0) {
    paired_long <- soma_samples |>
      select(patient, draw, surgery_group, all_of(top_av_apts$AptName)) |>
      pivot_longer(all_of(top_av_apts$AptName),
                   names_to = "AptName", values_to = "log2_rfu") |>
      left_join(top_av_apts |> select(AptName, panel_label), by = "AptName") |>
      mutate(draw = factor(draw, levels = c("Venous", "Arterial")))

    p_paired <- ggplot(paired_long,
                       aes(x = draw, y = log2_rfu, group = patient,
                           color = surgery_group)) +
      geom_line(alpha = 0.5, linewidth = 0.7) +
      geom_point(size = 2.2) +
      facet_wrap(~ panel_label, scales = "free_y", ncol = 4) +
      scale_color_manual(values = col_surg) +
      labs(
        title    = "EOS top hits: paired arterial vs. venous per patient",
        subtitle = "Lines connect the same patient across draw types.",
        x = NULL, y = "log2 RFU (ANML-normalized)", color = "Surgery"
      ) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom")

    ggsave(file.path(out_plots, "09_eos_paired_av.png"),
           p_paired,
           width = 12, height = ceiling(nrow(top_av_apts) / 4) * 3.2 + 2.5,
           limitsize = FALSE, dpi = 150)
    cat("Saved: 09_eos_paired_av.png\n")
  }

  # ── Plot E: Covariate dot plots for top EOS covariate hits ──────────────────

  top_cov_hits <- eos_cov_out |>
    filter(P.Value < 0.05, AptName %in% pass_analytes) |>
    group_by(covariate, EntrezGeneSymbol) |>
    slice_min(P.Value, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(panel_label = paste0(covariate, ": ",
                                if_else(is.na(Target) | Target == "",
                                        EntrezGeneSymbol, Target),
                                "\n(p=", signif(P.Value, 2), ")"))

  if (nrow(top_cov_hits) > 0) {
    ven_idx     <- soma_samples$draw == "Venous"
    soma_ven    <- soma_samples[ven_idx, ]

    cov_plot_data <- top_cov_hits |>
      select(covariate, AptName, panel_label) |>
      pmap_dfr(function(covariate, AptName, panel_label) {
        if (!AptName %in% colnames(expr_log2)) return(NULL)
        y   <- expr_log2[ven_idx, AptName]
        grp <- soma_ven[[covariate]]
        tibble(
          patient     = soma_ven$patient,
          log2_rfu    = y,
          group_val   = grp,
          panel_label = panel_label
        )
      }) |>
      filter(!is.na(group_val)) |>
      mutate(group_label = if_else(group_val, "Yes", "No"))

    if (nrow(cov_plot_data) > 0) {
      p_cov_dots <- ggplot(cov_plot_data,
                           aes(x = group_label, y = log2_rfu, color = group_label)) +
        geom_jitter(width = 0.12, size = 2.5, alpha = 0.85) +
        stat_summary(fun = median, geom = "crossbar", width = 0.4,
                     color = "black", linewidth = 0.5) +
        facet_wrap(~ panel_label, scales = "free_y",
                   ncol = min(4, length(unique(cov_plot_data$panel_label)))) +
        scale_color_manual(values = col_yn) +
        labs(
          title    = "EOS proteins: levels by covariate (venous samples)",
          subtitle = "Crossbar = median. p = nominal (N = 12; exploratory).",
          x = NULL, y = "log2 RFU", color = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(legend.position = "none")

      n_panels <- length(unique(cov_plot_data$panel_label))
      ggsave(file.path(out_plots, "09_eos_covariate_dots.png"),
             p_cov_dots,
             width  = min(4, n_panels) * 3.2,
             height = ceiling(n_panels / 4) * 3.5 + 2,
             limitsize = FALSE, dpi = 150)
      cat("Saved: 09_eos_covariate_dots.png\n")
    }
  }

  # ── Plot F: A-V delta dot plots for top EOS delta hits ──────────────────────

  # Rebuild delta matrix (log2 arterial − log2 venous per patient)
  art_soma <- soma_samples |> filter(draw == "Arterial") |> arrange(patient)
  ven_soma <- soma_samples |> filter(draw == "Venous")   |> arrange(patient)
  stopifnot(all(art_soma$patient == ven_soma$patient))

  expr_art <- expr_log2[art_soma$SampleId, ]
  expr_ven <- expr_log2[ven_soma$SampleId, ]
  delta_mat <- expr_art - expr_ven
  rownames(delta_mat) <- art_soma$patient

  delta_meta <- art_soma |>
    select(patient, surgery_group, obese, hypothermia, hyperglycemia, male)

  top_dlt_hits <- eos_dlt_out |>
    filter(P.Value < 0.05, AptName %in% pass_analytes) |>
    group_by(covariate, EntrezGeneSymbol) |>
    slice_min(P.Value, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(panel_label = paste0(covariate, ": ",
                                if_else(is.na(Target) | Target == "",
                                        EntrezGeneSymbol, Target),
                                "\n(p=", signif(P.Value, 2), ")"))

  if (nrow(top_dlt_hits) > 0) {
    dlt_plot_data <- top_dlt_hits |>
      select(covariate, AptName, panel_label) |>
      pmap_dfr(function(covariate, AptName, panel_label) {
        if (!AptName %in% colnames(delta_mat)) return(NULL)
        grp <- delta_meta[[covariate]]
        tibble(
          patient     = delta_meta$patient,
          delta       = delta_mat[, AptName],
          group_val   = grp,
          panel_label = panel_label
        )
      }) |>
      filter(!is.na(group_val)) |>
      mutate(group_label = if_else(group_val, "Yes", "No"))

    if (nrow(dlt_plot_data) > 0) {
      p_dlt_dots <- ggplot(dlt_plot_data,
                           aes(x = group_label, y = delta, color = group_label)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
        geom_jitter(width = 0.12, size = 2.5, alpha = 0.85) +
        stat_summary(fun = median, geom = "crossbar", width = 0.4,
                     color = "black", linewidth = 0.5) +
        facet_wrap(~ panel_label, scales = "free_y",
                   ncol = min(4, length(unique(dlt_plot_data$panel_label)))) +
        scale_color_manual(values = col_yn) +
        labs(
          title    = "EOS proteins: A-V gradient magnitude by covariate",
          subtitle = "y = log2(arterial) − log2(venous). Positive = more arterial > venous in that subgroup.",
          x = NULL, y = "log2 A-V delta", color = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(legend.position = "none",
              strip.text = element_text(size = 8))

      n_panels_dlt <- length(unique(dlt_plot_data$panel_label))
      ggsave(file.path(out_plots, "09_eos_delta_dots.png"),
             p_dlt_dots,
             width  = min(4, n_panels_dlt) * 3.2,
             height = ceiling(n_panels_dlt / 4) * 3.5 + 2.5,
             limitsize = FALSE, dpi = 150)
      cat("Saved: 09_eos_delta_dots.png\n")
    }
  }
}

# ---------------------------------------------------------------------------
# 8. Print final summary to console
# ---------------------------------------------------------------------------

cat("\n========================================================\n")
cat("ENDOGENOUS OPIOID SYSTEM — PERIOPERATIVE A-V ANALYSIS\n")
cat("========================================================\n\n")

cat("── A-V gradient (nominal hits, p < 0.05) ──────────────\n")
eos_summary |> filter(analysis == "A-V gradient") |>
  select(short_name, eos_category, logFC, p_value, direction) |>
  mutate(across(c(logFC, p_value), ~ signif(., 3))) |>
  print(n = Inf)

cat("\n── Covariate hits (nominal, venous samples) ───────────\n")
eos_summary |> filter(analysis == "Covariate") |>
  select(covariate, short_name, eos_category, logFC, p_value, direction) |>
  mutate(across(c(logFC, p_value), ~ signif(., 3))) |>
  arrange(covariate, p_value) |>
  print(n = Inf)

cat("\n── A-V delta hits (gradient modified by covariate) ────\n")
eos_summary |> filter(analysis == "A-V delta") |>
  select(covariate, short_name, eos_category, logFC, p_value, direction) |>
  mutate(across(c(logFC, p_value), ~ signif(., 3))) |>
  arrange(covariate, p_value) |>
  print(n = Inf)

cat("\nAll outputs in outputs/tables/09_eos_*.csv and outputs/plots/09_eos_*\n")
