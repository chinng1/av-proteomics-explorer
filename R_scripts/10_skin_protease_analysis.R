# 10_skin_protease_analysis.R
# Question: Do skin-expressed serine proteases (KLK family) and their epidermal
# inhibitors (SPINK family) show a systematic arteriovenous gradient during surgery?
#
# Motivation: The top two individual protein hits in the A-V analysis are
# SPINK9 (rank 1, p = 1.74e-4) and KLK7 (rank 3, p = 1.86e-3), both
# skin-specific. A cluster of related family members (SPINK1, KLK8, KLK10,
# SERPINA10) also nominally significant. These govern epidermal desquamation
# and their venous enrichment suggests surgical skin incision releases them
# into the venous return.
#
# Reads pre-computed output CSVs for statistics; sources 00_setup.R for
# per-patient expression values (requires SomaDataIO).
#
# Outputs: outputs/tables/10_*.csv, outputs/plots/10_*.png

data_dir <- tryCatch(
  normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)
source(file.path(data_dir, "R_scripts", "00_setup.R"))

out_tables <- file.path(data_dir, "outputs", "tables")
out_plots  <- file.path(data_dir, "outputs", "plots")

col_av_dir <- c("Arterial > Venous" = "#c0392b", "Venous > Arterial" = "#2980b9")

# ---------------------------------------------------------------------------
# 1. Define skin serine protease / inhibitor gene universe
# ---------------------------------------------------------------------------

skin_gene_table <- tribble(
  ~gene_symbol,  ~family,   ~skin_role,
  # SPINK — epidermal serine protease inhibitors (Kazal-type)
  "SPINK1",   "SPINK",  "Serine protease inhibitor; skin, pancreas",
  "SPINK2",   "SPINK",  "Serine protease inhibitor; skin, testis",
  "SPINK4",   "SPINK",  "Gastrointestinal and epidermal inhibitor",
  "SPINK5",   "SPINK",  "LEKTI: major epidermal barrier inhibitor; KLK5/7 substrate",
  "SPINK6",   "SPINK",  "Epidermal inhibitor of KLK5/14",
  "SPINK7",   "SPINK",  "Epidermis-specific inhibitor",
  "SPINK8",   "SPINK",  "Epidermis-associated inhibitor",
  "SPINK9",   "SPINK",  "Skin-specific inhibitor of KLK1; top A-V hit (rank 1)",
  "SPINK13",  "SPINK",  "SPINK family member; skin-associated",
  "SPINK14",  "SPINK",  "SPINK family member; skin-associated",
  # KLK — secreted serine proteases expressed in skin epidermis
  "KLK1",     "KLK",   "Tissue kallikrein; skin, kidney, pancreas",
  "KLK5",     "KLK",   "Skin-specific; cleaves SPINK5/LEKTI; drives desquamation cascade",
  "KLK6",     "KLK",   "Neurosin; brain and skin; cleaves desmoglein",
  "KLK7",     "KLK",   "Skin-specific chymase; cleaves corneodesmosin; top A-V hit (rank 3)",
  "KLK8",     "KLK",   "Neuropsin; skin cornified envelope; A-V hit",
  "KLK10",    "KLK",   "Kallikrein-10; epidermal and epithelial; A-V hit",
  "KLK11",    "KLK",   "Trypsin-like; expressed in skin",
  "KLK12",    "KLK",   "Kallikrein-12; skin and salivary gland",
  "KLK13",    "KLK",   "Kallikrein-13; skin and seminal plasma",
  "KLK14",    "KLK",   "Skin-specific; activates pro-KLK5/7; desquamation cascade",
  # Serpins with epidermal protease-inhibitor activity
  "SERPINA10","SERPIN", "Protein Z-dependent protease inhibitor; A-V hit",
  "SERPINB3", "SERPIN", "SCCA1: squamous cell carcinoma antigen; inhibits KLK and chymase",
  "SERPINB4", "SERPIN", "SCCA2: skin; inhibits KLK and chymase"
)

cat("Skin protease gene universe:", nrow(skin_gene_table), "entries\n")

# ---------------------------------------------------------------------------
# 2. A-V results (from 03_av_limma.R output)
# ---------------------------------------------------------------------------

av_soma <- read_csv(file.path(out_tables, "03_av_results_full.csv"),
                    show_col_types = FALSE)

skin_av <- av_soma |>
  inner_join(skin_gene_table, by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  arrange(P.Value) |>
  select(Target, EntrezGeneSymbol, family, skin_role,
         logFC, AveExpr, t, P.Value, adj.P.Val, direction)

cat("\n--- SomaScan A-V: skin protease genes ---\n")
cat("Detected in SomaScan:", nrow(skin_av), "\n")
cat("Nominal p < 0.05:", sum(skin_av$P.Value < 0.05, na.rm = TRUE), "\n\n")
print(skin_av, n = Inf)

write_csv(skin_av, file.path(out_tables, "10_skin_av.csv"))

# ---------------------------------------------------------------------------
# 3. Covariate effects (from 05_covariate_analysis.R output)
# ---------------------------------------------------------------------------

cov_soma <- read_csv(file.path(out_tables, "05_covariate_results.csv"),
                     show_col_types = FALSE)

skin_cov <- cov_soma |>
  inner_join(skin_gene_table, by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  arrange(P.Value) |>
  select(covariate, Target, EntrezGeneSymbol, family, skin_role,
         logFC, t, P.Value, adj.P.Val)

cat("\n--- Covariate effects: skin protease genes ---\n")
cat("Nominal p < 0.05:", sum(skin_cov$P.Value < 0.05, na.rm = TRUE), "\n\n")
skin_cov |> dplyr::filter(P.Value < 0.05) |> print(n = Inf)

write_csv(skin_cov, file.path(out_tables, "10_skin_covariate.csv"))

# ---------------------------------------------------------------------------
# 4. A-V delta × covariate (from 06_av_delta.R output)
# ---------------------------------------------------------------------------

delta_soma <- read_csv(file.path(out_tables, "06_delta_results.csv"),
                       show_col_types = FALSE)

skin_delta <- delta_soma |>
  inner_join(skin_gene_table, by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  arrange(P.Value) |>
  select(covariate, Target, EntrezGeneSymbol, family, skin_role,
         logFC, t, P.Value, adj.P.Val)

cat("\n--- A-V delta × covariate: skin protease genes ---\n")
cat("Nominal p < 0.05:", sum(skin_delta$P.Value < 0.05, na.rm = TRUE), "\n\n")
skin_delta |> dplyr::filter(P.Value < 0.05) |> print(n = Inf)

write_csv(skin_delta, file.path(out_tables, "10_skin_delta.csv"))

# ---------------------------------------------------------------------------
# 5. Surgery-type A-V gradient
# ---------------------------------------------------------------------------
# Does the skin-protease A-V gradient differ across surgery groups?
# Hypothesis: Head/Neck (more skin surface) > Laparoscopic (minimal skin) ≈ Spine.
# Test: Kruskal-Wallis on per-patient A-V delta, stratified by surgery group.
# Groups: Head/Neck n=4, Laparoscopic n=3, Spine n=5 — very small; exploratory only.

art_meta <- soma_samples |> dplyr::filter(draw == "Arterial") |> arrange(patient)
ven_meta <- soma_samples |> dplyr::filter(draw == "Venous")   |> arrange(patient)
stopifnot(all(art_meta$patient == ven_meta$patient))

skin_nom_apts <- av_soma |>
  inner_join(skin_gene_table |> dplyr::filter(family %in% c("SPINK", "KLK")),
             by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  dplyr::filter(P.Value < 0.05) |>
  group_by(EntrezGeneSymbol) |>
  slice_min(P.Value, n = 1) |>
  ungroup() |>
  arrange(P.Value)

if (nrow(skin_nom_apts) > 0) {
  nom_apts <- skin_nom_apts$AptName

  delta_skin <- expr_all[art_meta$SampleId, nom_apts, drop = FALSE] -
                expr_all[ven_meta$SampleId, nom_apts, drop = FALSE]
  rownames(delta_skin) <- art_meta$patient

  surg_delta_long <- as_tibble(delta_skin, rownames = "patient") |>
    pivot_longer(-patient, names_to = "AptName", values_to = "av_delta") |>
    left_join(art_meta |> select(patient, surgery_group), by = "patient") |>
    left_join(skin_nom_apts |> select(AptName, EntrezGeneSymbol, P.Value, logFC),
              by = "AptName") |>
    rename(gene_symbol = EntrezGeneSymbol)

  kw_results <- surg_delta_long |>
    group_by(gene_symbol) |>
    summarise(
      av_p          = first(P.Value),
      kw_p          = tryCatch(kruskal.test(av_delta ~ surgery_group)$p.value,
                               error = function(e) NA_real_),
      mean_HN       = mean(av_delta[surgery_group == "Head/Neck"],    na.rm = TRUE),
      mean_Lap      = mean(av_delta[surgery_group == "Laparoscopic"], na.rm = TRUE),
      mean_Spine    = mean(av_delta[surgery_group == "Spine"],        na.rm = TRUE),
      .groups       = "drop"
    ) |>
    arrange(kw_p)

  cat("\n--- Surgery-type K-W test: skin protease A-V delta (exploratory, N=3-5/group) ---\n")
  print(kw_results, n = Inf)

  write_csv(kw_results, file.path(out_tables, "10_skin_surgery_kw.csv"))
} else {
  surg_delta_long <- tibble()
  kw_results      <- tibble()
  cat("No nominally significant SPINK/KLK aptamers for surgery comparison.\n")
}

# ---------------------------------------------------------------------------
# 6. Plots
# ---------------------------------------------------------------------------

# ── Plot A: Forest plot — all detected family members ────────────────────────
if (nrow(skin_av) > 0) {
  p_forest <- skin_av |>
    mutate(
      label = paste0(EntrezGeneSymbol,
                     if_else(!is.na(Target) & Target != "" & Target != EntrezGeneSymbol,
                             paste0(" (", Target, ")"), "")),
      label     = fct_reorder(label, logFC),
      sig       = P.Value < 0.05,
      point_col = if_else(logFC > 0, "Arterial > Venous", "Venous > Arterial")
    ) |>
    ggplot(aes(x = logFC, y = label, color = point_col, size = sig)) +
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = if_else(sig, sprintf("p=%.3f", P.Value), "")),
              hjust = -0.15, size = 2.8, show.legend = FALSE) +
    scale_color_manual(values = col_av_dir) +
    scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4), guide = "none") +
    facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(
      title    = "Skin serine proteases & inhibitors: SomaScan A-V gradient (N = 12 paired)",
      subtitle = "Negative logFC = higher in venous blood (tissue releasing). Larger points = p < 0.05.",
      x        = "log2FC (Arterial / Venous, paired limma)",
      y        = NULL, color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text.y    = element_text(size = 9, angle = 0),
      strip.placement = "outside",
      panel.spacing   = unit(0.3, "lines")
    )

  ggsave(file.path(out_plots, "10_skin_forest.png"),
         p_forest, width = 9, height = 0.5 * nrow(skin_av) + 3,
         limitsize = FALSE, dpi = 150)
  cat("Saved: 10_skin_forest.png\n")
}

# ── Plot B: Paired trajectories — nominally significant SPINK/KLK members ───
if (nrow(skin_nom_apts) > 0) {
  nom_labels <- setNames(
    paste0(skin_nom_apts$EntrezGeneSymbol,
           "\n(logFC = ", round(skin_nom_apts$logFC, 2),
           ", p = ", signif(skin_nom_apts$P.Value, 2), ")"),
    skin_nom_apts$AptName
  )

  p_traj <- soma_samples |>
    select(patient, draw, surgery_group, all_of(nom_apts)) |>
    mutate(across(all_of(nom_apts), log2)) |>
    pivot_longer(all_of(nom_apts), names_to = "AptName", values_to = "log2_rfu") |>
    mutate(
      draw        = factor(draw, levels = c("Venous", "Arterial")),
      panel_label = factor(nom_labels[AptName], levels = nom_labels)
    ) |>
    ggplot(aes(x = draw, y = log2_rfu, group = patient, color = surgery_group)) +
    geom_line(alpha = 0.55, linewidth = 0.7) +
    geom_point(size = 2.5) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = min(3L, nrow(skin_nom_apts))) +
    scale_color_manual(values = col_surg) +
    labs(
      title    = "Skin serine proteases: per-patient A-V trajectories (p < 0.05)",
      subtitle = "Consistent downward slope = venous > arterial (tissue releasing into venous return).",
      x        = NULL, y = "log2 RFU", color = "Surgery"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", strip.text = element_text(size = 9))

  n_cols <- min(3L, nrow(skin_nom_apts))
  n_rows <- ceiling(nrow(skin_nom_apts) / n_cols)
  ggsave(file.path(out_plots, "10_skin_paired.png"),
         p_traj, width = 4 * n_cols, height = 3.5 * n_rows + 1.5,
         limitsize = FALSE, dpi = 200)
  cat("Saved: 10_skin_paired.png\n")
}

# ── Plot C: Surgery-type A-V delta for nominally significant hits ────────────
if (nrow(surg_delta_long) > 0 && nrow(kw_results) > 0) {
  kw_lookup <- setNames(kw_results$kw_p, kw_results$gene_symbol)
  av_lookup <- setNames(kw_results$av_p,  kw_results$gene_symbol)

  p_surg <- surg_delta_long |>
    mutate(
      panel_label = paste0(gene_symbol,
                           "\n(KW p = ", signif(kw_lookup[gene_symbol], 2), ")"),
      panel_label = fct_reorder(panel_label, av_lookup[gene_symbol])
    ) |>
    ggplot(aes(x = surgery_group, y = av_delta, color = surgery_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_jitter(width = 0.15, size = 3, alpha = 0.8) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4,
                 color = "black", linewidth = 0.5) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = min(3L, nrow(skin_nom_apts))) +
    scale_color_manual(values = col_surg) +
    labs(
      title    = "Skin protease A-V gradient by surgery type (exploratory)",
      subtitle = "Negative delta = venous > arterial. Black bar = median. KW = Kruskal-Wallis p across 3 groups.",
      x        = NULL, y = "log2 A-V delta (arterial − venous)", color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "none",
          strip.text       = element_text(size = 9),
          axis.text.x      = element_text(angle = 20, hjust = 1))

  n_cols <- min(3L, nrow(skin_nom_apts))
  n_rows <- ceiling(nrow(skin_nom_apts) / n_cols)
  ggsave(file.path(out_plots, "10_skin_surgery.png"),
         p_surg, width = 4 * n_cols, height = 3.5 * n_rows + 1.5,
         limitsize = FALSE, dpi = 200)
  cat("Saved: 10_skin_surgery.png\n")
}

# ---------------------------------------------------------------------------
# 7. Summary
# ---------------------------------------------------------------------------

skin_summary <- bind_rows(
  skin_av |>
    dplyr::filter(P.Value < 0.05) |>
    transmute(analysis = "A-V gradient", family,
              gene_symbol = EntrezGeneSymbol, target = Target,
              logFC = round(logFC, 3), p_value = signif(P.Value, 3), direction),
  skin_cov |>
    dplyr::filter(P.Value < 0.05) |>
    transmute(analysis = paste0("Covariate: ", covariate), family,
              gene_symbol = EntrezGeneSymbol, target = Target,
              logFC = round(logFC, 3), p_value = signif(P.Value, 3),
              direction = if_else(logFC > 0, "Up in Yes", "Up in No")),
  skin_delta |>
    dplyr::filter(P.Value < 0.05) |>
    transmute(analysis = paste0("Delta: ", covariate), family,
              gene_symbol = EntrezGeneSymbol, target = Target,
              logFC = round(logFC, 3), p_value = signif(P.Value, 3),
              direction = if_else(logFC > 0, "Larger in Yes", "Smaller in Yes"))
) |>
  arrange(p_value)

cat("\n========================================\n")
cat("SUMMARY: Skin protease nominal hits (p < 0.05)\n")
cat("========================================\n")
print(skin_summary, n = Inf)

write_csv(skin_summary, file.path(out_tables, "10_skin_summary.csv"))
cat("\nDone. Outputs in outputs/tables/10_*.csv and outputs/plots/10_*.png\n")
