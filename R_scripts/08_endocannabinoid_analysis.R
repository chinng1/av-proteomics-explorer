# 08_endocannabinoid_analysis.R
# Question: Are any proteins in the SomaScan A-V or covariate analyses related
# to endocannabinoid physiology or chronic exogenous cannabinoid exposure?
#
# Reads pre-computed output CSVs — does not require SomaDataIO.
# Luminex cytokine results are re-derived from AP_65plex_byJMR.csv.
#
# Outputs: outputs/tables/08_*.csv, outputs/plots/08_*.png

library(tidyverse)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(readr)

data_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path) |> dirname(),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)

out_tables <- file.path(data_dir, "outputs", "tables")
out_plots  <- file.path(data_dir, "outputs", "plots")

# ---------------------------------------------------------------------------
# 1. Define ECS gene universe
# ---------------------------------------------------------------------------

# Grouped by functional category for later annotation
ecs_gene_table <- tribble(
  ~gene_symbol,   ~category,                        ~ecs_role,
  # Core receptors
  "CNR1",         "Receptor",                       "CB1 receptor",
  "CNR2",         "Receptor",                       "CB2 receptor",
  "GPR55",        "Receptor",                       "Putative CB3 / LPI receptor",
  "GPR119",       "Receptor",                       "OEA / 2-OG receptor",
  "GPR18",        "Receptor",                       "NAGly receptor",
  "TRPV1",        "Receptor",                       "Anandamide / capsaicin receptor",
  # Synthesis
  "NAPEPLD",      "Synthesis",                      "NAPE-PLD: AEA synthesis",
  "DAGLA",        "Synthesis",                      "DAGLα: 2-AG synthesis",
  "DAGLB",        "Synthesis",                      "DAGLβ: 2-AG synthesis",
  "GDE1",         "Synthesis",                      "Lysophospholipase: AEA precursor",
  "PTPN22",       "Synthesis",                      "Regulates 2-AG synthesis (indirect)",
  # Degradation
  "FAAH",         "Degradation",                    "FAAH: AEA hydrolysis (CBD inhibits)",
  "FAAH2",        "Degradation",                    "FAAH2: alternate AEA hydrolysis",
  "MGLL",         "Degradation",                    "MAGL: 2-AG hydrolysis",
  "ABHD6",        "Degradation",                    "ABHD6: 2-AG hydrolysis",
  "ABHD12",       "Degradation",                    "ABHD12: 2-AG hydrolysis",
  # Reuptake / transport
  "FABP1",        "Transport",                      "Intracellular AEA transport",
  "FABP3",        "Transport",                      "Intracellular AEA transport",
  "FABP5",        "Transport",                      "Intracellular AEA transport",
  "FABP7",        "Transport",                      "Intracellular AEA transport",
  # Nuclear receptors / downstream
  "PPARG",        "Nuclear receptor",               "PPARγ: AEA / OEA / PEA target; CBD agonist",
  "PPARA",        "Nuclear receptor",               "PPARα: OEA target; anti-inflammatory",
  "PPARD",        "Nuclear receptor",               "PPARδ: endocannabinoid target",
  # Prostanoid / oxylipin crosstalk
  "PTGS1",        "Oxylipin crosstalk",             "COX-1: oxygenates AEA → prostamides",
  "PTGS2",        "Oxylipin crosstalk",             "COX-2: oxygenates AEA; suppressed by cannabinoids",
  "ALOX5",        "Oxylipin crosstalk",             "5-LOX: oxygenates AEA → leukotrienes",
  "ALOX12",       "Oxylipin crosstalk",             "12-LOX: oxygenates AEA",
  # Neurotrophins (ECS-regulated; relevant to developmental / periop neuroscience)
  "BDNF",         "Neurotrophin",                   "ECS retrograde signal at TrkB synapses",
  "NGF",          "Neurotrophin",                   "NGF drives TRPV1 sensitization; ECS feedback",
  "NTRK2",        "Neurotrophin",                   "TrkB: BDNF receptor; ECS-modulated",
  "NTRK1",        "Neurotrophin",                   "TrkA: NGF receptor",
  # Adipokines / metabolic (ECS tonically upregulated in obesity)
  "ADIPOQ",       "Metabolic",                      "Adiponectin; inverse correlation with ECS tone",
  "LEP",          "Metabolic",                      "Leptin; CB1 upregulation in obesity",
  "RETN",         "Metabolic",                      "Resistin; ECS-metabolic axis",
  # Inflammatory cytokines with strong ECS modulation
  "TNF",          "ECS-modulated cytokine",         "Suppressed by CB1/CB2 activation",
  "IL1B",         "ECS-modulated cytokine",         "Suppressed by CB2; NLRP3 inhibited by CBD",
  "IL6",          "ECS-modulated cytokine",         "Suppressed by CB2 agonism; CBD reduces IL-6",
  "IL10",         "ECS-modulated cytokine",         "Upregulated by cannabinoids (anti-inflam.)",
  "CXCL10",       "ECS-modulated cytokine",         "Suppressed by CB2 (IFN-γ–driven); top Luminex hit",
  "CCL2",         "ECS-modulated cytokine",         "MCP-1; CB2 suppresses monocyte recruitment",
  "MIF",          "ECS-modulated cytokine",         "Macrophage MIF suppressed by cannabinoids",
  "CX3CL1",       "ECS-modulated cytokine",         "Fractalkine; ECS modulates microglial CX3CR1",
  "FGF2",         "ECS-modulated cytokine",         "ECS / FGF2 overlap in neural development",
  "VEGFA",        "ECS-modulated cytokine",         "Angiogenesis; CB1 expressed on endothelium",
  # Serotonin axis (CBD acts at 5-HT1A)
  "HTR1A",        "Serotonin",                      "5-HT1A: direct CBD target",
  "SLC6A4",       "Serotonin",                      "SERT: serotonin transporter; CBD interaction",
  "TPH1",         "Serotonin",                      "Tryptophan hydroxylase",
  # TRP channels
  "TRPA1",        "TRP channel",                    "Activated by CBD; pain / inflammation",
  "TRPM8",        "TRP channel",                    "Inhibited by CBD; cold / pain",
  # Misc
  "SIGMAR1",      "Misc",                           "Sigma-1 receptor: CBD target; ER stress"
)

cat("ECS gene universe:", nrow(ecs_gene_table), "entries\n")

# ---------------------------------------------------------------------------
# 2. SomaScan A-V results
# ---------------------------------------------------------------------------

av_soma <- read_csv(file.path(out_tables, "03_av_results_full.csv"),
                    show_col_types = FALSE)

ecs_soma_av <- av_soma |>
  inner_join(ecs_gene_table, by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  arrange(P.Value) |>
  select(Target, EntrezGeneSymbol, category, ecs_role,
         logFC, AveExpr, t, P.Value, adj.P.Val, direction)

cat("\n--- SomaScan A-V: ECS genes found ---\n")
cat("Total ECS genes detected in SomaScan:", nrow(ecs_soma_av), "\n")
cat("Nominal p < 0.05:", sum(ecs_soma_av$P.Value < 0.05, na.rm = TRUE), "\n\n")
print(ecs_soma_av, n = Inf)

write_csv(ecs_soma_av, file.path(out_tables, "08_ecs_somascan_av.csv"))

# ---------------------------------------------------------------------------
# 3. SomaScan covariate results
# ---------------------------------------------------------------------------

cov_soma <- read_csv(file.path(out_tables, "05_covariate_results.csv"),
                     show_col_types = FALSE)

ecs_soma_cov <- cov_soma |>
  inner_join(ecs_gene_table, by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  arrange(P.Value) |>
  select(covariate, Target, EntrezGeneSymbol, category, ecs_role,
         logFC, t, P.Value, adj.P.Val)

cat("\n--- SomaScan covariate: ECS genes found ---\n")
cat("Nominal p < 0.05:", sum(ecs_soma_cov$P.Value < 0.05, na.rm = TRUE), "\n\n")
ecs_soma_cov |>
  filter(P.Value < 0.05) |>
  print(n = Inf)

write_csv(ecs_soma_cov, file.path(out_tables, "08_ecs_somascan_covariate.csv"))

# ---------------------------------------------------------------------------
# 4. SomaScan A-V delta (gradient × covariate)
# ---------------------------------------------------------------------------

delta_soma <- read_csv(file.path(out_tables, "06_delta_results.csv"),
                       show_col_types = FALSE)

ecs_soma_delta <- delta_soma |>
  inner_join(ecs_gene_table, by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  arrange(P.Value) |>
  select(covariate, Target, EntrezGeneSymbol, category, ecs_role,
         logFC, t, P.Value, adj.P.Val)

cat("\n--- SomaScan A-V delta × covariate: ECS genes ---\n")
cat("Nominal p < 0.05:", sum(ecs_soma_delta$P.Value < 0.05, na.rm = TRUE), "\n\n")
ecs_soma_delta |>
  filter(P.Value < 0.05) |>
  print(n = Inf)

write_csv(ecs_soma_delta, file.path(out_tables, "08_ecs_somascan_delta.csv"))

# ---------------------------------------------------------------------------
# 5. Luminex: ECS-relevant cytokines from the 65-plex
# ---------------------------------------------------------------------------

plate_map <- tribble(
  ~unknown_id,  ~patient, ~draw,
  "Unknown1",   "0179",   "Venous",
  "Unknown2",   "0193",   "Venous",
  "Unknown3",   "0209",   "Venous",
  "Unknown4",   "0215",   "Venous",
  "Unknown5",   "0179",   "Arterial",
  "Unknown6",   "0193",   "Arterial",
  "Unknown7",   "0209",   "Arterial",
  "Unknown8",   "0215",   "Arterial",
  "Unknown9",   "0223",   "Venous",
  "Unknown10",  "0245",   "Venous",
  "Unknown11",  "0250",   "Venous",
  "Unknown12",  "0254",   "Venous",
  "Unknown13",  "0223",   "Arterial",
  "Unknown14",  "0245",   "Arterial",
  "Unknown15",  "0250",   "Arterial",
  "Unknown16",  "0254",   "Arterial",
  "Unknown17",  "0268",   "Venous",
  "Unknown18",  "0275",   "Venous",
  "Unknown19",  "0284",   "Venous",
  "Unknown20",  "0291",   "Venous",
  "Unknown21",  "0268",   "Arterial",
  "Unknown22",  "0275",   "Arterial",
  "Unknown23",  "0284",   "Arterial",
  "Unknown24",  "0291",   "Arterial"
)

# Luminex name → gene symbol for ECS-relevant analytes on the 65-plex
lx_ecs_map <- tribble(
  ~luminex_name,         ~gene_symbol,
  "TNFa",                "TNF",
  "IL-1b",               "IL1B",
  "IL-6",                "IL6",
  "IL-10",               "IL10",
  "IP-10(CXCL10)",       "CXCL10",
  "MCP-1(CCL2)",         "CCL2",
  "MIF",                 "MIF",
  "Fractalkine(CX3CL1)", "CX3CL1",
  "FGF-2",               "FGF2",
  "VEGF-A",              "VEGFA",
  "NGFb",                "NGF",
  "HGF",                 "HGF"
)

lx_raw <- read_csv(
  file.path(data_dir, "AP_65plex_byJMR.csv"),
  skip = 54, show_col_types = FALSE
) |>
  filter(str_starts(Sample, "Unknown"))

lx_cols <- lx_raw |> select(-Location, -Sample, -`Total Events`) |> names()

lx_avg <- lx_raw |>
  mutate(across(all_of(lx_cols), ~ suppressWarnings(as.numeric(.)))) |>
  left_join(plate_map, by = c("Sample" = "unknown_id")) |>
  group_by(patient, draw) |>
  summarise(across(all_of(lx_cols), ~ mean(., na.rm = TRUE)), .groups = "drop")

# Restrict to ECS cytokines present in this run
ecs_lx_names <- lx_ecs_map |> filter(luminex_name %in% lx_cols) |> pull(luminex_name)

art_lx <- lx_avg |> filter(draw == "Arterial") |> arrange(patient)
ven_lx <- lx_avg |> filter(draw == "Venous")   |> arrange(patient)

lx_ecs_results <- map_dfr(ecs_lx_names, function(cyt) {
  a <- art_lx[[cyt]]; v <- ven_lx[[cyt]]
  test <- tryCatch(wilcox.test(a, v, paired = TRUE, exact = FALSE),
                   error = function(e) list(p.value = NA))
  tibble(
    luminex_name    = cyt,
    mean_arterial   = mean(a, na.rm = TRUE),
    mean_venous     = mean(v, na.rm = TRUE),
    log2FC          = log2((mean(a, na.rm = TRUE) + 0.1) / (mean(v, na.rm = TRUE) + 0.1)),
    p_value         = test$p.value,
    direction       = if_else(log2FC > 0, "Arterial > Venous", "Venous > Arterial")
  )
}) |>
  left_join(lx_ecs_map, by = "luminex_name") |>
  left_join(ecs_gene_table, by = "gene_symbol") |>
  arrange(p_value)

cat("\n--- Luminex: ECS-relevant cytokines (A-V paired Wilcoxon) ---\n")
lx_ecs_results |>
  select(luminex_name, gene_symbol, category, log2FC, mean_arterial, mean_venous,
         p_value, direction) |>
  mutate(across(c(log2FC, mean_arterial, mean_venous, p_value), ~ signif(., 3))) |>
  print(n = Inf)

write_csv(lx_ecs_results, file.path(out_tables, "08_ecs_luminex_av.csv"))

# ---------------------------------------------------------------------------
# 6. Plots
# ---------------------------------------------------------------------------

col_av    <- c("Arterial > Venous" = "#c0392b", "Venous > Arterial" = "#2980b9")
col_cat   <- setNames(
  RColorBrewer::brewer.pal(min(length(unique(ecs_gene_table$category)), 8), "Set2"),
  unique(ecs_gene_table$category)[seq_len(min(length(unique(ecs_gene_table$category)), 8))]
)

# ── Plot A: SomaScan A-V forest plot for ECS genes ──────────────────────────
if (nrow(ecs_soma_av) > 0) {
  p_forest <- ecs_soma_av |>
    mutate(
      label     = if_else(is.na(Target) | Target == "", EntrezGeneSymbol, Target),
      label     = fct_reorder(label, logFC),
      sig       = P.Value < 0.05,
      point_col = if_else(logFC > 0, "Arterial > Venous", "Venous > Arterial")
    ) |>
    ggplot(aes(x = logFC, y = label, color = point_col, size = sig)) +
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = if_else(sig, sprintf("p=%.3f", P.Value), "")),
              hjust = -0.15, size = 2.8, show.legend = FALSE) +
    scale_color_manual(values = col_av) +
    scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4), guide = "none") +
    facet_grid(category ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    labs(
      title    = "ECS-related proteins: SomaScan A-V gradient (N = 12 paired)",
      subtitle = "Positive logFC = higher in arterial blood. Larger points = p < 0.05.",
      x        = "log2FC (Arterial / Venous, paired limma)",
      y        = NULL,
      color    = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      strip.text.y     = element_text(size = 7, angle = 0),
      strip.placement  = "outside",
      panel.spacing    = unit(0.2, "lines")
    )

  ggsave(file.path(out_plots, "08_ecs_soma_forest.png"),
         p_forest, width = 9, height = 0.45 * nrow(ecs_soma_av) + 3,
         limitsize = FALSE, dpi = 150)
  cat("\nSaved: 08_ecs_soma_forest.png\n")
}

# ── Plot B: Luminex ECS cytokines dot plot ──────────────────────────────────
if (nrow(lx_ecs_results) > 0 && length(ecs_lx_names) > 0) {

  ecs_lx_long <- lx_avg |>
    select(patient, draw, all_of(ecs_lx_names)) |>
    pivot_longer(all_of(ecs_lx_names), names_to = "luminex_name", values_to = "mfi") |>
    left_join(lx_ecs_results |> select(luminex_name, gene_symbol, p_value, direction),
              by = "luminex_name") |>
    mutate(
      draw        = factor(draw, levels = c("Venous", "Arterial")),
      panel_label = paste0(gene_symbol, "\n(p=", signif(p_value, 2), ")")
    )

  p_lx_dots <- ggplot(ecs_lx_long,
                       aes(x = draw, y = mfi, group = patient, color = direction)) +
    geom_line(alpha = 0.4, linewidth = 0.7) +
    geom_point(size = 2.2, alpha = 0.85) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = col_av) +
    labs(
      title    = "ECS-relevant Luminex cytokines: paired arterial vs. venous (N = 12)",
      subtitle = "Lines connect the same patient. p = paired Wilcoxon.",
      x        = NULL, y = "MFI", color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  ggsave(file.path(out_plots, "08_ecs_luminex_paired.png"),
         p_lx_dots, width = 12, height = ceiling(length(ecs_lx_names) / 4) * 3.2 + 2,
         limitsize = FALSE, dpi = 150)
  cat("Saved: 08_ecs_luminex_paired.png\n")
}

# ── Plot C: Summary table — ECS hits ranked by evidence ─────────────────────
# Combine Luminex and SomaScan nominal hits into one ranked summary

soma_hits <- ecs_soma_av |>
  filter(P.Value < 0.05) |>
  transmute(
    platform     = "SomaScan",
    gene_symbol  = EntrezGeneSymbol,
    category,
    ecs_role,
    logFC,
    p_value      = P.Value,
    direction
  )

lx_hits <- lx_ecs_results |>
  filter(p_value < 0.05) |>
  transmute(
    platform     = "Luminex",
    gene_symbol,
    category,
    ecs_role,
    logFC        = log2FC,
    p_value,
    direction
  )

summary_hits <- bind_rows(soma_hits, lx_hits) |>
  arrange(p_value)

cat("\n========================================\n")
cat("SUMMARY: Nominal ECS hits (p < 0.05)\n")
cat("========================================\n")
if (nrow(summary_hits) > 0) {
  summary_hits |>
    mutate(across(c(logFC, p_value), ~ signif(., 3))) |>
    print(n = Inf)
} else {
  cat("No ECS gene reaches nominal significance on either platform.\n")
  cat("Full ECS gene list with p-values:\n")
  bind_rows(
    ecs_soma_av |> transmute(platform = "SomaScan", gene_symbol = EntrezGeneSymbol,
                              category, logFC, p_value = P.Value, direction),
    lx_ecs_results |> transmute(platform = "Luminex", gene_symbol,
                                 category, logFC = log2FC, p_value, direction)
  ) |>
    arrange(p_value) |>
    mutate(across(c(logFC, p_value), ~ signif(., 3))) |>
    print(n = 30)
}

write_csv(summary_hits, file.path(out_tables, "08_ecs_summary_hits.csv"))

# ── Plot D: Dot plot of all ECS genes in SomaScan ranked by –log10(p) ───────
if (nrow(ecs_soma_av) > 0) {
  p_pval <- ecs_soma_av |>
    mutate(
      label       = if_else(is.na(Target) | Target == "", EntrezGeneSymbol, Target),
      neg_log10_p = -log10(P.Value),
      sig         = P.Value < 0.05
    ) |>
    arrange(neg_log10_p) |>
    mutate(label = fct_inorder(label)) |>
    ggplot(aes(x = neg_log10_p, y = label, fill = category, alpha = sig)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1), guide = "none") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title    = "ECS genes in SomaScan A-V analysis: evidence ranking",
      subtitle = "Dashed = p < 0.05. Color = ECS functional category.",
      x        = "-log10(p-value, paired limma)",
      y        = NULL,
      fill     = "ECS category"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")

  ggsave(file.path(out_plots, "08_ecs_pval_rank.png"),
         p_pval, width = 10, height = 0.35 * nrow(ecs_soma_av) + 3,
         limitsize = FALSE, dpi = 150)
  cat("Saved: 08_ecs_pval_rank.png\n")
}

# ---------------------------------------------------------------------------
# 7. Endogenous opioid system (EOS): endorphins, enkephalins, dynorphins
# ---------------------------------------------------------------------------
# Note: β-endorphin (POMC-derived), met/leu-enkephalin (PENK-derived), and
# dynorphins (PDYN-derived) are peptides processed post-translationally.
# SomaScan v4.1 aptamers target the precursor proteins and soluble forms;
# the mature peptides themselves are typically below plasma detection limits.
# The degrading enzyme MME/neprilysin has a well-characterized circulating form.

eos_gene_table <- tribble(
  ~gene_symbol, ~eos_category,       ~eos_role,
  # Precursor proteins → mature peptides
  "POMC",       "Precursor",         "Pro-opiomelanocortin → β-endorphin, ACTH, α-MSH",
  "PENK",       "Precursor",         "Proenkephalin → met-enkephalin, leu-enkephalin",
  "PDYN",       "Precursor",         "Prodynorphin → dynorphin A/B, α-neoendorphin",
  "PNOC",       "Precursor",         "Pronociceptin → nociceptin/orphanin FQ (OFQ)",
  # Opioid receptors (mostly membrane-bound; circulating ectodomains rare but probed)
  "OPRM1",      "Receptor",          "μ-opioid receptor (MOR): β-endorphin, enkephalin target",
  "OPRD1",      "Receptor",          "δ-opioid receptor (DOR): enkephalin-preferring",
  "OPRK1",      "Receptor",          "κ-opioid receptor (KOR): dynorphin target",
  "OPRL1",      "Receptor",          "Nociceptin/OFQ receptor (NOP/ORL1)",
  # Enkephalin-degrading enzymes (soluble forms circulate in plasma)
  "MME",        "Degradation",       "Neprilysin (NEP): primary enkephalin-degrading endopeptidase",
  "ANPEP",      "Degradation",       "Aminopeptidase N (APN/CD13): N-terminal enkephalin cleavage",
  "ACE",        "Degradation",       "Angiotensin-converting enzyme: cleaves enkephalin C-terminus",
  "PREP",       "Degradation",       "Prolyl endopeptidase: neuropeptide degradation",
  # POMC-processing enzymes
  "PCSK1",      "Processing",        "PC1/3: POMC → β-endorphin / ACTH cleavage",
  "PCSK2",      "Processing",        "PC2: POMC → α-MSH / β-endorphin processing",
  "CPE",        "Processing",        "Carboxypeptidase E: opioid peptide maturation; circulates",
  # ECS–EOS crosstalk nodes
  "FAAH",       "ECS-EOS crosstalk", "FAAH degrades both AEA and some opioid-related lipids",
  "CNR1",       "ECS-EOS crosstalk", "CB1 and MOR co-expressed; heterodimer modulates opioid signaling",
  "BDNF",       "ECS-EOS crosstalk", "BDNF upregulated by opioids and ECS; pain sensitization",
  "SIGMAR1",    "ECS-EOS crosstalk", "Sigma-1 receptor: modulates both opioid and cannabinoid signaling",
  # Stress-axis context (perioperative opioid release is HPA-mediated)
  "CRH",        "HPA axis",          "Corticotropin-releasing hormone: triggers ACTH/β-endorphin from pituitary",
  "CGA",        "HPA axis",          "Glycoprotein hormone alpha chain: pituitary secretory marker"
)

cat("\n=== ENDOGENOUS OPIOID SYSTEM ANALYSIS ===\n")
cat("EOS gene universe:", nrow(eos_gene_table), "entries\n")

# ── 7a. SomaScan A-V ──────────────────────────────────────────────────────────

eos_soma_av <- av_soma |>
  inner_join(eos_gene_table, by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  arrange(P.Value) |>
  select(Target, EntrezGeneSymbol, eos_category, eos_role,
         logFC, AveExpr, t, P.Value, adj.P.Val, direction)

cat("\n--- SomaScan A-V: EOS genes found ---\n")
cat("Total EOS genes detected:", nrow(eos_soma_av), "\n")
cat("Nominal p < 0.05:", sum(eos_soma_av$P.Value < 0.05, na.rm = TRUE), "\n\n")
print(eos_soma_av, n = Inf)

write_csv(eos_soma_av, file.path(out_tables, "08_eos_somascan_av.csv"))

# ── 7b. SomaScan covariate ────────────────────────────────────────────────────

eos_soma_cov <- cov_soma |>
  inner_join(eos_gene_table, by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  arrange(P.Value) |>
  select(covariate, Target, EntrezGeneSymbol, eos_category, eos_role,
         logFC, t, P.Value, adj.P.Val)

cat("\n--- SomaScan covariate: EOS genes ---\n")
cat("Nominal p < 0.05:", sum(eos_soma_cov$P.Value < 0.05, na.rm = TRUE), "\n\n")
eos_soma_cov |> filter(P.Value < 0.05) |> print(n = Inf)

write_csv(eos_soma_cov, file.path(out_tables, "08_eos_somascan_covariate.csv"))

# ── 7c. SomaScan A-V delta × covariate ───────────────────────────────────────

eos_soma_delta <- delta_soma |>
  inner_join(eos_gene_table, by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  arrange(P.Value) |>
  select(covariate, Target, EntrezGeneSymbol, eos_category, eos_role,
         logFC, t, P.Value, adj.P.Val)

cat("\n--- SomaScan A-V delta × covariate: EOS genes ---\n")
cat("Nominal p < 0.05:", sum(eos_soma_delta$P.Value < 0.05, na.rm = TRUE), "\n\n")
eos_soma_delta |> filter(P.Value < 0.05) |> print(n = Inf)

write_csv(eos_soma_delta, file.path(out_tables, "08_eos_somascan_delta.csv"))

# ── 7d. Combined EOS summary ──────────────────────────────────────────────────

eos_summary <- bind_rows(
  eos_soma_av |> filter(P.Value < 0.05) |>
    transmute(analysis = "A-V gradient", platform = "SomaScan",
              gene_symbol = EntrezGeneSymbol, eos_category, eos_role,
              logFC, p_value = P.Value, direction),
  eos_soma_cov |> filter(P.Value < 0.05) |>
    transmute(analysis = paste0("Covariate: ", covariate), platform = "SomaScan",
              gene_symbol = EntrezGeneSymbol, eos_category, eos_role,
              logFC, p_value = P.Value, direction = if_else(logFC > 0, "Up in Yes", "Up in No")),
  eos_soma_delta |> filter(P.Value < 0.05) |>
    transmute(analysis = paste0("Delta: ", covariate), platform = "SomaScan",
              gene_symbol = EntrezGeneSymbol, eos_category, eos_role,
              logFC, p_value = P.Value, direction = if_else(logFC > 0, "Larger in Yes", "Smaller in Yes"))
) |>
  arrange(p_value)

cat("\n========================================\n")
cat("SUMMARY: EOS nominal hits (p < 0.05)\n")
cat("========================================\n")
if (nrow(eos_summary) > 0) {
  eos_summary |>
    mutate(across(c(logFC, p_value), ~ signif(., 3))) |>
    print(n = Inf)
} else {
  cat("No EOS gene reaches nominal significance in any analysis.\n")
  cat("\nFull EOS A-V results (all genes detected, sorted by p):\n")
  eos_soma_av |>
    select(Target, EntrezGeneSymbol, eos_category, logFC, P.Value, direction) |>
    mutate(across(c(logFC, P.Value), ~ signif(., 3))) |>
    print(n = Inf)
}

write_csv(eos_summary, file.path(out_tables, "08_eos_summary_hits.csv"))

# ── 7e. Forest plot: EOS genes in SomaScan A-V ───────────────────────────────

if (nrow(eos_soma_av) > 0) {
  p_eos_forest <- eos_soma_av |>
    mutate(
      label     = if_else(is.na(Target) | Target == "", EntrezGeneSymbol, Target),
      label     = fct_reorder(label, logFC),
      sig       = P.Value < 0.05,
      point_col = if_else(logFC > 0, "Arterial > Venous", "Venous > Arterial")
    ) |>
    ggplot(aes(x = logFC, y = label, color = point_col, size = sig)) +
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = if_else(sig, sprintf("p=%.3f", P.Value), "")),
              hjust = -0.15, size = 2.8, show.legend = FALSE) +
    scale_color_manual(values = col_av) +
    scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4), guide = "none") +
    facet_grid(eos_category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(
      title    = "Endogenous opioid system proteins: SomaScan A-V gradient (N = 12 paired)",
      subtitle = "Positive logFC = higher in arterial. Larger points = p < 0.05.",
      x        = "log2FC (Arterial / Venous, paired limma)",
      y        = NULL, color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text.y    = element_text(size = 7, angle = 0),
      strip.placement = "outside",
      panel.spacing   = unit(0.2, "lines")
    )

  ggsave(file.path(out_plots, "08_eos_soma_forest.png"),
         p_eos_forest,
         width = 9, height = 0.45 * nrow(eos_soma_av) + 3,
         limitsize = FALSE, dpi = 150)
  cat("\nSaved: 08_eos_soma_forest.png\n")
}

cat("\nDone. Outputs in outputs/tables/08_*.csv and outputs/plots/08_*.png\n")
