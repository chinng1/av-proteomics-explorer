# 04_gsea.R
# Gene-set enrichment analysis (GSEA) on the SomaScan A-V t-statistic ranking.
# Tests Hallmark and KEGG gene sets, plus a custom 65-cytokine panel gene set.
# Outputs: plots/04_gsea_barplot.png, 04_panel_enrichment.png, 04_panel_rank.png
#          tables/04_gsea_hallmark.csv, 04_gsea_kegg.csv, 04_gsea_panel.csv

data_dir <- tryCatch(
  normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)
source(file.path(data_dir, "R_scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(fgsea); library(msigdbr) })

# ── Re-run limma to get av_results (needed for ranking) ───────────────────────
source(file.path(data_dir, "R_scripts", "03_av_limma.R"))

# ── Build ranked vector (Entrez gene ID → t-statistic) ────────────────────────
rank_df <- av_results |>
  dplyr::filter(!is.na(EntrezGeneID), EntrezGeneID != "") |>
  mutate(EntrezGeneID = as.character(EntrezGeneID)) |>
  group_by(EntrezGeneID) |>
  slice_max(abs(t), n = 1, with_ties = FALSE) |>
  ungroup()

ranked_vec <- setNames(rank_df$t, rank_df$EntrezGeneID) |> sort(decreasing = TRUE)
cat("Unique Entrez IDs in ranking:", length(ranked_vec), "\n")

set.seed(42)

# ── Hallmark ───────────────────────────────────────────────────────────────────
h_sets <- msigdbr(species = "Homo sapiens", category = "H") |>
  select(gs_name, entrez_gene) |>
  mutate(entrez_gene = as.character(entrez_gene)) |>
  group_by(gs_name) |> summarise(genes = list(entrez_gene), .groups = "drop")
hallmark_list <- setNames(h_sets$genes, h_sets$gs_name)

gsea_hallmark <- fgsea(pathways = hallmark_list, stats = ranked_vec,
                        minSize = 10, maxSize = 500, nPermSimple = 10000) |>
  as_tibble() |>
  mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";")) |>
  arrange(pval)

cat("Hallmark nominal (p<0.05):", sum(gsea_hallmark$pval < 0.05), "\n")
write_csv(gsea_hallmark, file.path(table_dir, "04_gsea_hallmark.csv"))

# ── KEGG ───────────────────────────────────────────────────────────────────────
kegg_sets <- msigdbr(species = "Homo sapiens", category = "C2",
                      subcategory = "CP:KEGG_LEGACY") |>
  select(gs_name, entrez_gene) |>
  mutate(entrez_gene = as.character(entrez_gene)) |>
  group_by(gs_name) |> summarise(genes = list(entrez_gene), .groups = "drop")
kegg_list <- setNames(kegg_sets$genes, kegg_sets$gs_name)

gsea_kegg <- fgsea(pathways = kegg_list, stats = ranked_vec,
                    minSize = 10, maxSize = 500, nPermSimple = 10000) |>
  as_tibble() |>
  mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";")) |>
  arrange(pval)

cat("KEGG nominal (p<0.05):", sum(gsea_kegg$pval < 0.05), "\n")
write_csv(gsea_kegg, file.path(table_dir, "04_gsea_kegg.csv"))

# ── Combined barplot ───────────────────────────────────────────────────────────
all_gsea <- bind_rows(
  gsea_hallmark |> mutate(collection = "Hallmark"),
  gsea_kegg     |> mutate(collection = "KEGG")
) |>
  dplyr::filter(pval < 0.05) |>
  mutate(
    short_name = str_replace_all(pathway, "HALLMARK_|KEGG_", "") |>
                  str_replace_all("_", " ") |> str_to_title(),
    direction  = if_else(NES > 0, "Enriched in Arterial", "Enriched in Venous")
  ) |>
  arrange(NES)

if (nrow(all_gsea) > 0) {
  p_gsea <- ggplot(all_gsea,
         aes(x = NES, y = reorder(short_name, NES),
             fill = direction, alpha = -log10(pval))) +
    geom_col() +
    geom_vline(xintercept = 0, color = "gray30") +
    scale_fill_manual(values = c("Enriched in Arterial" = "#c0392b",
                                  "Enriched in Venous"   = "#2980b9")) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    facet_wrap(~ collection, scales = "free_y") +
    labs(
      title = "GSEA: A-V ranking — nominal hits (p < 0.05)",
      x = "Normalized Enrichment Score (positive = enriched in arterial)",
      y = NULL, fill = NULL
    ) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  ggsave(file.path(plot_dir, "04_gsea_barplot.png"), p_gsea,
         width = 14, height = 9, dpi = 300)
}

# ── Custom 65-cytokine panel gene set ─────────────────────────────────────────
panel_entrez <- av_results |>
  dplyr::filter(EntrezGeneSymbol %in% luminex_gene_map$gene_symbol,
                !is.na(EntrezGeneID), EntrezGeneID != "") |>
  pull(EntrezGeneID) |> as.character() |> unique()

cat("Panel Entrez IDs recovered:", length(panel_entrez), "\n")

gsea_panel <- fgsea(
  pathways    = list(Luminex_65plex_Cytokine_Panel = panel_entrez),
  stats       = ranked_vec, minSize = 5, maxSize = 500, nPermSimple = 100000
) |>
  as_tibble() |>
  mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))

cat("Panel GSEA: NES =", round(gsea_panel$NES, 3),
    "p =", signif(gsea_panel$pval, 3), "\n")
write_csv(gsea_panel, file.path(table_dir, "04_gsea_panel.csv"))

# ── Panel enrichment plot ──────────────────────────────────────────────────────
p_enrich <- plotEnrichment(
  list(Luminex_65plex_Cytokine_Panel = panel_entrez)[["Luminex_65plex_Cytokine_Panel"]],
  ranked_vec
) +
  labs(
    title    = "GSEA enrichment: Luminex 65-plex cytokine panel",
    subtitle = sprintf("NES = %.2f, p = %.3f — positive NES = enriched in arterial",
                       gsea_panel$NES, gsea_panel$pval),
    x = "Rank in SomaScan A-V t-statistic (left = more arterial)",
    y = "Enrichment score"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(plot_dir, "04_panel_enrichment.png"), p_enrich,
       width = 9, height = 5, dpi = 300)

# ── Panel rank plot ────────────────────────────────────────────────────────────
lx_raw <- read_csv(file.path(data_dir, "AP_65plex_byJMR.csv"),
                    skip = 54, show_col_types = FALSE) |>
  dplyr::filter(str_starts(Sample, "Unknown"))
lx_cols <- lx_raw |> select(-Location, -Sample, -`Total Events`) |> names()

plate_map <- data.frame(
  unknown_id = paste0("Unknown", 1:24),
  patient    = c("0179","0193","0209","0215","0179","0193","0209","0215",
                 "0223","0245","0250","0254","0223","0245","0250","0254",
                 "0268","0275","0284","0291","0268","0275","0284","0291"),
  draw       = c(rep("Venous",4),rep("Arterial",4),rep("Venous",4),rep("Arterial",4),
                 rep("Venous",4),rep("Arterial",4))
)
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
         log2FC_lx   = log2((mean(a, na.rm=TRUE)+0.1)/(mean(v, na.rm=TRUE)+0.1)),
         p_value_lx  = test$p.value)
})

panel_ranks <- av_results |>
  mutate(rank = row_number()) |>
  dplyr::filter(EntrezGeneSymbol %in% luminex_gene_map$gene_symbol) |>
  left_join(luminex_gene_map, by = c("EntrezGeneSymbol" = "gene_symbol")) |>
  left_join(lx_results |> select(luminex_name, log2FC_lx, p_value_lx),
            by = "luminex_name") |>
  mutate(
    sig_lx   = p_value_lx < 0.05,
    sig_soma = P.Value    < 0.05,
    label    = if_else(sig_lx | sig_soma, EntrezGeneSymbol, NA_character_)
  )

p_rank <- ggplot(panel_ranks,
       aes(x = rank, y = logFC, color = sig_lx, shape = sig_soma, label = label)) +
  geom_hline(yintercept = 0, color = "gray60") +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(size = 2.8, max.overlaps = 20, box.padding = 0.3,
                  na.rm = TRUE, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "#e67e22", "FALSE" = "gray60"),
                     labels = c("TRUE" = "Luminex p<0.05", "FALSE" = "Luminex p≥0.05")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                     labels = c("TRUE" = "SomaScan p<0.05", "FALSE" = "SomaScan p≥0.05")) +
  scale_x_continuous(sec.axis = sec_axis(~./nrow(av_results)*100, name = "Percentile")) +
  labs(
    title    = "Position of 65 Luminex cytokines in the SomaScan A-V ranking",
    subtitle = "Positive logFC = higher in arterial. Rank 1 = most arterial.",
    x = "Rank in SomaScan proteome (1 = most arterial)",
    y = "SomaScan log2FC", color = NULL, shape = NULL
  ) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "04_panel_rank.png"), p_rank, width = 11, height = 6, dpi = 300)

cat("\nOutputs saved to outputs/plots/ and outputs/tables/\n")
