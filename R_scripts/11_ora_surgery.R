# 11_ora_surgery.R
# Overrepresentation analysis (Fisher's exact) for:
#   (1) Surgery-type pairwise contrasts (Arterial / Venous / A-V Delta)
#   (2) Binary covariates (obese, hypothermia, hyperglycemia, sex) × draw
# Runs ORA separately for up- and down-regulated hits in each comparison,
# using Hallmark + KEGG gene sets from msigdbr.
# Outputs:
#   tables/11_ora_surgery_results.csv   — surgery ORA, all results
#   tables/11_ora_surgery_top.csv       — surgery ORA, p < 0.05
#   tables/11_ora_covariate_results.csv — covariate ORA, all results
#   tables/11_ora_covariate_top.csv     — covariate ORA, p < 0.05
#   plots/11_ora_surgery_heatmap.pdf
#   plots/11_ora_covariate_heatmap.pdf

data_dir <- tryCatch(
  normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)
source(file.path(data_dir, "R_scripts", "00_setup.R"))
library(msigdbr)

# ── Gene sets ──────────────────────────────────────────────────────────────────
h_sets <- msigdbr(species = "Homo sapiens", category = "H") |>
  select(gs_name, gene_symbol)
k_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY") |>
  select(gs_name, gene_symbol)
gene_sets <- bind_rows(
  h_sets |> mutate(collection = "Hallmark"),
  k_sets |> mutate(collection = "KEGG")
)

# ── ORA function ───────────────────────────────────────────────────────────────
run_ora <- function(hit_syms, background_syms, gene_sets_df, min_overlap = 3) {
  if (length(hit_syms) == 0) return(tibble())
  sets <- split(gene_sets_df$gene_symbol, gene_sets_df$gs_name)
  coll <- gene_sets_df |> distinct(gs_name, collection) |> deframe()

  map(names(sets), function(gs) {
    in_set  <- background_syms %in% sets[[gs]]
    is_hit  <- background_syms %in% hit_syms
    a <- sum( is_hit &  in_set)
    b <- sum( is_hit & !in_set)
    c <- sum(!is_hit &  in_set)
    d <- sum(!is_hit & !in_set)
    if (a < min_overlap) return(NULL)
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    tibble(
      pathway    = gs,
      collection = coll[[gs]],
      n_hit_in_set  = a,
      n_set      = a + c,
      n_hit      = a + b,
      n_background  = length(background_syms),
      odds_ratio = as.numeric(ft$estimate),
      p_value    = ft$p.value
    )
  }) |>
    bind_rows() -> out
  if (nrow(out) == 0) return(tibble())
  out |>
    mutate(padj = p.adjust(p_value, method = "BH")) |>
    arrange(p_value)
}

# ── Load pairwise results ──────────────────────────────────────────────────────
pw_art   <- read_csv(file.path(table_dir, "05_arterial_surgery_pairwise.csv"),
                     show_col_types = FALSE)
pw_ven   <- read_csv(file.path(table_dir, "05_venous_surgery_pairwise.csv"),
                     show_col_types = FALSE)
pw_delta <- read_csv(file.path(table_dir, "06_delta_surgery_pairwise.csv"),
                     show_col_types = FALSE)

all_pw <- bind_rows(
  pw_art   |> mutate(draw = "Arterial"),
  pw_ven   |> mutate(draw = "Venous"),
  pw_delta |> mutate(draw = "AV_Delta")
)

contrasts <- unique(all_pw$contrast)
draws     <- unique(all_pw$draw)

# ── Run ORA for each contrast × draw × direction ───────────────────────────────
contrast_labels <- c(
  Spine_vs_Laparoscopic    = "Spine / Laparoscopic",
  Spine_vs_HeadNeck        = "Spine / Head-Neck",
  Laparoscopic_vs_HeadNeck = "Laparoscopic / Head-Neck"
)

ora_all <- map(draws, function(dr) {
  map(contrasts, function(ct) {
    df <- all_pw |> dplyr::filter(draw == dr, contrast == ct)
    bg <- df$EntrezGeneSymbol

    map(c(up = "up", down = "down"), function(dir) {
      hits <- if (dir == "up") {
        df |> dplyr::filter(P.Value < 0.05, logFC > 0) |> pull(EntrezGeneSymbol)
      } else {
        df |> dplyr::filter(P.Value < 0.05, logFC < 0) |> pull(EntrezGeneSymbol)
      }
      if (length(hits) < 3) return(tibble())

      run_ora(hits, bg, gene_sets) |>
        mutate(
          draw      = dr,
          contrast  = ct,
          direction = dir,
          n_sig_hits = length(hits),
          label     = paste0(dr, " · ", contrast_labels[ct],
                             " · ", if (dir == "up") "higher in 1st" else "higher in 2nd")
        )
    }) |> bind_rows()
  }) |> bind_rows()
}) |> bind_rows()

cat("\nORA hits (p < 0.05) per draw × contrast × direction:\n")
print(ora_all |> dplyr::filter(p_value < 0.05) |>
      count(draw, contrast, direction) |> arrange(draw, contrast, direction))

write_csv(ora_all, file.path(table_dir, "11_ora_surgery_results.csv"))
write_csv(
  ora_all |> dplyr::filter(p_value < 0.05) |> arrange(draw, contrast, direction, p_value),
  file.path(table_dir, "11_ora_surgery_top.csv")
)

# ── Heatmap: –log10(p) across contrast × direction columns ────────────────────
ora_sig <- ora_all |>
  dplyr::filter(p_value < 0.05) |>
  mutate(pathway_clean = pathway |>
           str_remove("^HALLMARK_") |>
           str_remove("^KEGG_") |>
           str_replace_all("_", " ") |>
           str_to_title())

save_ora_heatmap <- function(ora_df, out_path, title) {
  sig <- ora_df |>
    dplyr::filter(p_value < 0.05) |>
    mutate(pathway_clean = pathway |>
             str_remove("^HALLMARK_") |>
             str_remove("^KEGG_") |>
             str_replace_all("_", " ") |>
             str_to_title())

  if (nrow(sig) == 0) {
    cat("No ORA hits at p < 0.05 for", title, "— skipping heatmap.\n")
    return(invisible(NULL))
  }

  hm_df <- sig |>
    mutate(neg_log10p = pmin(-log10(p_value), 4)) |>
    group_by(pathway_clean, label) |>
    summarise(neg_log10p = max(neg_log10p), .groups = "drop") |>
    pivot_wider(names_from = label, values_from = neg_log10p, values_fill = 0) |>
    column_to_rownames("pathway_clean") |>
    as.matrix()

  pheatmap(
    hm_df,
    color        = colorRampPalette(c("white", "#f7dc6f", "#c0392b"))(50),
    breaks       = seq(0, 4, length.out = 51),
    cluster_rows = TRUE, cluster_cols = FALSE,
    fontsize_row = 7, fontsize_col = 7,
    main         = paste0("ORA: ", title, " | -log10(p) (capped at 4)"),
    border_color = NA, angle_col = 45,
    filename     = out_path,
    width = 16, height = max(6, nrow(hm_df) * 0.2 + 3)
  )
  cat("Heatmap saved:", nrow(hm_df), "pathways ×", ncol(hm_df), "columns.\n")
}

save_ora_heatmap(ora_all, file.path(plot_dir, "11_ora_surgery_heatmap.pdf"),
                 "Surgery-type contrasts")

# ══════════════════════════════════════════════════════════════════════════════
# Binary covariate ORA
# ══════════════════════════════════════════════════════════════════════════════
cov_art_res   <- read_csv(file.path(table_dir, "05_arterial_results.csv"),
                           show_col_types = FALSE) |> mutate(draw = "Arterial")
cov_ven_res   <- read_csv(file.path(table_dir, "05_venous_results.csv"),
                           show_col_types = FALSE) |> mutate(draw = "Venous")
cov_delta_res <- read_csv(file.path(table_dir, "06_delta_results.csv"),
                           show_col_types = FALSE) |> mutate(draw = "AV_Delta")

all_cov <- bind_rows(cov_art_res, cov_ven_res, cov_delta_res)

cov_labels <- c(
  obese         = "Obese / Non-obese",
  hypothermia   = "Hypothermia / Normothermia",
  hyperglycemia = "Hyperglycemia / Normoglycemia",
  male          = "Male / Female"
)

cat("\nRunning binary covariate ORA...\n")

ora_cov <- map(unique(all_cov$draw), function(dr) {
  map(names(cov_labels), function(cov) {
    df <- all_cov |> dplyr::filter(draw == dr, covariate == cov)
    if (nrow(df) == 0) return(tibble())
    bg <- df$EntrezGeneSymbol

    map(c(up = "up", down = "down"), function(dir) {
      hits <- if (dir == "up") {
        df |> dplyr::filter(P.Value < 0.05, logFC > 0) |> pull(EntrezGeneSymbol)
      } else {
        df |> dplyr::filter(P.Value < 0.05, logFC < 0) |> pull(EntrezGeneSymbol)
      }
      if (length(hits) < 3) return(tibble())

      run_ora(hits, bg, gene_sets) |>
        mutate(
          draw      = dr,
          covariate = cov,
          direction = dir,
          n_sig_hits = length(hits),
          label     = paste0(dr, " · ", cov_labels[cov],
                             " · ", if (dir == "up") "higher in pos" else "higher in neg")
        )
    }) |> bind_rows()
  }) |> bind_rows()
}) |> bind_rows()

cat("\nORA hits (p < 0.05) per draw × covariate × direction:\n")
print(ora_cov |> dplyr::filter(p_value < 0.05) |>
      count(draw, covariate, direction) |> arrange(draw, covariate, direction))

write_csv(ora_cov, file.path(table_dir, "11_ora_covariate_results.csv"))
write_csv(
  ora_cov |> dplyr::filter(p_value < 0.05) |> arrange(draw, covariate, direction, p_value),
  file.path(table_dir, "11_ora_covariate_top.csv")
)

save_ora_heatmap(ora_cov, file.path(plot_dir, "11_ora_covariate_heatmap.pdf"),
                 "Binary covariates")

cat("\nAll outputs saved to outputs/plots/ and outputs/tables/\n")
