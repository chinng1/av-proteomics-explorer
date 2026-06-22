# Perioperative A-V Proteomics Explorer
# Run with: shiny::runApp("shiny_explorer")  from project root
# Or open this file in RStudio and click "Run App"

pkgs <- c("shiny", "bslib", "tidyverse", "plotly", "reactable")
missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install missing packages first:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(tidyverse)
  library(plotly)
  library(reactable)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
# Shiny sets the working directory to the app folder (shiny_explorer/) at launch
proj_dir <- normalizePath("..")
tbl_dir  <- file.path(proj_dir, "outputs", "tables")

# ── Data loading ──────────────────────────────────────────────────────────────
av        <- read_csv(file.path(tbl_dir, "03_av_results_full.csv"),    show_col_types = FALSE)
loadings  <- read_csv(file.path(tbl_dir, "02_pca_loadings.csv"),        show_col_types = FALSE)
pca_var   <- read_csv(file.path(tbl_dir, "02_pca_variance.csv"),        show_col_types = FALSE)
gsea_hall <- read_csv(file.path(tbl_dir, "04_gsea_hallmark.csv"),       show_col_types = FALSE)
gsea_kegg <- read_csv(file.path(tbl_dir, "04_gsea_kegg.csv"),           show_col_types = FALSE)
ecs_hits  <- read_csv(file.path(tbl_dir, "08_ecs_summary_hits.csv"),    show_col_types = FALSE)
eos8_hits <- read_csv(file.path(tbl_dir, "08_eos_summary_hits.csv"),    show_col_types = FALSE)
eos9_hits <- read_csv(file.path(tbl_dir, "09_eos_summary_hits.csv"),    show_col_types = FALSE)
cov_res   <- read_csv(file.path(tbl_dir, "05_covariate_results.csv"),   show_col_types = FALSE)

scores_path <- file.path(tbl_dir, "02_pca_scores.csv")
pca_scores  <- if (file.exists(scores_path)) read_csv(scores_path, show_col_types = FALSE) else NULL

skin_av  <- if (file.exists(file.path(tbl_dir, "10_skin_av.csv")))
              read_csv(file.path(tbl_dir, "10_skin_av.csv"),         show_col_types = FALSE) else NULL
skin_kw  <- if (file.exists(file.path(tbl_dir, "10_skin_surgery_kw.csv")))
              read_csv(file.path(tbl_dir, "10_skin_surgery_kw.csv"), show_col_types = FALSE) else NULL

# ── Normalize A-V results ──────────────────────────────────────────────────────
av <- av |>
  mutate(
    neg_log10p   = -log10(P.Value),
    EntrezGeneID = as.character(EntrezGeneID),
    sig_label    = case_when(
      sig_fdr ~ "FDR < 0.05",
      sig_nom ~ "p < 0.05",
      TRUE    ~ "NS"
    )
  )

# ── Normalize ECS / EOS into a unified table ───────────────────────────────────
ecs_norm <- ecs_hits |>
  mutate(system = "ECS", analysis = "A-V") |>
  select(system, analysis, platform, gene_symbol, category,
         role = ecs_role, logFC, p_value, direction)

eos8_norm <- eos8_hits |>
  mutate(system = "EOS") |>
  select(system, analysis, platform, gene_symbol,
         category = eos_category, role = eos_role, logFC, p_value, direction)

eos9_norm <- eos9_hits |>
  mutate(system = "EOS", platform = "SomaScan",
         analysis = paste0(analysis, ": ", covariate)) |>
  select(system, analysis, platform,
         gene_symbol = EntrezGeneSymbol,
         category = eos_category, role = eos_role, logFC, p_value, direction)

ecs_eos <- bind_rows(ecs_norm, eos8_norm, eos9_norm) |>
  distinct(system, analysis, platform, gene_symbol, .keep_all = TRUE)

# ── GSEA combined ─────────────────────────────────────────────────────────────
gsea_all <- bind_rows(
  mutate(gsea_hall, collection = "Hallmark"),
  mutate(gsea_kegg, collection = "KEGG")
) |>
  arrange(padj)

parse_leading_edge <- function(le_str) {
  if (is.na(le_str) || le_str == "") return(character(0))
  strsplit(le_str, ";")[[1]]
}

# ── Shared reactable style ─────────────────────────────────────────────────────
rtbl <- function(df, ...) {
  reactable(df,
    highlight = TRUE, striped = TRUE, compact = TRUE,
    defaultPageSize = 12, searchable = TRUE,
    theme = reactableTheme(
      stripedColor = "#f7f9fc",
      highlightColor = "#e8f0fe"
    ),
    ...
  )
}

# ── UI ────────────────────────────────────────────────────────────────────────
gene_search_ui <- div(
  style = "display:flex; align-items:center; gap:6px; padding-right:8px;",
  tags$span(style = "color:#aaa; font-size:12px; white-space:nowrap;", "Gene search:"),
  textInput("gene_search", NULL,
            placeholder = "e.g. IL6  TNF  CNR1",
            width = "290px"),
  actionButton("clear_search", "×",
               class = "btn-sm btn-outline-secondary",
               style = "padding:2px 8px; line-height:1.4;")
)

ui <- page_navbar(
  title = tags$span(
    tags$b("A-V Proteomics"),
    tags$span(style = "font-size:13px; color:#888;", " · 7,481 proteins · 12 patients")
  ),
  theme = bs_theme(bootswatch = "flatly", primary = "#2c7bb6"),
  nav_spacer(),
  nav_item(gene_search_ui),

  # ── Tab 1: A-V Volcano ─────────────────────────────────────────────────────
  nav_panel(
    "A-V Volcano",
    layout_sidebar(
      sidebar = sidebar(
        width = 265, open = TRUE,
        h6("Display options"),
        selectInput("sig_filter", "Show proteins",
                    choices = c("All (7,481)" = "all",
                                "Nominal p < 0.05" = "nom",
                                "FDR < 0.05" = "fdr")),
        selectInput("dir_filter", "Direction",
                    choices = c("Both", "Arterial > Venous", "Venous > Arterial")),
        hr(),
        h6("Searched gene(s) — A-V results"),
        verbatimTextOutput("gene_av_summary"),
        hr(),
        h6("Also appears in:"),
        verbatimTextOutput("gene_cross_summary")
      ),
      card(plotlyOutput("volcano_plot", height = "520px")),
      card(
        card_header("Nominally significant proteins (p < 0.05) — or gene search results"),
        reactableOutput("av_table")
      )
    )
  ),

  # ── Tab 2: PCA ─────────────────────────────────────────────────────────────
  nav_panel(
    "PCA",
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Protein loadings — PC1 vs PC2"),
        plotlyOutput("pca_loadings_plot", height = "460px"),
        p(class = "text-muted small mt-1",
          "Each point = one aptamer. Searched genes highlighted in red.")
      ),
      card(
        card_header("Top PC1 contributors"),
        reactableOutput("pca_top_table")
      )
    ),
    if (!is.null(pca_scores)) {
      card(
        card_header("Sample scores (PC1 vs PC2)"),
        plotlyOutput("pca_scores_plot", height = "380px"),
        p(class = "text-muted small", "Lines connect A-V pairs per patient.")
      )
    } else {
      card(
        card_header("Sample scores"),
        p(class = "text-warning",
          "Re-run R_scripts/02_pca.R to enable the interactive sample PCA.",
          " Static image: outputs/plots/02_pca_pc1_pc2.png")
      )
    }
  ),

  # ── Tab 3: GSEA ────────────────────────────────────────────────────────────
  nav_panel(
    "GSEA",
    layout_sidebar(
      sidebar = sidebar(
        width = 265,
        h6("Filter pathways"),
        selectInput("gsea_collection", "Collection",
                    choices = c("All", "Hallmark", "KEGG")),
        selectInput("gsea_dir", "Direction",
                    choices = c("Both",
                                "Suppressed in arterial (NES < 0)" = "neg",
                                "Enriched in arterial (NES > 0)" = "pos")),
        sliderInput("gsea_padj", "Max adj.p", min = 0.001, max = 1,
                    value = 0.25, step = 0.005),
        hr(),
        h6("Leading-edge genes"),
        verbatimTextOutput("leading_edge_genes")
      ),
      card(
        card_header("Pathways — click a row to highlight its leading-edge genes on the volcano below"),
        reactableOutput("gsea_table")
      ),
      card(
        card_header("Leading-edge genes on A-V volcano"),
        plotlyOutput("gsea_volcano", height = "420px")
      )
    )
  ),

  # ── Tab 4: Skin Proteases ──────────────────────────────────────────────────
  nav_panel(
    "Skin Proteases",
    layout_sidebar(
      sidebar = sidebar(
        width = 265,
        h6("Filter"),
        selectInput("skin_family", "Gene family",
                    choices = c("All", "SPINK", "KLK", "SERPIN")),
        hr(),
        p(class = "text-muted small",
          tags$b("SPINK9"), " — top A-V hit (rank 1, p = 1.7×10⁻⁴)", br(),
          tags$b("KLK7"),   " — rank 3 (p = 1.9×10⁻³)", br(), br(),
          "All nominally significant members are venous > arterial, consistent with
          wound-site release of epidermal proteins into venous return.")
      ),
      card(
        card_header("Skin serine protease family: A-V effect sizes (larger dot = p < 0.05)"),
        plotlyOutput("skin_forest", height = "520px")
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(
          card_header("All family members detected in SomaScan"),
          reactableOutput("skin_table")
        ),
        card(
          card_header("Surgery-type gradient (Kruskal-Wallis, exploratory)"),
          if (!is.null(skin_kw))
            reactableOutput("skin_kw_table")
          else
            p(class = "text-muted",
              "Run R_scripts/10_skin_protease_analysis.R to generate this table.")
        )
      )
    )
  ),

  # ── Tab 5: ECS / EOS ───────────────────────────────────────────────────────
  nav_panel(
    "ECS / EOS",
    layout_sidebar(
      sidebar = sidebar(
        width = 265,
        h6("Filter"),
        selectInput("ecs_system", "System",
                    choices = c("Both" = "both", "ECS", "EOS")),
        selectInput("ecs_platform", "Platform",
                    choices = c("All" = "all", "Luminex", "SomaScan")),
        hr(),
        p(class = "text-muted small",
          "Click a row to locate that gene on the volcano.",
          br(),
          tags$b("Orange:"), " all ECS/EOS genes", br(),
          tags$b("Red:"), " selected row", br(),
          tags$b("Purple:"), " gene search")
      ),
      card(
        card_header("ECS / EOS hits across analyses — click a row to locate on volcano"),
        reactableOutput("ecs_table")
      ),
      card(
        card_header("ECS / EOS gene positions on A-V volcano"),
        plotlyOutput("ecs_volcano", height = "440px")
      )
    )
  )
)

# ── Server ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Gene search ─────────────────────────────────────────────────────────────
  observeEvent(input$clear_search, {
    updateTextInput(session, "gene_search", value = "")
  })

  searched_genes <- reactive({
    raw <- trimws(input$gene_search)
    if (raw == "") return(character(0))
    toupper(strsplit(raw, "[,;\\s]+", perl = TRUE)[[1]])
  })

  av_marked <- reactive({
    genes <- searched_genes()
    av |> mutate(is_hit = toupper(EntrezGeneSymbol) %in% genes)
  })

  # ── Sidebar cross-analysis summaries ──────────────────────────────────────
  output$gene_av_summary <- renderText({
    genes <- searched_genes()
    if (length(genes) == 0) return("(type gene symbol above)")
    hits <- av |> filter(toupper(EntrezGeneSymbol) %in% genes) |>
      arrange(P.Value)
    if (nrow(hits) == 0) return(paste("Not found:", paste(genes, collapse = ", ")))
    paste(apply(hits, 1, function(r) {
      paste0(r["EntrezGeneSymbol"],
             "\n  logFC = ", round(as.numeric(r["logFC"]), 3),
             "\n  p = ", signif(as.numeric(r["P.Value"]), 3),
             "\n  ", r["direction"])
    }), collapse = "\n\n")
  })

  output$gene_cross_summary <- renderText({
    genes <- searched_genes()
    if (length(genes) == 0) return("")
    lines <- character(0)
    cov_hit <- cov_res |> filter(toupper(EntrezGeneSymbol) %in% genes)
    if (nrow(cov_hit) > 0)
      lines <- c(lines, paste0("Covariate: ", paste(unique(cov_hit$covariate), collapse = ", ")))
    ecs_hit <- ecs_eos |> filter(toupper(gene_symbol) %in% genes)
    if (nrow(ecs_hit) > 0)
      lines <- c(lines, paste0(paste(unique(ecs_hit$system), collapse = "/"), " panel"))
    gsea_le <- gsea_all |> rowwise() |> mutate(
      in_le = any(av$EntrezGeneID[toupper(av$EntrezGeneSymbol) %in% genes] %in%
                    parse_leading_edge(leadingEdge))
    ) |> filter(in_le) |> pull(pathway) |> head(3)
    if (length(gsea_le) > 0)
      lines <- c(lines, paste0("GSEA leading edge:\n  ",
                               paste(gsub("HALLMARK_|KEGG_", "", gsea_le) |>
                                       tolower() |> gsub("_", " ", x = _),
                                     collapse = "\n  ")))
    if (length(lines) == 0) return("(no other analysis hits)")
    paste(lines, collapse = "\n")
  })

  # ── A-V Volcano ─────────────────────────────────────────────────────────────
  output$volcano_plot <- renderPlotly({
    df <- av_marked()
    if (input$sig_filter == "nom") df <- filter(df, sig_nom)
    if (input$sig_filter == "fdr") df <- filter(df, sig_fdr)
    if (input$dir_filter == "Arterial > Venous") df <- filter(df, direction == "Arterial > Venous")
    if (input$dir_filter == "Venous > Arterial") df <- filter(df, direction == "Venous > Arterial")

    if (nrow(df) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers") |>
        layout(title = "No proteins pass the current filters"))
    }

    df <- df |>
      mutate(
        pt_col = case_when(
          is_hit  ~ "Searched",
          sig_fdr ~ "FDR < 0.05",
          sig_nom ~ "p < 0.05",
          TRUE    ~ "NS"
        ),
        pt_size  = case_when(is_hit ~ 9L, sig_fdr ~ 6L, sig_nom ~ 5L, TRUE ~ 3L),
        pt_alpha = if_else(pt_col == "NS", 0.25, 0.82)
      ) |>
      arrange(pt_col == "NS")

    col_map <- c("Searched" = "#e31a1c", "FDR < 0.05" = "#ff7f00",
                 "p < 0.05" = "#1f78b4", "NS" = "#bbbbbb")

    p <- plot_ly(df,
      x = ~logFC, y = ~neg_log10p,
      color = ~pt_col, colors = col_map,
      type = "scatter", mode = "markers",
      text = ~paste0(
        "<b>", Target, "</b> (", EntrezGeneSymbol, ")<br>",
        "logFC: ", round(logFC, 3), "<br>",
        "p = ", signif(P.Value, 3), " | adj.p = ", signif(adj.P.Val, 3), "<br>",
        direction
      ),
      hoverinfo = "text",
      marker = list(opacity = ~pt_alpha, size = ~pt_size)
    ) |>
      add_segments(x = 0, xend = 0, y = 0, yend = max(df$neg_log10p, na.rm = TRUE),
        line = list(color = "gray70", dash = "dot", width = 1),
        showlegend = FALSE, hoverinfo = "none", inherit = FALSE) |>
      add_segments(
        x = min(df$logFC, na.rm = TRUE), xend = max(df$logFC, na.rm = TRUE),
        y = -log10(0.05), yend = -log10(0.05),
        line = list(color = "#1f78b4", dash = "dot", width = 1),
        showlegend = FALSE, hoverinfo = "none", inherit = FALSE) |>
      layout(
        xaxis = list(title = "log2 Fold Change (Arterial / Venous)", zeroline = FALSE),
        yaxis = list(title = "−log10(p-value)"),
        legend = list(title = list(text = "Significance")),
        hovermode = "closest"
      )

    hits <- filter(df, is_hit)
    if (nrow(hits) > 0) {
      p <- p |> add_annotations(
        x = hits$logFC, y = hits$neg_log10p, text = hits$EntrezGeneSymbol,
        showarrow = TRUE, arrowhead = 2, arrowsize = 0.6, arrowwidth = 1.2,
        font = list(size = 12, color = "#c0392b"),
        bgcolor = "white", bordercolor = "#c0392b", borderwidth = 1.2,
        borderpad = 3, opacity = 0.9
      )
    }
    p
  })

  output$av_table <- renderReactable({
    genes <- searched_genes()
    df    <- av_marked()
    if (length(genes) > 0) {
      df <- filter(df, is_hit)
    } else {
      df <- filter(df, sig_nom) |> arrange(P.Value) |> head(300)
    }
    df |>
      mutate(logFC = round(logFC, 3),
             P.Value = signif(P.Value, 3),
             adj.P.Val = signif(adj.P.Val, 3)) |>
      select(Gene = EntrezGeneSymbol, Protein = Target,
             logFC, `p` = P.Value, `adj.p` = adj.P.Val,
             Direction = direction, Nominal = sig_nom, FDR = sig_fdr) |>
      rtbl(defaultSorted = "p",
           columns = list(
             Nominal = colDef(cell = function(v) if (v) "✓" else ""),
             FDR     = colDef(cell = function(v) if (v) "✓" else "")
           ))
  })

  # ── PCA ─────────────────────────────────────────────────────────────────────
  output$pca_loadings_plot <- renderPlotly({
    genes <- searched_genes()
    pct   <- pca_var$var_pct

    df <- loadings |>
      mutate(is_hit = toupper(EntrezGeneSymbol) %in% genes,
             col    = if_else(is_hit, "Searched", "Other")) |>
      arrange(col == "Other")

    p <- plot_ly(df, x = ~PC1, y = ~PC2,
      color = ~col, colors = c("Searched" = "#e31a1c", "Other" = "#aaaaaa"),
      type = "scatter", mode = "markers",
      text = ~paste0("<b>", Target, "</b> (", EntrezGeneSymbol, ")<br>",
                     "PC1: ", round(PC1, 4), "  PC2: ", round(PC2, 4)),
      hoverinfo = "text",
      marker = list(size = ~if_else(is_hit, 8, 4),
                    opacity = ~if_else(is_hit, 1, 0.4))
    ) |>
      layout(xaxis = list(title = paste0("PC1 (", pct[1], "%)")),
             yaxis = list(title = paste0("PC2 (", pct[2], "%)")))

    hits <- filter(df, is_hit)
    if (nrow(hits) > 0) {
      p <- p |> add_annotations(
        x = hits$PC1, y = hits$PC2, text = hits$EntrezGeneSymbol,
        showarrow = TRUE, arrowhead = 2,
        font = list(size = 11, color = "#c0392b"),
        bgcolor = "white", bordercolor = "#c0392b"
      )
    }
    p
  })

  output$pca_top_table <- renderReactable({
    loadings |>
      mutate(abs_PC1 = abs(PC1)) |>
      arrange(desc(abs_PC1)) |>
      mutate(across(c(PC1, PC2, PC3, PC4), ~round(.x, 4))) |>
      select(Gene = EntrezGeneSymbol, Protein = Target, PC1, PC2, PC3, PC4) |>
      head(100) |>
      rtbl(defaultSorted = list(PC1 = "desc"))
  })

  if (!is.null(pca_scores)) {
    output$pca_scores_plot <- renderPlotly({
      # Pair the A and V rows per patient for line segments
      art <- filter(pca_scores, draw == "Arterial")
      ven <- filter(pca_scores, draw == "Venous")
      paired <- inner_join(art, ven, by = "patient", suffix = c("_a", "_v"))

      p <- plot_ly(pca_scores,
        x = ~PC1, y = ~PC2,
        color = ~draw, colors = c("Arterial" = "#d73027", "Venous" = "#4575b4"),
        symbol = ~surgery_group,
        symbols = c("Head/Neck" = "triangle-up", "Laparoscopic" = "square", "Spine" = "circle"),
        type = "scatter", mode = "markers+text",
        text = ~patient, textposition = "top center",
        marker = list(size = 11, opacity = 0.9)
      )
      for (i in seq_len(nrow(paired))) {
        p <- add_segments(p,
          x = paired$PC1_a[i], xend = paired$PC1_v[i],
          y = paired$PC2_a[i], yend = paired$PC2_v[i],
          line = list(color = "gray70", width = 1),
          showlegend = FALSE, hoverinfo = "none", inherit = FALSE
        )
      }
      p |> layout(
        xaxis = list(title = paste0("PC1 (", pca_var$var_pct[1], "%)")),
        yaxis = list(title = paste0("PC2 (", pca_var$var_pct[2], "%)"))
      )
    })
  }

  # ── GSEA ────────────────────────────────────────────────────────────────────
  gsea_reactive <- reactive({
    df <- gsea_all
    if (input$gsea_collection != "All") df <- filter(df, collection == input$gsea_collection)
    if (input$gsea_dir == "neg") df <- filter(df, NES < 0)
    if (input$gsea_dir == "pos") df <- filter(df, NES > 0)
    filter(df, padj <= input$gsea_padj) |> arrange(padj)
  })

  clean_pathway <- function(x) {
    gsub("HALLMARK_|KEGG_", "", x) |> tolower() |> gsub("_", " ", x = _)
  }

  output$gsea_table <- renderReactable({
    gsea_reactive() |>
      mutate(pathway = clean_pathway(pathway),
             NES  = round(NES, 3),
             pval = signif(pval, 3),
             padj = signif(padj, 3)) |>
      select(Collection = collection, Pathway = pathway,
             NES, `p-val` = pval, `adj.p` = padj, Size = size) |>
      rtbl(selection = "single", onClick = "select",
           defaultSorted = "adj.p")
  })

  selected_pathway_row <- reactive({
    sel <- getReactableState("gsea_table", "selected")
    if (is.null(sel)) return(NULL)
    gsea_reactive()[sel, ]
  })

  leading_edge_symbols <- reactive({
    pw <- selected_pathway_row()
    if (is.null(pw)) return(character(0))
    le_ids <- parse_leading_edge(pw$leadingEdge)
    av |> filter(EntrezGeneID %in% le_ids) |> arrange(P.Value) |> pull(EntrezGeneSymbol)
  })

  output$leading_edge_genes <- renderText({
    syms <- leading_edge_symbols()
    if (length(syms) == 0) return("(click a pathway row)")
    paste(syms, collapse = ", ")
  })

  output$gsea_volcano <- renderPlotly({
    le_genes <- leading_edge_symbols()
    searched <- searched_genes()
    pw       <- selected_pathway_row()

    df <- av |>
      mutate(
        in_le     = toupper(EntrezGeneSymbol) %in% toupper(le_genes),
        is_search = toupper(EntrezGeneSymbol) %in% searched,
        col = case_when(
          is_search ~ "Searched",
          in_le     ~ "Leading edge",
          sig_nom   ~ "p < 0.05",
          TRUE      ~ "NS"
        )
      ) |>
      arrange(col == "NS")

    col_map <- c("Searched" = "#e31a1c", "Leading edge" = "#ff7f00",
                 "p < 0.05" = "#1f78b4", "NS" = "#dddddd")

    plot_ly(df, x = ~logFC, y = ~neg_log10p,
      color = ~col, colors = col_map,
      type = "scatter", mode = "markers",
      text = ~paste0("<b>", Target, "</b> (", EntrezGeneSymbol, ")<br>",
                     "logFC: ", round(logFC, 3), "  p=", signif(P.Value, 3)),
      hoverinfo = "text",
      marker = list(
        size    = ~if_else(col %in% c("Searched", "Leading edge"), 6, 3),
        opacity = ~if_else(col == "NS", 0.2, 0.75)
      )
    ) |>
      layout(
        xaxis = list(title = "log2FC (Arterial/Venous)"),
        yaxis = list(title = "−log10(p)"),
        title = if (!is.null(pw)) clean_pathway(pw$pathway) else "Select a pathway above"
      )
  })

  # ── Skin Proteases ──────────────────────────────────────────────────────────
  skin_filtered <- reactive({
    req(!is.null(skin_av))
    if (input$skin_family == "All") skin_av else filter(skin_av, family == input$skin_family)
  })

  output$skin_forest <- renderPlotly({
    req(!is.null(skin_av))
    df <- skin_filtered() |>
      mutate(
        label     = paste0(EntrezGeneSymbol,
                           if_else(!is.na(Target) & Target != "" & Target != EntrezGeneSymbol,
                                   paste0(" (", Target, ")"), "")),
        sig       = P.Value < 0.05,
        point_col = if_else(logFC > 0, "Arterial > Venous", "Venous > Arterial")
      ) |>
      arrange(logFC) |>
      mutate(label = factor(label, levels = label))

    col_map <- c("Arterial > Venous" = "#c0392b", "Venous > Arterial" = "#2980b9")

    plot_ly(df,
      x = ~logFC, y = ~label,
      color = ~point_col, colors = col_map,
      symbol = ~family, symbols = c(SPINK = "circle", KLK = "square", SERPIN = "diamond"),
      type = "scatter", mode = "markers",
      text = ~paste0(
        "<b>", EntrezGeneSymbol, "</b> (", family, ")<br>",
        if_else(!is.na(Target) & Target != "", paste0(Target, "<br>"), ""),
        "logFC = ", round(logFC, 3), "<br>",
        "p = ", signif(P.Value, 3), " | ", direction, "<br>",
        skin_role
      ),
      hoverinfo = "text",
      marker = list(
        size    = ~if_else(sig, 13, 7),
        opacity = ~if_else(sig, 0.9, 0.45),
        line    = list(color = "white", width = 0.5)
      )
    ) |>
      add_segments(
        x = 0, xend = 0, y = 0.5, yend = nrow(df) + 0.5,
        line = list(color = "gray55", dash = "dot", width = 1),
        showlegend = FALSE, hoverinfo = "none", inherit = FALSE
      ) |>
      layout(
        xaxis = list(title = "log2FC (Arterial / Venous)   ← venous higher | arterial higher →",
                     zeroline = FALSE),
        yaxis = list(title = "", tickfont = list(size = 10)),
        legend = list(title = list(text = "Direction / Shape = Family")),
        margin = list(l = 180)
      )
  })

  output$skin_table <- renderReactable({
    req(!is.null(skin_av))
    skin_filtered() |>
      mutate(logFC     = round(logFC, 3),
             P.Value   = signif(P.Value, 3),
             sig       = P.Value < 0.05) |>
      arrange(P.Value) |>
      select(Family = family, Gene = EntrezGeneSymbol,
             Protein = Target, logFC, p = P.Value,
             Direction = direction, `p<0.05` = sig) |>
      rtbl(
        defaultSorted = "p",
        columns = list(
          logFC = colDef(
            style = function(val) {
              list(color = if (!is.na(val) && val < 0) "#2980b9" else "#c0392b",
                   fontWeight = "bold")
            }
          ),
          `p<0.05` = colDef(cell = function(v) if (isTRUE(v)) "✓" else "")
        )
      )
  })

  if (!is.null(skin_kw)) {
    output$skin_kw_table <- renderReactable({
      skin_kw |>
        mutate(across(where(is.numeric), \(x) signif(x, 3))) |>
        select(Gene = gene_symbol,
               `A-V p`   = av_p,
               `KW p`    = kw_p,
               `Δ Head/Neck` = mean_HN,
               `Δ Lap`   = mean_Lap,
               `Δ Spine` = mean_Spine) |>
        rtbl(defaultSorted = "A-V p")
    })
  }

  # ── ECS / EOS ───────────────────────────────────────────────────────────────
  ecs_filtered <- reactive({
    df <- ecs_eos
    if (input$ecs_system != "both")  df <- filter(df, system == input$ecs_system)
    if (input$ecs_platform != "all") df <- filter(df, platform == input$ecs_platform)
    df
  })

  output$ecs_table <- renderReactable({
    ecs_filtered() |>
      mutate(logFC = round(logFC, 3), p_value = signif(p_value, 3)) |>
      select(System = system, Analysis = analysis, Platform = platform,
             Gene = gene_symbol, Category = category,
             logFC, p = p_value, Direction = direction) |>
      rtbl(selection = "single", onClick = "select",
           defaultSorted = "p")
  })

  selected_ecs_gene <- reactive({
    sel <- getReactableState("ecs_table", "selected")
    if (is.null(sel)) return(character(0))
    toupper(ecs_filtered()$gene_symbol[sel])
  })

  output$ecs_volcano <- renderPlotly({
    all_ecs_genes <- toupper(ecs_eos$gene_symbol)
    searched      <- searched_genes()
    sel_gene      <- selected_ecs_gene()

    df <- av |>
      mutate(
        sym_up    = toupper(EntrezGeneSymbol),
        is_sel    = sym_up %in% sel_gene,
        is_search = sym_up %in% searched & !is_sel,
        is_ecs    = sym_up %in% all_ecs_genes & !is_sel & !is_search,
        col = case_when(
          is_sel    ~ "Selected",
          is_search ~ "Searched",
          is_ecs    ~ "ECS/EOS",
          sig_nom   ~ "p < 0.05",
          TRUE      ~ "NS"
        )
      ) |>
      arrange(col == "NS")

    col_map <- c("Selected" = "#e31a1c", "Searched" = "#9400d3",
                 "ECS/EOS" = "#ff7f00", "p < 0.05" = "#1f78b4", "NS" = "#dddddd")

    p <- plot_ly(df, x = ~logFC, y = ~neg_log10p,
      color = ~col, colors = col_map,
      type = "scatter", mode = "markers",
      text = ~paste0("<b>", Target, "</b> (", EntrezGeneSymbol, ")<br>",
                     "logFC: ", round(logFC, 3), "  p=", signif(P.Value, 3)),
      hoverinfo = "text",
      marker = list(
        size    = ~case_when(is_sel ~ 10, is_ecs | is_search ~ 7, TRUE ~ 3),
        opacity = ~if_else(col == "NS", 0.2, 0.8)
      )
    ) |>
      layout(
        xaxis = list(title = "log2FC (Arterial/Venous)"),
        yaxis = list(title = "−log10(p)"),
        title = "ECS / EOS gene positions on A-V volcano"
      )

    hits <- filter(df, is_sel | is_search)
    if (nrow(hits) > 0) {
      p <- p |> add_annotations(
        x = hits$logFC, y = hits$neg_log10p, text = hits$EntrezGeneSymbol,
        showarrow = TRUE, arrowhead = 2, arrowsize = 0.6,
        font = list(size = 11,
                    color = if_else(hits$is_sel, "#c0392b", "#6a0dad")),
        bgcolor = "white",
        bordercolor = if_else(hits$is_sel, "#c0392b", "#6a0dad")
      )
    }
    p
  })
}

shinyApp(ui, server)
