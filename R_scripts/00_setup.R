# 00_setup.R
# Shared data loading and preparation sourced by all analysis scripts.
# Run any script standalone; it will source this file first.

suppressPackageStartupMessages({
  library(SomaDataIO)
  library(limma)
  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# ── Directory paths ────────────────────────────────────────────────────────────
if (!exists("data_dir")) {
  data_dir <- tryCatch(
    normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
    error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
  )
}
plot_dir  <- file.path(data_dir, "outputs", "plots")
table_dir <- file.path(data_dir, "outputs", "tables")
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ── Colour palettes ────────────────────────────────────────────────────────────
col_av   <- c("Arterial" = "#c0392b", "Venous" = "#2980b9")
col_yn   <- c("Yes" = "#c0392b", "No" = "#2980b9")
col_surg <- setNames(RColorBrewer::brewer.pal(3, "Set2"),
                     c("Head/Neck", "Laparoscopic", "Spine"))

# ── Surgery map ────────────────────────────────────────────────────────────────
surgery_map <- tribble(
  ~patient, ~surgery_type,            ~surgery_group,
  "0179",   "Spine",                  "Spine",
  "0193",   "Lap. cholecystectomy",   "Laparoscopic",
  "0209",   "Spine",                  "Spine",
  "0215",   "Spine",                  "Spine",
  "0223",   "Craniectomy",            "Head/Neck",
  "0245",   "Lap. adrenalectomy",     "Laparoscopic",
  "0250",   "Lap. nephrectomy",       "Laparoscopic",
  "0254",   "Carotid endarterectomy", "Head/Neck",
  "0268",   "Thyroidectomy",          "Head/Neck",
  "0275",   "Transsphenoidal",        "Head/Neck",
  "0284",   "Spine",                  "Spine",
  "0291",   "Spine",                  "Spine"
)

# ── Luminex panel gene map (used in GSEA and platform concordance) ─────────────
luminex_gene_map <- tribble(
  ~luminex_name,           ~gene_symbol,
  "MCP-2(CCL8)",           "CCL8",
  "IL-2R",                 "IL2RA",
  "MIP-1a(CCL3)",          "CCL3",
  "SDF-1a",                "CXCL12",
  "IL-27",                 "IL27",
  "LIF",                   "LIF",
  "IL-1b",                 "IL1B",
  "IL-2",                  "IL2",
  "IL-4",                  "IL4",
  "IL-5",                  "IL5",
  "IP-10(CXCL10)",         "CXCL10",
  "IL-6",                  "IL6",
  "IL-7",                  "IL7",
  "IL-8(CXCL8)",           "CXCL8",
  "IL-10",                 "IL10",
  "BLC(CXCL13)",           "CXCL13",
  "EOTAXIN-2(CCL24)",      "CCL24",
  "EOTAXIN(CCL11)",        "CCL11",
  "IL-12p70",              "IL12A",
  "IL-13",                 "IL13",
  "IL-17A(CTLA-8)",        "IL17A",
  "IL-31",                 "IL31",
  "SCF",                   "KITLG",
  "G-CSF(CSF-3)",          "CSF3",
  "IFNg",                  "IFNG",
  "GM-CSF",                "CSF2",
  "TNFa",                  "TNF",
  "HGF",                   "HGF",
  "MIP-1b(CCL4)",          "CCL4",
  "IFNa",                  "IFNA1",
  "EOTAXIN-3(CCL26)",      "CCL26",
  "MCP-1(CCL2)",           "CCL2",
  "IL-9",                  "IL9",
  "MIF",                   "MIF",
  "TNFb",                  "LTA",
  "NGFb",                  "NGF",
  "MIP-3a(CCL20)",         "CCL20",
  "I-TAC(CXCL11)",         "CXCL11",
  "TRAIL",                 "TNFSF10",
  "Fractalkine(CX3CL1)",   "CX3CL1",
  "GROa(CXCL1)",           "CXCL1",
  "IL-1a",                 "IL1A",
  "IL-23",                 "IL23A",
  "MMP-1",                 "MMP1",
  "IL-15",                 "IL15",
  "IL-18",                 "IL18",
  "M-CSF",                 "CSF1",
  "MCP-3(CCL7)",           "CCL7",
  "MIG(CXCL9)",            "CXCL9",
  "IL-16",                 "IL16",
  "IL-21",                 "IL21",
  "IL-3",                  "IL3",
  "CD40L",                 "CD40LG",
  "FGF-2",                 "FGF2",
  "IL-22",                 "IL22",
  "VEGF-A",                "VEGFA",
  "TSLP",                  "TSLP",
  "IL-20",                 "IL20",
  "ENA-78(CXCL5)",         "CXCL5",
  "CD30",                  "TNFRSF8",
  "TNF-RII",               "TNFRSF1B",
  "BAFF",                  "TNFSF13B",
  "MDC",                   "CCL22",
  "APRIL",                 "TNFSF13",
  "TWEAK",                 "TNFSF12"
)

# ── Load SomaScan ──────────────────────────────────────────────────────────────
adat_file <- file.path(
  data_dir,
  "SS-229807_v4.1_EDTA_Plasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat"
)
soma         <- read_adat(adat_file)
analyte_info <- getAnalyteInfo(soma)
pass_analytes <- analyte_info |> dplyr::filter(ColCheck == "PASS") |> pull(AptName)

soma_samples <- soma |>
  dplyr::filter(SampleType == "Sample") |>
  mutate(
    patient = str_extract(SampleId, "\\d{4}"),
    draw    = if_else(str_ends(SampleId, "-2"), "Venous", "Arterial"),
    flagged = RowCheck == "FLAG"
  ) |>
  left_join(surgery_map |> select(patient, surgery_group), by = "patient")

# ── Covariates ─────────────────────────────────────────────────────────────────
clean_yn <- function(x) {
  case_when(str_to_upper(str_trim(x)) == "Y" ~ TRUE,
            str_to_upper(str_trim(x)) == "N" ~ FALSE,
            TRUE ~ NA)
}
covariate_raw <- read_xlsx(
  file.path(data_dir,
    "SS-229807_v4.1_EDTA_Plasma - Art vs Venous and surgery type and other covariates after AP2 PK copy.xlsx")
) |>
  rename(sample_id = SampleId,
         obese         = `Obese (BMI>30)`,
         hypothermia   = `Hypothermia (Temp <36C)`,
         hyperglycemia = `Hyperglycemia (BG>150)`,
         male          = `Male (Y/N)`) |>
  mutate(patient = str_extract(sample_id, "\\d{4}"),
         across(c(obese, hypothermia, hyperglycemia, male), clean_yn)) |>
  group_by(patient) |>
  summarise(across(c(obese, hypothermia, hyperglycemia, male), first), .groups = "drop")

soma_samples <- soma_samples |> left_join(covariate_raw, by = "patient")

# ── log2-transformed expression matrix ────────────────────────────────────────
expr_all <- soma_samples |>
  select(all_of(pass_analytes)) |>
  as.matrix() |>
  log2()
rownames(expr_all) <- soma_samples$SampleId

binary_covariates <- c("obese", "hypothermia", "hyperglycemia", "male")

cat("Setup complete —", nrow(soma_samples), "samples,", length(pass_analytes), "PASS analytes\n")
