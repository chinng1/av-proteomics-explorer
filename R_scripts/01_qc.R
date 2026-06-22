# 01_qc.R
# SomaScan sample and analyte QC flags.
# Outputs: outputs/tables/01_sample_qc.csv, 01_analyte_qc.csv

data_dir <- tryCatch(
  normalizePath(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")),
  error = function(e) "/Users/gregchinn/Desktop/coding projects/Re_ Metabolomics"
)
source(file.path(data_dir, "R_scripts", "00_setup.R"))

# ── Analyte QC ─────────────────────────────────────────────────────────────────
analyte_qc <- analyte_info |>
  select(AptName, Target, EntrezGeneSymbol, ColCheck, UniProt) |>
  arrange(ColCheck, Target)

cat("Analyte ColCheck summary:\n")
print(table(analyte_qc$ColCheck))

write_csv(analyte_qc, file.path(table_dir, "01_analyte_qc.csv"))

# ── Sample QC ──────────────────────────────────────────────────────────────────
sample_qc <- soma_samples |>
  select(SampleId, patient, draw, surgery_group, flagged,
         obese, hypothermia, hyperglycemia, male) |>
  arrange(patient, draw)

cat("\nSample RowCheck FLAG:\n")
print(sample_qc |> dplyr::filter(flagged) |> select(SampleId, patient, draw))

write_csv(sample_qc, file.path(table_dir, "01_sample_qc.csv"))

cat("\nOutputs saved to outputs/tables/\n")
