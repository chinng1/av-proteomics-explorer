# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**PI:** Arun Prakash, M.D. | **Co-investigator:** Philip Kurien, M.D.  
**Analyst:** Gregory A. Chinn, M.D. Ph.D.  
**Question:** Do plasma protein levels differ between simultaneously collected arterial and venous blood during elective surgery?  
**Cohort:** 12 patients; one simultaneous A+V draw each at a single intraoperative time point.

The arteriovenous (A-V) gradient is the unit of inference: arterial > venous = peripheral tissue *consuming* the mediator; venous > arterial = tissue *releasing* it.

## Running the Analysis

**Primary entry point — modular R pipeline** (preferred for reproduction or extension):

```bash
Rscript R_scripts/00_setup.R   # must be sourced first; all others source it automatically
Rscript R_scripts/01_qc.R
Rscript R_scripts/02_pca.R
Rscript R_scripts/03_av_limma.R
Rscript R_scripts/04_gsea.R
Rscript R_scripts/05_covariate_analysis.R
Rscript R_scripts/06_av_delta.R
Rscript R_scripts/07_platform_concordance.R
Rscript R_scripts/08_endocannabinoid_analysis.R  # reads CSV outputs; no SomaDataIO needed
Rscript R_scripts/09_opioid_analysis.R           # requires SomaDataIO
Rscript R_scripts/10_skin_protease_analysis.R    # requires SomaDataIO; run after 03/05/06
```

**Monolithic Quarto documents** (self-contained; alternative to the pipeline):

```bash
quarto render arterial_venous_analysis.qmd   # SomaScan + Luminex analysis
quarto render somascan_analysis.qmd          # SomaScan only
quarto render summary_report.qmd            # narrative summary
```

**Shiny explorer** (requires pipeline outputs in `outputs/tables/`):

```r
shiny::runApp("shiny_explorer")  # run from project root
```

## Architecture

### Modular pipeline (`R_scripts/`)

| Script | Purpose | Key output |
|--------|---------|------------|
| `00_setup.R` | Shared data load; sourced automatically by all other scripts | `soma_samples`, `expr_all`, `covariate_raw` objects |
| `01_qc.R` | SomaScan ColCheck/RowCheck flags | `outputs/tables/QC_*.csv` |
| `02_pca.R` | PCA on 7,481 PASS proteins | `02_pca_loadings.csv`, `02_pca_variance.csv` |
| `03_av_limma.R` | Paired limma A vs V, blocking on patient | `03_av_results_full.csv` (used by Shiny) |
| `04_gsea.R` | GSEA on limma t-statistic (Hallmark + KEGG + Luminex custom set) | `04_gsea_hallmark.csv`, `04_gsea_kegg.csv` |
| `05_covariate_analysis.R` | Limma: venous-only, effects of obesity/hypothermia/hyperglycemia/sex | `Covariate_*.csv` |
| `06_av_delta.R` | Limma on per-patient A-V delta log2 by covariate | `AV_gradient_*.csv` |
| `07_platform_concordance.R` | Luminex vs. SomaScan fold-change comparison | `Platform_*.csv` |
| `08_endocannabinoid_analysis.R` | ECS protein sub-analysis from CSV outputs; no `.adat` needed | `08_*.csv/png` |
| `09_opioid_analysis.R` | Endogenous opioid system sub-analysis; requires `.adat` | `09_eos_*.csv/png` |
| `10_skin_protease_analysis.R` | KLK/SPINK skin protease cluster; top A-V finding; requires `.adat` | `10_skin_*.csv/png` |

`00_setup.R` exposes shared globals: `soma_samples`, `expr_all` (log2 matrix, rows = samples, cols = PASS aptamers), `pass_analytes`, `covariate_raw`, `surgery_map`, `luminex_gene_map`, `col_av`, `col_yn`, `col_surg`.

### Shiny explorer (`shiny_explorer/app.R`)

Reads pre-computed CSVs from `outputs/tables/` (does not re-run any analysis). Tabs: volcano, PCA, GSEA pathway browser, per-patient A-V dot plots. Dependencies: `shiny`, `bslib`, `plotly`, `reactable`.

## Data Files

| File | Contents |
|------|----------|
| `AP_65plex_byJMR.csv` | 65-plex Luminex; skip first 54 rows (instrument header) |
| `AP_Cytokines_12_2022.docx` | Plate layout map for Luminex run |
| `SS-229807_v4.1_...anmlSMP.adat` | SomaScan + ANML sample normalization — **use this for all differential analysis** |
| `SS-229807_v4.1_...anmlQC.qcCheck.adat` | SomaScan QC flags only |
| `SS-229807_v4.1_..._PK copy.xlsx` | Covariate table (surgery type, obesity, hypothermia, hyperglycemia, sex) |

## Sample ID Convention

```
000_XXXX-2  =  patient XXXX, Venous draw
000_XXXX-4  =  patient XXXX, Arterial draw
```

## Key Data Notes

- **Labeling flag:** Patient 0245's arterial Luminex well (F4) is transcribed as "245-2" in the plate map doc — likely should be "245-4". The `plate_map` tribble in `arterial_venous_analysis.qmd` carries an inline `# NOTE`. Verify against original tube labels before reporting 0245-specific results.
- **Power:** N=12 is underpowered for 65-analyte or ~7,500-protein BH correction. FDR failure reflects sample size, not absence of biology. Nominal p-values are the interpretively relevant threshold.
- `00_setup.R` sets `data_dir` via `rstudioapi::getActiveDocumentContext()$path` with a hard-coded fallback. If moving the directory, update the fallback path in `00_setup.R` (and in the `.qmd` files' `setup` chunks).
- SomaScan `.adat` file has all normalization already applied — do not re-normalize.
- GSEA uses `set.seed(42)`; all other scripts are deterministic.

## R Conventions

- Luminex values are MFI (positive scale); use log2FC = `log2((mean_A + 0.1) / (mean_V + 0.1))` with pseudocount.
- SomaScan values are log2-transformed before limma.
- `surgery_group` collapses individual surgery types into three categories (Spine, Laparoscopic, Head/Neck) via `surgery_map` tribble in `00_setup.R`.
- Color palettes are defined in `00_setup.R`: `col_av` (Arterial = red, Venous = blue), `col_yn`, `col_surg`.

## Dependencies

```r
install.packages(c(
  "tidyverse", "readxl", "ggplot2", "ggrepel",
  "pheatmap", "RColorBrewer", "limma", "fgsea", "msigdbr"
))
install.packages("SomaDataIO")          # required for scripts 00–07, 09
install.packages(c("shiny", "bslib", "plotly", "reactable"))  # Shiny explorer only
```
