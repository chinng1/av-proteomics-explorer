# Arterial vs. Venous Perioperative Proteomics — Analysis Guide

**Study:** Arteriovenous gradient of plasma proteins during elective surgery  
**PI:** Arun Prakash, M.D. | **Co-investigator:** Philip Kurien, M.D.  
**Analyst:** Gregory A. Chinn, M.D. Ph.D., UCSF Anesthesia & Perioperative Care

---

## What this analysis does

12 patients undergoing elective surgery each had one simultaneous arterial and venous blood draw intraoperatively. Two platforms measured proteins in these samples:

- **SomaScan v4.1** — aptamer-based proteomics, ~7,500 proteins
- **Luminex 65-plex** — antibody-based cytokine panel

The central question is whether protein concentrations differ between arterial and venous blood — and if so, in which direction. A **venous > arterial** gradient means peripheral tissues are *releasing* that protein into the venous circulation. **Arterial > venous** means tissues are *consuming* it.

---

## Data files required to reproduce the analysis

| File | Description |
|------|-------------|
| `AP_65plex_byJMR.csv` | Luminex raw export (54-row instrument header; data starts row 55) |
| `SS-229807_v4.1_EDTA_Plasma...anmlSMP.adat` | SomaScan — use this normalized file for all analyses |
| `SS-229807_v4.1_EDTA_Plasma...anmlQC.qcCheck.adat` | SomaScan QC-only file (used for QC flags, not differential analysis) |
| `SS-229807_v4.1_..._PK copy.xlsx` | Covariate table (surgery type, obesity, hypothermia, hyperglycemia, sex) |

---

## R package dependencies

```r
install.packages(c(
  "tidyverse", "readxl", "ggplot2", "ggrepel",
  "pheatmap", "RColorBrewer", "limma", "fgsea", "msigdbr"
))

# SomaScan file reader — required for scripts 01–07
install.packages("SomaDataIO")
```

---

## Running the pipeline

Scripts in `R_Scripts/` are numbered and run in order. Each reads inputs and writes
outputs to `outputs/tables/` and `outputs/plots/`. Run from the project root directory
(the folder containing this file and the data files).

```r
# In RStudio or Positron: open each script and source it, or run from terminal:
Rscript R_Scripts/00_setup.R
Rscript R_Scripts/01_qc.R
Rscript R_Scripts/02_pca.R
Rscript R_Scripts/03_av_limma.R
Rscript R_Scripts/04_gsea.R
Rscript R_Scripts/05_covariate_analysis.R
Rscript R_Scripts/06_av_delta.R
Rscript R_Scripts/07_platform_concordance.R
```

**Important:** Each script hard-codes the project directory path. If you move the
folder, update `data_dir` at the top of each script.

---

## Script overview

| Script | What it does | Key output |
|--------|-------------|------------|
| `00_setup.R` | Loads libraries, sets paths, defines color palettes | — |
| `01_qc.R` | SomaScan aptamer (ColCheck) and sample (RowCheck) QC flags | `QC_*.csv` |
| `02_pca.R` | PCA on all 7,481 PASS proteins; tests whether draw type, surgery, or covariates drive variation | `PCA_*.csv/png` |
| `03_av_limma.R` | Paired limma linear model: arterial vs. venous, blocking on patient | `SomaScan_AV_*.csv`, volcano, MA, heatmap |
| `04_gsea.R` | Gene-set enrichment on the A-V limma t-statistic ranking (Hallmark + KEGG + custom 65-plex set) | `GSEA_*.csv/png` |
| `05_covariate_analysis.R` | Limma: protein levels in venous samples by obesity, hypothermia, hyperglycemia, sex | `Covariate_*.csv/png` |
| `06_av_delta.R` | Limma on per-patient A-V delta (log2 arterial − log2 venous): does gradient magnitude differ by covariate? | `AV_gradient_*.csv/png` |
| `07_platform_concordance.R` | Cross-platform comparison: Luminex A-V fold-changes vs. SomaScan; per-patient Spearman correlations | `Platform_*.csv/png` |

---

## Key findings (summary)

- **Surgery type** is the dominant source of proteomic variation (PC1, PC2). Arterial vs. venous separation only emerges on PC3/PC4 — the A-V gradient is real but small relative to between-patient variation.
- **556 of 7,481 proteins** show a nominal A-V difference (p < 0.05); none survive Benjamini-Hochberg FDR correction at N = 12. This reflects sample size, not absence of biology.
- **Overall directional bias:** 62% of nominal hits are arterial > venous (binomial p = 1.5 × 10⁻⁸), suggesting the peripheral tissue bed is broadly consuming circulating proteins intraoperatively.
- **GSEA main finding:** Inflammatory response, IL-6/JAK-STAT3, coagulation, and complement pathways are enriched in the venous return (tissue releasing). PI3K-AKT-mTOR and protein secretion pathways are enriched in arterial blood.
- **Strongest single-protein finding:** IP-10/CXCL10 — consistently venous > arterial on both Luminex (p = 0.03) and SomaScan, suggesting active release of this IFN-γ–driven chemokine from the peripheral surgical bed.
- **Covariate modifier:** Hypothermia has the largest effect on the A-V gradient — hypothermic patients show broader and more extreme A-V differences across many proteins.

---

## Using Claude Code to extend the analysis

This project was analyzed using [Claude Code](https://claude.com/claude-code), Anthropic's
AI coding assistant. If you want to run additional analyses or ask questions about the data,
Claude Code can read R scripts, interpret results, and write new analyses directly.

### Getting started

1. Install Claude Code: `npm install -g @anthropic-ai/claude-code`
2. Open a terminal in the project folder and run: `claude`
3. Claude Code will read the project files and can answer questions or run new analyses

### Example questions you can ask

```
"Which proteins show the largest arterial-venous difference in spine surgery patients?"

"Are there any proteins associated with obesity that also show an A-V gradient?"

"Run a sensitivity analysis excluding the four flagged samples and compare to the main results."

"Plot the A-V gradient for [protein name] across all 12 patients."

"What do the GSEA results mean biologically for perioperative inflammation?"
```

### How Claude Code works with this project

- It reads the numbered R scripts to understand the analysis pipeline
- It reads the output CSVs in `outputs/tables/` to work with results without re-running the full pipeline
- It can write new R scripts and run them directly
- For questions about SomaScan data specifically, it will load the `.adat` file using the `SomaDataIO` package

### Notes for reproducibility

- All scripts write outputs deterministically except GSEA (uses `set.seed(42)` for permutation testing)
- The SomaScan `.adat` file contains all normalization steps already applied — do not re-normalize
- Luminex concentrations are in MFI (median fluorescence intensity), not calibrated concentrations
- Patient 0245's arterial Luminex well was labeled "245-2" in the plate map (likely a transcription error for "245-4") — verify against original tube records before reporting 0245-specific results

---

*Analysis completed June 2026. Contact: gregorychinn@gmail.com*
