# Sarcopenia–Obesity

This repository contains R scripts to generate figures and analysis outputs for a **Sarcopenia–Obesity (SO)** project, including microbiome community typing (enterotypes), alpha/beta diversity, a top-feature heatmap, and a multiclass XGBoost model with OOF (out-of-fold) evaluation and SHAP importance.

## Contents

- `Fig1.abcsup.R`  
  **Fig1a** Enterotype (PAM clustering on Bray–Curtis) + PCoA  
  **Fig1b** Top genera boxplots (pseudo-log scale) + pairwise tests  
  **Fig1c** Enterotype proportions across groups (SO / Sarcopenia / Obesity / Healthy)  
  Also exports multiple intermediate CSVs for supplements.

- `Fig2.alpha.R`  
  Alpha diversity metrics + group comparisons (Wilcoxon, BH/FDR) + a GLM “forest plot style” panel (uses filtered taxonomic features).

- `Fig3.PCoA.R`  
  Bray–Curtis PCoA + **confounder-adjusted PERMANOVA (adonis2)** + PERMDISP + multi-panel layout; exports design and test tables.

- `Fig5.heatmap.R`  
  Top-30 features selected by Kruskal–Wallis (FDR) → group means → row Z-score heatmap; adds “*” for each group vs Healthy (Wilcoxon + FDR < 0.05).

- `Fig6.XGBOOSTSHAP.R`  
  **4-class XGBoost** (SO / S / O / H) with **stratified 5-fold CV**, OOF predictions, confusion matrices, One-vs-Rest AUC, and **pairwise ROC + 95% CI** for **SO vs O/S/H**. Trains a final model on all samples (using CV-derived nrounds) and outputs feature importance + SHAP mean(|SHAP|).

---

## Requirements

- R >= 4.1 recommended (4.2/4.3+ preferred)
- Key packages used across scripts include:
  - Data I/O & wrangling: `data.table`, `readr`, `dplyr`, `tidyr`, `stringr`
  - Microbiome ecology: `vegan`, `ape`, `cluster`
  - Plotting: `ggplot2`, `scales`, `patchwork`, `ggpubr`, `rstatix`, `svglite`
  - Stats helpers: `broom`, `multcomp`
  - Heatmap: `ComplexHeatmap`, `circlize`, `grid`
  - ML: `xgboost`, `pROC`

> Note: `Fig1.abcsup.R` auto-installs missing packages. Others may require you to install packages manually.

---

## Data Inputs

### 1) Genus relative abundance (used by Fig1/Fig2/Fig3)
**File:** `trans.genus.relative.csv`  
**Format:** rows = samples, columns = genera; **first column = SampleID**.

### 2) Sample metadata / mapping (used by Fig2/Fig3; Fig1c uses grouping too)
**File:** `mappingV1.csv`  
Must include:
- `ID` (matching SampleID in abundance file)
- `group_code` with the coding:
  - `1 = SO`
  - `2 = Sarcopenia`
  - `3 = Obesity`
  - `4 = Healthy`

**Optional covariates (Fig3 confounder screening):**
- Numeric candidates: `age`, `bmi`, `mmse`, `phq9`, `SPPB`, `grip_strength`
- Categorical candidates: `sex`, `hearing_impaired`, `teeth_loss`, `constipation`, `yogurt`, `eatprob`, `diarrhea`, `laxatives`  
Sex is harmonized if coded as 1/2 → Male/Female.

### 3) Feature table for heatmap (Fig5)
**File:** `abund.csv`  
Must include a `SampleID` column and numeric feature columns.

Group labels are inferred from `SampleID` prefix:
- `SO...` → SO  
- `S...`  → Sarcopenia  
- `O...`  → Obesity  
- `H...`  → Healthy  

### 4) Feature matrix for XGBoost (Fig6)
**File:** `V1abund_plus_species_merged_by_SampleID.csv`  
Must include `SampleID` and numeric features.

Class label is parsed from `SampleID` and expects patterns like:
- `SO-xxx` or `SO_xxx`
- `S-xxx` or `S_xxx`
- `O-xxx` or `O_xxx`
- `H-xxx` or `H_xxx`

---

## Before Running (Important)

The scripts currently contain **absolute paths** such as `setwd("D:/...")` or `workdir <- "D:/..."`.

You must update:
- `setwd(...)` / `workdir`
- input file paths (`abund_file`, `mapping_file`, `infile`, etc.)

Recommended folder structure (optional but convenient):
