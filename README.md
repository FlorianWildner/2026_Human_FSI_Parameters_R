# Human FSI Electrophysiology & IMARIS Analysis

This repository contains the R code used for the analysis and visualization of electrophysiological and anatomical data in the study.

All data have been pseudonymized to ensure protection of personal information.

---

## Repository structure

data/
Human_IN_Properties.csv
WFA_HUMAN_QUANT_IMARIS_PatA.csv
WFA_HUMAN_QUANT_IMARIS_PatB.csv
WFA_HUMAN_QUANT_IMARIS_PatC.csv

FSIN_intrinsic.Rmd
WFA_IMARIS.R


---

## Electrophysiological analysis

**Script:** `FSIN_intrinsic.Rmd`

This script:
- Loads the pseudonymized dataset `data/Human_IN_Properties.csv`
- Performs statistical analysis of intrinsic electrophysiological properties
- Generates all analyses and figures related to electrophysiological data presented in the main manuscript

---

## Anatomical / IMARIS analysis

**Script:** `WFA_IMARIS.R`

This script:
- Loads IMARIS-derived anatomical datasets:
  - `PatA`
  - `PatB`
  - `PatC`
- Performs preprocessing and statistical analysis of WFA staining data
- Generates analyses and figures corresponding to:
  - Figure 4
  - Supplementary Figure S4

---

## Requirements

The code was developed in R and uses common statistical and visualization packages (e.g., `tidyverse`, `lme4`, `ggplot2`, `emmeans`).

---

## Data availability

All data included in this repository are pseudonymized to ensure privacy and compliance with data protection regulations.

---

## Notes

- Paths are defined relative to the project structure to ensure reproducibility.
- Figures are automatically saved to an output directory when scripts are executed.
