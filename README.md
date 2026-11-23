# MCPH_Figures
This repository contains the complete R workflow and datasets required to reproduce the visual analyses of the MCPH-associated genes study. The pipeline generates publication-ready figures covering functional pathway involvement, structural brain trait associations, evolutionary metrics (dN/dS, adaptation rates), and High-Frequency Variant (HFV) profiling.
## 🌌 Visual Abstract
<img width="2000" height="1600" alt="Figures_Preview" src="https://github.com/user-attachments/assets/5c4b3c04-8d70-430c-8744-5988aec86da8" />

## 🧬 Repository Structure
* **`Data/`**: Contains the input datasets for each generated figure.
* **`Scripts/`**: R scripts used for data pre-processing, categorization, and plotting.
* **`Figures/`**: Final exported high-resolution figures (PDF for vector formats and PNG for web preview).
## 🧠 Dependencies
To reproduce the figures, you will need R version 4.6.2 installed along with the following packages:
### Figure 1
* `ggplot2`, `dplyr`, `readxl`, `forcats`, `ggh4x`, `ragg`, `showtext`.
### Figure 2
* `ggplot2`, `dplyr`, `readxl`, `ragg`, `showtext`, `patchwork`.
### Figure 3
* `ggplot2`, `dplyr`, `readxl`, `ragg`, `showtext`, `patchwork`, `ggExtra`, `ggrepel`.
### Figure 4
* `ggplot2`, `dplyr`, `readxl`, `ragg`, `showtext`, `patchwork`.
  
## 🔮 Figures
### Figure 1: MCPH Gene Functions (Binary Dot Plot) (`Script_fig1.R`)
This script creates a binary dot plot to visualize which specific functional pathways each MCPH-associated gene is involved in, dynamically grouping them into two main blocks ("Cell Division-Related Processes" and "Other Functions & Pathways").
<div align="center">
  <img src="https://github.com/user-attachments/assets/ceb2a311-947c-425c-9292-af93703fb328" width="600">
</div>

### Figure 2: HFV Density & Total Counts (`Script_fig2.R`)
This script generates a two-panel slim bar plot comparing the amino acid-normalized density (`HFVs / aa (%)`) and the raw total count (`Total HFVs`) of High-Frequency Variants across all MCPH-associated genes, dynamically sorted by descending variant density.
<div align="center">
  <img src="https://github.com/user-attachments/assets/34759fd2-b483-4918-a4c9-f255598e5dba" width="600">
</div>

### Figure 3: Scatterplot & Violin Plot (`Script_fig3.R`)
This script produces a two-panel composite plot:
* Panel A: Scatterplot of divergence (dN/dS) vs. polymorphism (pN/pS) with marginal density distributions. Highlights genes under positive selection based on BUSTED, BS-REL, or empirical dN/dS thresholds.
* Panel B: Violin plot displaying the Adaptation rate (Alpha) distribution, with staggered direct labeling for genes showing positive adaptation.
