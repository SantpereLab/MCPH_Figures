# TOP 5% CUT-OFF CALCULATION
# Data source: Supplementary Table (Sheet E) from:
# Enard, D., Cai, L., Gwennap, C., & Petrov, D. A. (2016). Viruses are a 
# dominant driver of protein adaptation in mammals. eLife, 5, e12469.
# https://doi.org/10.7554/eLife.12469

# ------------------------------------------------------------------------------
# STEP 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------

library(dplyr)  # Data manipulation package enabling mutate() and the pipe operator (%>%)
library(readxl) # Excel spreadsheet import

# ------------------------------------------------------------------------------
# STEP 2. LOAD SUPPLEMENTARY DATA (SHEET E)
# ------------------------------------------------------------------------------

# Sheet E contains human-chimpanzee divergence counts, specifically providing the 
# raw number of synonymous (DS) and non-synonymous (DN) substitutions per gene.

sheet_E <- read_excel("../../../../../../../Data/Supplementary_data_enard.xlsx", sheet = "E")

# ------------------------------------------------------------------------------
# STEP 3. APPLY HALDANE CORRECTION AND COMPUTE DN/DS
# ------------------------------------------------------------------------------

# Add +0.5 when DS == 0 to avoid division by zero without discarding genes,
# retaining 100% of genome-wide genes in the percentile distribution

sheet_E_haldane <- sheet_E %>%
  mutate(
    # Add 0.5 to DS if equal to 0; keep original value otherwise
    DS_adj = ifelse(DS == 0, DS + 0.5, DS),
    DN_DS = DN / DS_adj
  )

# ------------------------------------------------------------------------------
# STEP 4. COMPUTE 95TH PERCENTILE THRESHOLD (TOP 5% CUT-OFF)
# ------------------------------------------------------------------------------

# $: extracts the DN_DS column as a numeric vector
# quantile(..., 0.95): computes the 95th percentile (value above which the top 5% of genome-wide genes fall)
# na.rm = TRUE: ignores missing (NA) values during calculation to prevent errors

threshold_p95 <- quantile(sheet_E_haldane$DN_DS, 0.95, na.rm = TRUE)

print(threshold_p95) 

# Note: The +0.5 Haldane correction on DS=0 retains previously discarded genes, yielding high 
# DN/DS ratios (DN/0.5) that shift the 95th percentile threshold up from 3.0 to 4.0.
