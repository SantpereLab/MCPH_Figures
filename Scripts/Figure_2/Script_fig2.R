# FIGURE 2: SLIM BAR PLOT HFVs by Cristina Torredemer 

# ------------------------------------------------------------------------------
# STEP 1. LOAD LIBRARIES AND CONFIGURE VECTOR FONTS
# ------------------------------------------------------------------------------

library(readxl)    # Excel file import
library(dplyr)     # Data manipulation
library(ggplot2)   # Core data visualization
library(ragg)      # High-definition raster rendering engine
library(showtext)  # Vector font embedding for PDF output
library(patchwork) # Multi-panel layout composition

# Enable showtext for crisp vector metrics and typographic fidelity

showtext_auto()
showtext_opts(dpi = 600)
update_geom_defaults("text", list(family = "sans"))

# ------------------------------------------------------------------------------
# STEP 2. IMPORT AND CLEAN DATA
# ------------------------------------------------------------------------------

df_raw <- read_excel("Data/HFVs_summary.xlsx")

df_clean <- df_raw %>%
  filter(!is.na(Gene)) %>% # Exclude summary/empty rows
  rename(
    MCPH_id  = MCPH,
    Gene     = Gene,
    n_HFV    = `# of HFV`,
    Ratio_aa = `Ratio HFV/aa (%)`
  ) %>%
  mutate(
    Gene     = trimws(gsub("\\s+", "", Gene)),       # Remove whitespace in gene symbols
    Ratio_aa = ifelse(n_HFV == 0, 0, Ratio_aa)       # Ensure zero-count consistency
  ) %>%
  arrange(desc(Ratio_aa), desc(n_HFV)) %>%           # Sort by density, then total count
  mutate(Gene = factor(Gene, levels = unique(Gene))) # Lock factor levels for plotting

# ------------------------------------------------------------------------------
# STEP 3. DEFINE CUSTOM THEME
# ------------------------------------------------------------------------------

theme_perfect <- function() {
  theme_bw(base_size = 7, base_family = "sans") +
    theme(
      # Panel tags ("a", "b") styling
      plot.tag = element_text(
        face = "bold", 
        size = 8.5, 
        color = "#374151", 
        margin = margin(l = 1, b = 2, unit = "pt")
      ),
      # Axis text formatting (vertical italicized gene labels on X-axis)
      axis.text.y = element_text(size = 5.5, color = "#374151"),
      axis.text.x = element_text(
        angle = 90, 
        vjust = 0.5, 
        hjust = 1, 
        face = "italic", 
        size = 5.5,
        color = "#374151",
        margin = margin(t = 2.5, unit = "pt")
      ),
      axis.title = element_text(size = 6.2, color = "#374151"),
      # Gridlines: keep subtle horizontal guides, remove vertical ones
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "#f3f4f6", linewidth = 0.2),
      panel.grid.minor.y = element_blank(),
      # Panel border, axis ticks, and outer margins
      panel.border       = element_rect(color = "#d1d5db", linewidth = 0.3, fill = NA),
      axis.ticks         = element_line(color = "#9ca3af", linewidth = 0.2),
      axis.ticks.length  = unit(1.2, "pt"),
      legend.position    = "none",
      plot.margin        = margin(t = 1, r = 2, b = 1, l = 2, unit = "pt")
    )
}

# ------------------------------------------------------------------------------
# STEP 4. BUILD PANELS
# ------------------------------------------------------------------------------

# Panel A: HFV density normalized by protein length (Purple palette)

p_a <- ggplot(df_clean, aes(x = Gene, y = Ratio_aa)) +
  geom_col(fill = "#a594c8", color = "#685489", linewidth = 0.2, width = 0.55) + # Slim bars
  geom_hline(yintercept = 0, color = "#d1d5db", linewidth = 0.3) +                 # Baseline at y = 0
  labs(tag = "a", x = NULL, y = "HFVs / aa (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +                      # Flush bottom, 5% top headroom
  theme_perfect() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())             # Hide X-axis (shared with Panel B)

# Panel B: Raw total HFV count per gene (Pink palette)

p_b <- ggplot(df_clean, aes(x = Gene, y = n_HFV)) +
  geom_col(fill = "#db8a9a", color = "#8e4757", linewidth = 0.2, width = 0.55) + # Slim bars
  geom_hline(yintercept = 0, color = "#d1d5db", linewidth = 0.3) +                 # Baseline at y = 0
  labs(tag = "b", x = "Gene", y = "Total HFVs") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +                      # Flush bottom, 5% top headroom
  theme_perfect()

# ------------------------------------------------------------------------------
# STEP 5. COMBINE PANELS AND EXPORT
# ------------------------------------------------------------------------------

# Stack panels vertically with equal height ratio (1:1)

figure_2 <- p_a / p_b + 
  plot_layout(heights = c(1, 1))

print(figure_2)

# Export vector PDF (Cairo) and high-resolution raster PNG (600 DPI via ragg)

ggsave("Figure_2.pdf", plot = figure_2, width = 180, height = 72, units = "mm", device = cairo_pdf)
ggsave("Figure_2.png", plot = figure_2, width = 180, height = 72, units = "mm", dpi = 600, device = agg_png)
