# FIGURE 4: POLYMORPHIC SITES BAR PLOTS (PANEL A AND PANEL B) by Cristina Torredemer

# ------------------------------------------------------------------------------
# STEP 1. LOAD LIBRARIES AND CONFIGURE VECTOR FONTS
# ------------------------------------------------------------------------------

library(readxl)     # Excel data import
library(dplyr)      # Data manipulation
library(ggplot2)    # Core visualization
library(patchwork)  # Panel assembly
library(ragg)       # High-definition raster rendering engine
library(showtext)   # Vector font embedding for PDF exports

# Enable showtext for strict typographic fidelity and vector metrics

showtext_auto()
showtext_opts(dpi = 600)
update_geom_defaults("text", list(family = "sans"))

# ------------------------------------------------------------------------------
# STEP 2. LOAD AND CLEAN DATASET
# ------------------------------------------------------------------------------

df_raw <- read_excel("./Data/polymorphic_positions.xlsx", na = c("", "NA", "N/A", "-"))

df_clean <- df_raw %>%
  filter(!is.na(Gene)) %>%
  rename(
    MCPH_id  = MCPH,
    Gene     = Gene,
    n_poly   = `Polymorphic positions all species`,
    Ratio_aa = `Polymorphism all species/alignment size (%)`
  ) %>%
  mutate(
    Gene     = trimws(gsub("\\s+", "", Gene)),
    Ratio_aa = suppressWarnings(as.numeric(Ratio_aa)),
    n_poly   = suppressWarnings(as.numeric(n_poly))
  ) %>%
  # Sort descending by Ratio_aa and n_poly; NAs are automatically placed at the end
  arrange(desc(Ratio_aa), desc(n_poly)) %>%
  mutate(Gene = factor(Gene, levels = unique(Gene)))

# ------------------------------------------------------------------------------
# STEP 3. CUSTOM VISUAL THEME (OVERLAP-FREE)
# ------------------------------------------------------------------------------

theme_perfect <- function() {
  theme_bw(base_size = 7, base_family = "sans") +
    theme(
      # Panel tags 'a' and 'b' positioned cleanly in the top-left corner
      plot.tag = element_text(
        face = "bold",
        size = 8.5,
        color = "#374151",
        margin = margin(b = 6, l = -14, unit = "pt")
      ),
      
      axis.text.y    = element_text(size = 5.5, color = "#374151"),
      # Y-axis title separated from tick labels (6pt)
      axis.title.y   = element_text(size = 6.2, color = "#374151", margin = margin(r = 6, unit = "pt")),
      
      axis.text.x    = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        face = "italic",
        size = 5.5,
        color = "#374151",
        margin = margin(t = 2.5, unit = "pt")
      ),
      axis.title.x   = element_text(size = 6.2, color = "#374151", margin = margin(t = 4, unit = "pt")),
      
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "#f3f4f6", linewidth = 0.2),
      panel.grid.minor.y = element_blank(),
      panel.border       = element_rect(color = "#d1d5db", linewidth = 0.3, fill = NA),
      axis.ticks         = element_line(color = "#9ca3af", linewidth = 0.2),
      axis.ticks.length  = unit(1.2, "pt"),
      legend.position    = "none",
      
      # Expanded left margin (l = 22pt) to accommodate large numbers, Y-title, and panel tag
      plot.margin        = margin(t = 8, r = 6, b = 2, l = 22, unit = "pt")
    )
}

# ------------------------------------------------------------------------------
# STEP 4. BUILD PANELS
# ------------------------------------------------------------------------------

# Panel A: Polymorphic sites / aa (%) ~ Purple palette

p_a <- ggplot(df_clean, aes(x = Gene, y = Ratio_aa)) +
  geom_col(fill = "#a594c8", color = "#685489", linewidth = 0.2, width = 0.55) +
  geom_hline(yintercept = 0, color = "#d1d5db", linewidth = 0.3) +
  # Add subtle "NA" label at the baseline for genes with missing Ratio_aa
  geom_text(
    data     = filter(df_clean, is.na(Ratio_aa)),
    aes(y = 0, label = "NA"),
    angle    = 90,
    hjust    = -0.2,
    vjust    = 0.5,
    size     = 1.6,
    color    = "#374151",
    fontface = "italic"
  ) +
  labs(tag = "a", x = NULL, y = "Polymorphic sites / aa (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_perfect() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Panel B: Total count of polymorphic sites ~ Pink palette

p_b <- ggplot(df_clean, aes(x = Gene, y = n_poly)) +
  geom_col(fill = "#db8a9a", color = "#8e4757", linewidth = 0.2, width = 0.55) +
  geom_hline(yintercept = 0, color = "#d1d5db", linewidth = 0.3) +
  # Add subtle "NA" label at the baseline for genes with missing n_poly
  geom_text(
    data     = filter(df_clean, is.na(n_poly)),
    aes(y = 0, label = "NA"),
    angle    = 90,
    hjust    = -0.2,
    vjust    = 0.5,
    size     = 1.6,
    color    = "#374151",
    fontface = "italic"
  ) +
  labs(tag = "b", x = "Gene", y = "No. of polymorphic sites") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_perfect()

# ------------------------------------------------------------------------------
# STEP 5. PANEL COMBINATION AND EXPORT
# ------------------------------------------------------------------------------

figura_4_final <- p_a / p_b + 
  plot_layout(heights = c(1, 1.25))

print(figura_4_final)

# Export vector PDF using Cairo and high-res raster PNG using ragg (agg_png)

ggsave("Figure_4.pdf", plot = figura_4_final, width = 180, height = 76, units = "mm", device = cairo_pdf)
ggsave("Figure_4.png", plot = figura_4_final, width = 180, height = 76, units = "mm", dpi = 600, device = agg_png)
