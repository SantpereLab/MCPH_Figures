# FIGURE 3: SCATTERPLOT (CUTOFF) AND VIOLIN PLOT (ALL GENES)

# ------------------------------------------------------------------------------
# STEP 1. LOAD LIBRARIES AND CONFIGURE VECTOR FONTS
# ------------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggExtra)   
library(patchwork) 
library(ragg)      
library(showtext)  

showtext_auto()
showtext_opts(dpi = 600)
update_geom_defaults("text", list(family = "sans"))

# ------------------------------------------------------------------------------
# STEP 2. LOAD AND CLEAN MAIN DATASET
# ------------------------------------------------------------------------------

df_raw <- read_excel("../../../../../../../Data/Figure_generation/final_table.xlsx", na = c("", "NA", "N/A", "#N/A", "NaN"))

df_clean <- df_raw %>%
  filter(!is.na(GeneName)) %>%
  rename(
    Gene         = GeneName,
    pN_pS        = `pN/pS`,
    DN_DS        = `DN/DS`,
    Alpha        = Alphas,
    perce_pN_pS  = perce_pNpS,
    perce_DN_DS  = perce_DNDS,
    perce_Alpha  = perce_Alphas
  ) %>%
  mutate(
    Gene        = trimws(gsub("[\\*†]+", "", gsub("\\s+", "", Gene))),
    pN_pS       = as.numeric(pN_pS),
    DN_DS       = as.numeric(DN_DS),
    Alpha       = as.numeric(Alpha),
    perce_pN_pS = as.numeric(perce_pN_pS),
    perce_DN_DS = as.numeric(perce_DN_DS),
    perce_Alpha = as.numeric(perce_Alpha)
  )

set.seed(42)

# ------------------------------------------------------------------------------
# STEP 3. PUBLICATION-READY VISUAL THEMES 
# ------------------------------------------------------------------------------

theme_scatter <- function() {
  theme_bw(base_size = 7, base_family = "sans") +
    theme(
      axis.text          = element_text(size = 5.5, color = "#374151", face = "plain"),
      axis.title.x       = element_text(size = 6.2, color = "#374151", face = "plain", margin = margin(t = 5, unit = "pt")),
      axis.title.y       = element_text(size = 6.2, color = "#374151", face = "plain", margin = margin(r = 5, unit = "pt")),
      panel.background   = element_rect(fill = "white", color = NA),
      panel.grid.major   = element_line(color = "#EFEAEF", linewidth = 0.2),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "#d1d5db", linewidth = 0.3, fill = NA),
      axis.ticks         = element_line(color = "#9ca3af", linewidth = 0.2),
      axis.ticks.length  = unit(1.2, "pt"),
      legend.position    = "none",
      plot.margin        = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
    )
}

theme_violin <- function() {
  theme_bw(base_size = 7, base_family = "sans") +
    theme(
      axis.text.y        = element_text(size = 5.8, color = "#374151", face = "plain"),
      axis.title.y       = element_text(size = 6.8, color = "#374151", face = "plain", margin = margin(r = 5, unit = "pt")),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.title.x       = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "#EFEAEF", linewidth = 0.2),
      panel.background   = element_rect(fill = "white", color = NA),
      panel.border       = element_rect(color = "#d1d5db", linewidth = 0.3, fill = NA),
      axis.ticks.y       = element_line(color = "#9ca3af", linewidth = 0.2),
      axis.ticks.length  = unit(1.2, "pt"),
      legend.position    = "none",
      plot.margin        = margin(t = 18, r = 6, b = 2, l = 6, unit = "pt")
    )
}

# ------------------------------------------------------------------------------
# PANEL A: SCATTERPLOT (dN/dS vs pN/pS) WITH CUTOFF
# ------------------------------------------------------------------------------

df_scatter <- df_clean %>%
  filter(!is.na(pN_pS) & !is.na(DN_DS)) %>%
  mutate(
    BUSTED_p_value = as.numeric(BUSTED_p_value),
    Human          = as.numeric(Human),
    
    DN_DS_plot = case_when(
      Gene == "NCAPD2" ~ DN_DS - 0.03,
      Gene == "RRP7A"  ~ DN_DS + 0.03,
      TRUE             ~ DN_DS
    ),
    
    asterisks = case_when(
      !is.na(BUSTED_p_value) & BUSTED_p_value < 0.001 ~ "***",
      !is.na(BUSTED_p_value) & BUSTED_p_value < 0.01  ~ "**",
      !is.na(BUSTED_p_value) & BUSTED_p_value < 0.05  ~ "*",
      TRUE                                           ~ ""
    ),
    Gene_label = paste0(Gene, asterisks),
    
    is_bsrel_human = !is.na(Human) & Human > 0 & Human < 0.05,
    
    # Cutoff applied: Only highlight significant or dN/dS >= 4
    
    is_significant = (!is.na(BUSTED_p_value) & BUSTED_p_value < 0.05) | 
      is_bsrel_human |
      (!is.na(DN_DS) & DN_DS >= 4),
    
    point_shape = factor(if_else(is_bsrel_human, "23", "21"), levels = c("21", "23")),
    
    nudge_x = case_when(
      Gene == "CENPJ"    ~  0.22,
      Gene == "CEP152"   ~  0.00,
      Gene == "COPB2"    ~  0.22,
      Gene == "STIL"     ~  0.00,
      Gene == "CDK5RAP2" ~ -0.22,
      Gene == "NCAPD2"   ~ -0.28,
      Gene == "RRP7A"    ~  0.25,
      Gene == "WDFY3"    ~  0.34,
      Gene == "MFSD2A"   ~  0.00,
      Gene == "ZNF335"   ~  0.00,
      Gene == "AKNA"     ~  0.30,
      Gene == "ANKLE2"   ~  0.00,
      Gene == "WDR62"    ~  0.00,
      TRUE               ~  0.10
    ),
    nudge_y = case_when(
      Gene == "CENPJ"    ~  0.12,
      Gene == "CEP152"   ~  0.15,
      Gene == "COPB2"    ~  0.12,
      Gene == "STIL"     ~ -0.15,
      Gene == "CDK5RAP2" ~  0.08,
      Gene == "NCAPD2"   ~  0.12,
      Gene == "RRP7A"    ~  0.00,
      Gene == "WDFY3"    ~  0.00,
      Gene == "MFSD2A"   ~ -0.14,
      Gene == "ZNF335"   ~ -0.12,
      Gene == "AKNA"     ~  0.00,
      Gene == "ANKLE2"   ~ -0.15,
      Gene == "WDR62"    ~ -0.15,
      TRUE               ~  0.10
    )
  )

df_scatter_bg <- filter(df_scatter, !is_significant)
df_scatter_hl <- filter(df_scatter, is_significant)

p_scatter <- ggplot(df_scatter, aes(x = DN_DS_plot, y = pN_pS)) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "#d1d5db", linewidth = 0.3) +
  
  # Background points (muted)
  
  geom_point(
    data   = df_scatter_bg,
    aes(shape = point_shape),
    fill   = "#E6DCF0",
    color  = "#8A71B5",
    stroke = 0.3,
    size   = 1.1,
    alpha  = 0.85
  ) +
  
  # Highlighted points (vivid)
  
  geom_point(
    data   = df_scatter_hl,
    aes(shape = point_shape),
    fill   = "#937BBA",
    color  = "#5C4182",
    stroke = 0.3,
    size   = 1.4,
    alpha  = 0.95
  ) +
  
  scale_shape_manual(values = c("21" = 21, "23" = 23), limits = c("21", "23"), drop = FALSE) +
  
  # Labels only for highlighted points
  
  geom_text_repel(
    data          = df_scatter_hl,
    aes(label     = Gene_label),
    fontface      = "italic",
    nudge_x       = df_scatter_hl$nudge_x,
    nudge_y       = df_scatter_hl$nudge_y,
    size          = 1.6,
    color         = "#374151",
    bg.color      = "white",
    bg.r          = 0.10,
    segment.color = NA,
    box.padding   = 0.15,
    point.padding = 0.15,
    max.overlaps  = Inf
  ) +
  scale_x_continuous(limits = c(-0.2, 4.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.15, 3.55), expand = c(0, 0)) +
  labs(x = "dN/dS (Divergence)", y = "pN/pS (Polymorphism)") +
  theme_scatter()

p_scatter_marginal <- ggMarginal(
  p_scatter,
  type      = "density",
  margins   = "both",
  size      = 8,
  fill      = "#E6DCF0",
  color     = "#8A71B5",
  linewidth = 0.3,
  alpha     = 0.4
)

p_a_wrapped <- wrap_elements(p_scatter_marginal)

# ------------------------------------------------------------------------------
# PANEL B: ADAPTATION RATE (ALPHA) VIOLIN PLOT (ALL GENES LABELED)
# ------------------------------------------------------------------------------

df_c <- filter(df_clean, !is.na(Alpha)) %>%
  mutate(x_dummy = 0)

# No cutoff: All genes are arranged and labeled

df_c_hl <- df_c %>%
  arrange(desc(Alpha)) %>%
  mutate(
    Gene_label = gsub("[\\*†]+", "", Gene),
    rank_idx   = row_number(),
    side       = if_else(rank_idx %% 2 == 1, "right", "left")
  )

spread_y_positions <- function(y_vals, min_gap = 0.28) {
  if (length(y_vals) <= 1) return(y_vals)
  y_out <- y_vals
  for (i in 2:length(y_out)) {
    if ((y_out[i - 1] - y_out[i]) < min_gap) {
      y_out[i] <- y_out[i - 1] - min_gap
    }
  }
  return(y_out)
}

df_c_right <- df_c_hl %>%
  filter(side == "right") %>%
  mutate(y_text = spread_y_positions(Alpha, min_gap = 0.28))

df_c_left <- df_c_hl %>%
  filter(side == "left") %>%
  mutate(y_text = spread_y_positions(Alpha, min_gap = 0.28))

df_c_annot <- bind_rows(df_c_right, df_c_left) %>%
  mutate(
    x_pt        = 0.0,
    y_pt        = Alpha,
    x_seg_start = if_else(side == "right", 0.012, -0.012),
    y_seg_start = y_pt,
    x_seg_end   = if_else(side == "right", 0.07, -0.07),
    y_seg_end   = y_text,
    x_text      = if_else(side == "right", 0.08, -0.08),
    hjust       = if_else(side == "right", 0, 1)
  )

p_b <- ggplot(df_c, aes(x = x_dummy, y = Alpha)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#d1d5db", linewidth = 0.3) +
  geom_violin(fill = "#f4f0f9", color = "#d8cded", linewidth = 0.35, width = 0.45, alpha = 0.85) +
  
  # Highlighted top Alpha genes (ALL genes in this case)
  
  geom_point(data = df_c_annot, aes(x = x_pt, y = y_pt), shape = 21, fill = "#937BBA", color = "#5C4182", stroke = 0.3, size = 1.3, alpha = 0.95) +
  
  geom_segment(
    data = df_c_annot,
    aes(x = x_seg_start, y = y_seg_start, xend = x_seg_end, yend = y_seg_end),
    color = "#8A71B5", linewidth = 0.2, alpha = 0.60
  ) +
  
  geom_text(
    data = df_c_annot,
    aes(x = x_text, y = y_text, label = Gene_label, hjust = hjust),
    fontface = "italic", size = 2.0, color = "#374151"
  ) +
  labs(y = "Adaptation rate (α)") +
  theme_violin()

# ------------------------------------------------------------------------------
# STEP 4. PANEL COMBINATION AND EXPORT
# ------------------------------------------------------------------------------

p_combined <- (p_a_wrapped + p_b) + 
  plot_layout(widths = c(1.3, 0.7)) + 
  plot_annotation(tag_levels = 'a') & 
  theme(
    plot.tag          = element_text(face = "bold", size = 8.5, color = "#374151"),
    plot.tag.position = c(0.04, 0.99),
    plot.margin       = margin(t = 6, r = 6, b = 6, l = 4, unit = "pt")
  )

ggsave("Figure_3_All_violinplot_cutoff.pdf", plot = p_combined, width = 140, height = 80, units = "mm", device = cairo_pdf)
ggsave("Figure_3_All_violinplot_cutoff.png", plot = p_combined, width = 140, height = 80, units = "mm", dpi = 600, device = agg_png)