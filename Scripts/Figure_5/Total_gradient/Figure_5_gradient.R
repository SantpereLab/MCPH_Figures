# FIGURE 5: BRAIN TRAITS, CAAS DENSITY AND TOTAL COUNTS (PANELS A, B AND C) by Cristina Torredemer

# ------------------------------------------------------------------------------
# STEP 1. LOAD LIBRARIES AND CONFIGURE VECTOR FONTS
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ragg)
library(showtext)  # Vector font embedding for flawless PDF exports

# Enable showtext for strict typographic fidelity and vector metrics
showtext_auto()
showtext_opts(dpi = 600)
update_geom_defaults("text", list(family = "sans"))

# ------------------------------------------------------------------------------
# STEP 2. DEFINE MASTER GENE ORDER (33 GENES) AND CAAS DENSITIES
# ------------------------------------------------------------------------------
gene_order <- c(
  "CDK5RAP2", "AKNA", "KNL1", "MCPH1", "CENPJ", "CEP152", "KIF14", "ANKLE2",
  "STIL", "CENPE", "NCAPD3", "ASPM", "RRP7A", "WDR62", "CEP135",
  "NCAPH", "NCAPD2", "SASS6", "ZNF335", "MFSD2A", "LMNB2", "NUP37", "PHC1",
  "LMNB1", "WDFY3", "PDCD6IP", "CIT", "CDK6", "COPB2", "CDK4",
  "TRAPPC14", "SNRPE", "BUB1"
)

caas_kb_values <- c(
  MCPH1=130, WDR62=36, CDK5RAP2=3.5, KNL1=31, ASPM=86, CENPJ=173, STIL=0.5,
  CEP135=25, CEP152=161, ZNF335=0.2, PHC1=41, CDK6=5.1, CENPE=7.2, SASS6=56,
  MFSD2A=131, ANKLE2=12, CIT=13, WDFY3=3.4, COPB2=9.2, KIF14=137, NCAPD2=64,
  NCAPD3=55, NCAPH=93, NUP37=20, LMNB1=12, LMNB2=101, RRP7A=0, PDCD6IP=5,
  CDK4=2.1, AKNA=61, TRAPPC14=NA, SNRPE=NA, BUB1=NA
)

# ------------------------------------------------------------------------------
# STEP 3. LOAD EMPIRICAL P-VALUES FROM FINAL_TABLE.XLSX AND MAP ASTERISKS
# ------------------------------------------------------------------------------
table_path <- "../../Figure3/Data/Figure_generation/final_table.xlsx"
df_pvals <- read_excel(table_path) %>%
  filter(!is.na(GeneName)) %>%
  mutate(
    Gene_Clean = trimws(gsub("[*†]", "", GeneName)),
    p2 = suppressWarnings(as.numeric(CAAS_Pvalue_empirical_p2_Permutation)),
    Stars = case_when(
      is.na(p2) ~ "",
      p2 < 0.001 ~ "***",
      p2 < 0.01  ~ "**",
      p2 < 0.05  ~ "*",
      TRUE       ~ ""
    )
  )

# Create dynamic X-axis gene labels with asterisks BEFORE the gene name (so gene names sit closest to the heatmap)
gene_stars_map <- setNames(df_pvals$Stars, df_pvals$Gene_Clean)
gene_labels <- setNames(
  paste0(ifelse(is.na(gene_stars_map[gene_order]), "", gene_stars_map[gene_order]), gene_order),
  gene_order
)

# ------------------------------------------------------------------------------
# STEP 4. LOAD AND AGGREGATE MAIN EXCEL DATASET (AKNA & BUB1 CLEANING)
# ------------------------------------------------------------------------------
excel_path <- "../Data/Caas_MCPHgenes.xlsx"
sheets <- excel_sheets(excel_path)

brain_traits_full <- c(
  "Amygdala relative size", "Diencephalon relative size", "Hippocampus relative size",
  "Hypothalamus relative size", "Mesencephalon relative size", "Main olfactory bulb relative size",
  "Piriform lobe relative size", "Schizocortex relative size", "Septum relative size",
  "Striatum relative size", "Striatum including nucleus accens relative size",
  "Thalamus relative size", "Hypoglossal brainstem nuclei relative size"
)

brain_traits_short <- c(
  "Amygdala", "Diencephalon", "Hippocampus", "Hypothalamus", "Mesencephalon",
  "Main olfactory bulb", "Piriform lobe", "Schizocortex", "Septum",
  "Striatum", "Striatum inc. nucleus", "Thalamus", "Hypoglossal brainstem"
)

data_list <- list()
for (s in sheets) {
  if (s %in% gene_order) {
    df_sheet <- read_excel(excel_path, sheet = s)
    if ("Brain traits" %in% colnames(df_sheet) && "Counts" %in% colnames(df_sheet)) {
      sub_df <- df_sheet %>%
        filter(!is.na(`Brain traits`)) %>%
        select(`Brain traits`, Counts)
      
      # For the AKNA tab, retain only the actual bottom block (13 regions)
      if (s == "AKNA") {
        sub_df <- sub_df %>% tail(13)
      }
      
      sub_df$Gene <- s
      data_list[[s]] <- sub_df
    }
  }
}

df_long <- bind_rows(data_list) %>%
  filter(!is.na(`Brain traits`)) %>%
  mutate(Counts = suppressWarnings(as.numeric(Counts))) %>%
  filter(!is.na(Counts)) %>%
  group_by(`Brain traits`, Gene) %>%
  summarise(Counts = sum(Counts, na.rm = TRUE), .groups = "drop")

na_genes <- c("TRAPPC14", "SNRPE", "BUB1")

df_grid_brain <- expand.grid(
  `Brain traits` = brain_traits_full,
  Gene = gene_order,
  stringsAsFactors = FALSE
) %>%
  left_join(df_long, by = c("Brain traits", "Gene")) %>%
  mutate(
    Counts = case_when(
      Gene %in% na_genes ~ NA_real_,
      is.na(Counts) ~ 0,
      TRUE ~ Counts
    ),
    Gene = factor(Gene, levels = gene_order),
    Trait_Label = factor(`Brain traits`, levels = rev(brain_traits_full), labels = rev(brain_traits_short)),
    Font_Face = ifelse(!is.na(Counts) & Counts >= 100, "bold", "plain"),
    Is_White_Text = (!is.na(Counts) & Counts >= 100)
  )

# ------------------------------------------------------------------------------
# STEP 5. PREPARE SEPARATE TABLES FOR PANEL B (DENSITY) AND PANEL C (TOTALS)
# ------------------------------------------------------------------------------
# Panel B dataset: CAAS/size (kb)
df_density <- data.frame(
  Gene = names(caas_kb_values),
  Value = as.numeric(caas_kb_values)
) %>%
  filter(Gene %in% gene_order) %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Metric = factor("CAAS/size (kb)"),
    Label_Text = case_when(
      is.na(Value) ~ "NA",
      Value %% 1 == 0 ~ sprintf("%.0f", Value),
      TRUE ~ sprintf("%.1f", Value)
    ),
    Font_Face = ifelse(!is.na(Value) & Value >= 100, "bold", "plain"),
    Is_White_Text = (!is.na(Value) & Value >= 100),
    Text_Size = ifelse(is.na(Value), 1.2, 1.8)
  )

# Panel C dataset: Total counts (slate-blue gradient & white bold text >= 500)
df_totals <- df_grid_brain %>%
  group_by(Gene) %>%
  summarise(Value = if(all(is.na(Counts))) NA_real_ else sum(Counts, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Metric = factor("Total"),
    Label_Text = case_when(
      is.na(Value) ~ "NA",
      TRUE ~ sprintf("%.0f", Value)
    ),
    Font_Face = ifelse(!is.na(Value) & Value >= 500, "bold", "plain"),
    Is_White_Text = (!is.na(Value) & Value >= 500),
    Text_Size = ifelse(is.na(Value), 1.2, 1.8)
  )

# ------------------------------------------------------------------------------
# STEP 6. GENERATE PUBLICATION-READY HEATMAP PANELS (A, B, AND C)
# ------------------------------------------------------------------------------
theme_nature <- function() {
  theme_minimal(base_size = 7, base_family = "sans") +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = 3, r = 6, b = 1, l = 12, unit = "pt")
    )
}

# Panel A (Top): 13 Brain traits (lilac gradient)
p_main <- ggplot(df_grid_brain, aes(x = Gene, y = Trait_Label)) +
  geom_tile(fill = "#ffffff", color = "#e5e7eb", linewidth = 0.25) +
  geom_tile(data = filter(df_grid_brain, is.na(Counts)), fill = "#f3f4f6", color = "#e5e7eb", linewidth = 0.25) +
  geom_tile(data = filter(df_grid_brain, Counts > 0), aes(fill = Counts), color = "#e5e7eb", linewidth = 0.25) +
  geom_text(
    data = filter(df_grid_brain, Counts > 0),
    aes(label = Counts, color = Is_White_Text, fontface = Font_Face),
    size = 1.8
  ) +
  scale_color_manual(values = c("FALSE" = "#374151", "TRUE" = "#ffffff"), guide = "none") +
  scale_fill_gradientn(
    colors = c("#f4f0f8", "#d8ceeb", "#b5a3d7", "#8c72b8", "#614493", "#3d2268"),
    name = "CAAS count",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(8, "pt"),
      barheight = unit(45, "pt")
    )
  ) +
  labs(title = "a", y = "Brain traits (relative size)") +
  theme_nature() +
  theme(
    plot.title   = element_text(size = 9, face = "bold", color = "#374151", hjust = -0.06, vjust = 1),
    axis.title.y = element_text(size = 6, color = "#374151", face = "plain", margin = margin(r = 4, unit = "pt")),
    axis.text.y  = element_text(size = 5.5, color = "#374151", face = "plain"),
    axis.text.x  = element_blank(),
    legend.title = element_text(size = 5.5, color = "#374151", face = "bold", hjust = 0.5, margin = margin(b = 3, unit = "pt")),
    legend.text  = element_text(size = 5, color = "#374151")
  )

# Panel B (Middle): CAAS/size (kb) density (pink/burgundy gradient)
p_density <- ggplot(df_density, aes(x = Gene, y = Metric)) +
  geom_tile(fill = "#ffffff", color = "#e5e7eb", linewidth = 0.25) +
  geom_tile(data = filter(df_density, is.na(Value)), fill = "#f3f4f6", color = "#e5e7eb", linewidth = 0.25) +
  geom_tile(data = filter(df_density, !is.na(Value) & Value > 0), aes(fill = Value), color = "#e5e7eb", linewidth = 0.25) +
  geom_text(
    aes(label = Label_Text, color = Is_White_Text, fontface = Font_Face, size = I(Text_Size))
  ) +
  scale_color_manual(values = c("FALSE" = "#374151", "TRUE" = "#ffffff"), guide = "none") +
  scale_fill_gradientn(
    colors = c("#fcf0f3", "#f4c2cb", "#db8a9a", "#b85c70", "#7a293a"),
    name = "CAAS/size (kb)",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(8, "pt"),
      barheight = unit(16, "pt")
    )
  ) +
  labs(title = "b") +
  theme_nature() +
  theme(
    plot.title   = element_text(size = 9, face = "bold", color = "#374151", hjust = -0.06, vjust = 1),
    axis.text.y  = element_text(size = 5.5, color = "#374151", face = "plain"),
    axis.text.x  = element_blank(),
    legend.title = element_text(size = 5.5, color = "#374151", face = "bold", hjust = 0.5, margin = margin(b = 2, unit = "pt")),
    legend.text  = element_text(size = 5, color = "#374151")
  )

# Panel C (Bottom): Total CAAS counts (muted slate-blue gradient, white for 0, soft grey for NA)
p_total <- ggplot(df_totals, aes(x = Gene, y = Metric)) +
  geom_tile(fill = "#ffffff", color = "#e5e7eb", linewidth = 0.25) +
  geom_tile(data = filter(df_totals, !is.na(Value) & Value > 0), aes(fill = Value), color = "#e5e7eb", linewidth = 0.25) +
  geom_tile(data = filter(df_totals, is.na(Value)), fill = "#f3f4f6", color = "#e5e7eb", linewidth = 0.25) +
  geom_text(
    aes(label = Label_Text, color = Is_White_Text, fontface = Font_Face, size = I(Text_Size))
  ) +
  scale_color_manual(values = c("FALSE" = "#374151", "TRUE" = "#ffffff"), guide = "none") +
  scale_fill_gradientn(
    colors = c("#f1f5f9", "#cbd5e1", "#94a3b8", "#64748b", "#334155", "#1e293b"),
    name = "Total CAAS",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(8, "pt"),
      barheight = unit(16, "pt")
    )
  ) +
  scale_x_discrete(labels = gene_labels) +
  labs(x = "Gene", title = "c") +
  theme_nature() +
  theme(
    plot.title   = element_text(size = 9, face = "bold", color = "#374151", hjust = -0.06, vjust = 1),
    axis.text.y  = element_text(size = 5.5, color = "#374151", face = "plain"),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic", size = 5.5, color = "#374151"),
    axis.title.x = element_text(size = 6, color = "#374151", face = "plain", margin = margin(t = 4, unit = "pt")),
    legend.title = element_text(size = 5.5, color = "#374151", face = "bold", hjust = 0.5, margin = margin(b = 2, unit = "pt")),
    legend.text  = element_text(size = 5, color = "#374151")
  )

# ------------------------------------------------------------------------------
# STEP 7. PANEL COMBINATION AND HIGH-RESOLUTION EXPORT (600 DPI)
# ------------------------------------------------------------------------------
figura_5_final <- p_main / p_density / p_total + 
  plot_layout(heights = c(13, 1.3, 1.3))

print(figura_5_final)

ggsave("Figure_5_fixed.pdf", plot = figura_5_final, width = 180, height = 108, units = "mm", device = cairo_pdf)
ggsave("Figure_5_fixed.png", plot = figura_5_final, width = 180, height = 108, units = "mm", dpi = 600, device = agg_png)