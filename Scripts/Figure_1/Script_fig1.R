# FIGURE 1: UNCLUSTERED BINARY DOT PLOT by Cristina Torredemer 

# ------------------------------------------------------------------------------
# STEP 1: DATA IMPORT AND PRE-PROCESSING
# ------------------------------------------------------------------------------

# Load required libraries 

library(ggplot2)    # Core package for data visualization and building the dot plot
library(dplyr)      # For data manipulation, filtering, and summarization
library(readxl)     # For importing data from the Gene_Function.xlsx file
library(forcats)    # For handling and reordering factor variables (e.g., sorting gene or function names)
library(ggh4x)      # For advanced facet customizations (useful for splitting the plot into the two functional blocks)
library(ragg)       # For rendering high-quality, anti-aliased graphics when exporting the plot
library(showtext)   # For using and customizing text fonts inside the plot

# Enable showtext for flawless vector rendering in PDF exports

showtext_auto()                                      # Automatically use showtext to render text in newly opened graphics devices
showtext_opts(dpi = 600)                             # Set the resolution to 600 DPI for high-quality text rendering
update_geom_defaults("text", list(family = "sans"))  # Set the default font family for all text geometries to 'sans' (sans-serif)

# Load and Prepare the Data 

data <- read_excel("./Data/Gene_Function.xlsx")      # Import the main dataset from the Excel file

# Rename the last column to remove the trailing space

data <- data %>%
  rename(Cell_Division_Related = `Cell Division-Related`)

# Create a new 'Block' column to categorize functions based on the 'Yes'/'No' values

data <- data %>%
  mutate(
    Block = ifelse(Cell_Division_Related == "Yes", 
                   "Cell Division-Related Processes", 
                   "Other Functions & Pathways")
  )

# Automatically calculate the counts for each block

counts <- data %>%
  group_by(Block) %>%
  summarise(
    n_pathways = n_distinct(Function),
    n_genes = n_distinct(Gene)
  )

# Dynamically create the subtitle text by extracting values from the 'counts' table

b1_p <- counts$n_pathways[counts$Block == "Cell Division-Related Processes"]
b1_g <- counts$n_genes[counts$Block == "Cell Division-Related Processes"]

b2_p <- counts$n_pathways[counts$Block == "Other Functions & Pathways"]
b2_g <- counts$n_genes[counts$Block == "Other Functions & Pathways"]

# Assigning the automated titles and subtitles

block_1_main <- "'Cell Division-Related Processes'"
block_1_sub  <- paste0("paste('(", b1_p, " pathways, ', italic(n), ' = ", b1_g, " genes)')")

block_2_main <- "'Other Functions & Pathways'"
block_2_sub  <- paste0("paste('(", b2_p, " pathways, ', italic(n), ' = ", b2_g, " genes)')")

# Assigning main and sub-block titles to each row based on their category

data <- data %>% 
  mutate(
    Block_Main = ifelse(Block == "Cell Division-Related Processes", block_1_main, block_2_main), # Set main title for the block
    Block_Sub  = ifelse(Block == "Cell Division-Related Processes", block_1_sub, block_2_sub)    # Set subtitle with parsed counts
  )

# Convert blocks to factors to ensure they are plotted in the correct order (Block 1 left, Block 2 right)

data$Block_Main <- factor(data$Block_Main, levels = c(block_1_main, block_2_main))
data$Block_Sub  <- factor(data$Block_Sub, levels = c(block_1_sub, block_2_sub))

# Calculate the order of functions (columns) sorted by the number of genes they involve within each block

func_order <- data %>% 
  group_by(Block_Main, Function) %>%                                  # Group by block and function name
  summarise(n = n_distinct(Gene), .groups = "drop") %>%               # Count unique genes per function
  arrange(Block_Main, desc(n)) %>%                                    # Sort by block first, then by highest gene count
  pull(Function)                                                      # Extract the sorted function names as a vector

# Apply the calculated order to the Function variable

data$Function <- factor(data$Function, levels = func_order)

# Calculate the order of genes (rows) sorted by their total functions and cell division-specific functions

gene_order <- data %>%
  group_by(Gene) %>%                                                  # Group by gene name
  summarise( 
    n_total = n_distinct(Function),                                   # Count total functions per gene
    n_cell_div = sum(Block == "Cell Division-Related Processes")      # Count how many functions belong to the first block
  ) %>% 
  arrange(n_total, n_cell_div, desc(Gene)) %>%                        # Sort by total functions, then block 1 functions, then alphabetically reverse
  pull(Gene)                                                          # Extract the sorted gene names as a vector

# Apply the calculated order to the Gene variable

data$Gene <- factor(data$Gene, levels = gene_order)

# ------------------------------------------------------------------------------
# STEP 2: GENERATE THE REFINED PLOT 
# ------------------------------------------------------------------------------

final_plot <- ggplot(data, aes(x = Function, y = Gene, fill = Block_Main, color = Block_Main)) + # Base aesthetics for x, y, and colors
  
  geom_point(shape = 21, size = 1.4, stroke = 0.35) +                 # Add dots with a filled center (shape 21) and colored border
  
  # Double faceting using ggh4x
  
  facet_nested(
    ~ Block_Main + Block_Sub,                                                        # Group columns first by Main Block, then by Sub-block
    scales = "free_x",                                                               # Let each block have its own independent x-axis labels
    space = "free_x",                                                                # Adjust panel width proportionally to the number of columns
    labeller = label_parsed,                                                         # Parse mathematical expressions (e.g., italic 'n') in titles
    strip = strip_nested(                                                            # Customize the visual appearance of the facet titles (strips)
      background_x = elem_list_rect(                                                 # Set background colors for the strip rectangles
        fill      = c("#f8f3ff", "#fbf0f5", "transparent", "transparent"),           # Soft purple and pink for top, transparent for bottom          
        color     = c("transparent", "transparent", "transparent", "transparent"),   # Remove borders from strips           
        linewidth = c(0, 0, 0, 0)                                                    # Set strip border width to zero
      ),
      text_x = elem_list_text(                                        # Customize text styling inside the strips
        color = c("#6b538e", "#a25a75", "#6b7280", "#6b7280"),        # Darker purple and pink for titles, gray for subtitles
        size  = c(6.5, 6.5, 5.6, 5.6),                                # Larger font for main titles, smaller for subtitles
        hjust = c(0.5, 0.5, 0.5, 0.5),                                # Center-align all text
        margin = list(                                                # Adjust vertical margins inside each strip box
          margin(t = 2.5, b = 2.0, unit = "pt"),                      # Top-left strip margin
          margin(t = 2.5, b = 2.0, unit = "pt"),                      # Top-right strip margin
          margin(t = 2.0, b = 2.5, unit = "pt"),                      # Bottom-left strip margin
          margin(t = 2.0, b = 2.5, unit = "pt")                       # Bottom-right strip margin
        )
      )
    )
  ) +
  
  # Harmonized color palettes (soft purple and cool pink/mauve)
  
  scale_fill_manual(values = c("#b1a1d4", "#e6b2c6")) +               # Set inner fill colors for the dots
  scale_color_manual(values = c("#756199", "#a25a75")) +              # Set outer border colors for the dots
  labs(tag = "a") +                                                   # Add panel tag 'a' in the top-left corner
  
  # Global theme and margin adjustments
  
  theme_bw(base_size = 7.5, base_family = "sans") +                   # Apply black-and-white base theme with sans-serif font
  theme(
    legend.position = "none",                                         # Hide the color legend since facets already explain the blocks
    strip.background = element_blank(),                               # Remove default background from strips to use custom nested ones
    
    # Position and styling of the 'a' tag
    
    plot.tag = element_text(
      face = "bold",                                                  # Make the tag bold
      size = 8.5,                                                     # Set tag font size
      color = "#374151",                                              # Set tag color to dark gray
      margin = margin(t = 8, l = 18, unit = "pt")                     # Adjust tag margins (top and left)
    ),
    
    axis.text.y = element_text(face = "italic", size = 5.6, color = "#4b5563"),                                     # Italicize gene names on y-axis and set color
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "plain", color = "#4b5563", size = 5.6),    # Rotate x-axis text 45 degrees
    panel.grid.major = element_line(color = "#f3f4f6", linewidth = 0.3),                                            # Set subtle color for major grid lines
    panel.grid.minor = element_blank(),                                                                             # Remove minor grid lines for a cleaner look
    panel.border = element_rect(color = "#e5e7eb", linewidth = 0.3, fill = NA),                                     # Soften the border around the plotting panels
    axis.title = element_blank(),                                                                                   # Remove 'Gene' and 'Function' axis titles to save space
    axis.ticks = element_line(color = "#e5e7eb", linewidth = 0.25),                                                 # Lighten the axis tick marks
    panel.spacing = unit(0.3, "lines"),                                                                             # Set spacing between the two facet panels
    
    # Extended right margin (r = 18pt) to prevent cutting off text at the edge
    
    plot.margin = margin(t = 2, r = 18, b = 2, l = 6, unit = "pt") # Overall plot margins (top, right, bottom, left)
  )

print(final_plot) # Render the plot to the graphics device

# ------------------------------------------------------------------------------
# STEP 3: HIGH-QUALITY EXPORT (600 DPI & VECTORIAL)
# ------------------------------------------------------------------------------

ggsave("Figure_1.pdf", plot = final_plot, width = 125, height = 108, units = "mm", device = cairo_pdf) # Export as a scalable vector PDF using cairo_pdf (ensures proper font embedding)
ggsave("Figure_1.png", plot = final_plot, width = 125, height = 108, units = "mm", dpi = 600, device = agg_png) # Export as a high-resolution 600 DPI PNG using agg_png for crisp anti-aliasing
