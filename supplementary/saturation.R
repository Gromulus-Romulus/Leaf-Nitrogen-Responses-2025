##' Create small multiples of different structural
##' and biochemical traits responding to nitrogen.
##' 
##' @author: Nathan D. Malamud
##' @date: 2021-11-11

library(tidyverse)
library(ggplot2)
library(ggpubr)

# Load data for traits
# REMINDER: Set Working Directory -> Source File Location
data <- read_csv("./data/traits.csv")

# Define colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Define the custom theme with a fixed aspect ratio
custom_theme <- theme_classic() +  # Start with a minimal theme
  theme(
    text = element_text(family = "sans", size = 12),  # Set font family and size
    axis.title.x = element_text(size = 12),  # Customize x-axis title
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # Add margin to y-axis title
    axis.text = element_text(size = 10),  # Customize axis text
    legend.text = element_text(size = 11),  # Customize legend text
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth=0.1),  # Customize major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and style the plot title
    aspect.ratio = 1  # Fix the aspect ratio (1:1)
  )

# Loop through each variable and create plots
# Define variables of interest
#variables_of_interest <- c("EWT", "LDMC", "LMA", "CHL", "CAR", "ANT", "Phi_PS2")
variables_of_interest <- c("LDMC", "LMA", "CHL", "Phi_PS2")

# Define y-axis labels with units
label_units <- c(
  #EWT = "EWT (g / m²)",
  Phi_PS2 = "Phi_PS2",
  LDMC = "LDMC (mg / g)",
  LMA = "LMA (g / m²)",
  CHL = "CHL (ug / cm²)"
  #CAR = "CAR (ug / cm²)",
  #ANT = "ANT (ug / cm²)"
)

# Initialize an empty list to store plots
plot_list <- list()

# Loop through each variable and create plots
for (variable in variables_of_interest) {
  # Filter data for the current variable
  data_filtered <- data %>%
    select(barcodeID, species, treatment_mmol, !!sym(variable)) %>%
    rename(value = !!sym(variable))
  
  # Create the plot with points and LOESS curves
  plot <- ggplot(data_filtered, aes(x = treatment_mmol, y = value, color = species, shape = species)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE, linetype = "solid") +  # LOESS curve in black
    scale_color_manual(values = josef_colors, name=NULL) +
    scale_shape_manual(values = c(16, 17, 18), name=NULL) +  # Set shapes for species
    labs(x = "N (mM)", y = label_units[[variable]]) +
    custom_theme + guides(color = guide_legend(override.aes = list(size = 5, alpha=1.0)))
  
  # Add the plot to the list
  plot_list[[variable]] <- plot
}

# Arrange the plots using ggarrange and assign to combined_plot
combined_plot <- ggarrange(
  plotlist = plot_list,
  nrow = 2, ncol=2,
  #labels = c("a", "b", "c", "d", "e", "f"),  # Add labels to each plot
  labels = c("a", "b", "c", "d"),  # Add labels to each plot
  common.legend = TRUE,
  legend = "bottom"
)

# Display the combined plot
print(combined_plot)

ggsave(filename = "./figures/traits_plot.tiff", 
       plot = combined_plot, device = "tiff", width = 5, height = 5)


