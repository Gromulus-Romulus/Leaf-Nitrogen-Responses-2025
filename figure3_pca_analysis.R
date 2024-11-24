# PCA Visualization of Plant Trait Responses to Nitrate Treatments
# 
# Author: Nathan D. Malamud
# Date: 2021.11.11

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(factoextra)

# Load data for traits
# REMINDER: Set Working Directory -> Source File Location
data <- read_csv("./data/traits.csv")

# Define a custom theme for all plots
custom_theme <- theme_classic() +  # Start with a minimal theme
  theme(
    text = element_text(family = "sans", size = 12),  # Set font family and size
    axis.title.x = element_text(size = 12),  # Customize x-axis title
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # Add margin to y-axis title
    axis.text = element_text(size = 10),  # Customize axis text
    legend.text = element_text(size = 11),  # Customize legend text
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.2),  # Customize major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background
    plot.title = element_blank(), # No title
    aspect.ratio = 1  # Fix the aspect ratio (1:1)
  )

# Define custom colors for species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Select relevant columns and preprocess the data
pca_data <- data %>%
  select(species, CHL, LMA, LDMC, EWT, dry_whole_g, treatment_mmol) %>%
  drop_na() %>%  # Remove any rows with NA values
  mutate(
    GRT = dry_whole_g / 6 * 7  # Calculate growth rate (example metric)
  ) %>%
  select(species, CHL, LMA, LDMC, GRT, EWT, treatment_mmol)

# Fit a PCA model
pca <- prcomp(pca_data %>% select(-species), center = TRUE, scale. = TRUE)

# Generate the PCA biplot with species ellipsoids
pca_plot <- fviz_pca_biplot(
  pca, 
  geom.ind = "point",  # Use points for individuals
  col.ind = pca_data$species,  # Color by species
  palette = josef_colors,  # Use custom colors for species
  addEllipses = TRUE,  # Add concentration ellipses for species
  legend.title = "Species", 
  repel = TRUE  # Avoid label overlap
) + custom_theme

# Display the plot
print(pca_plot)
