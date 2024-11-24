# Author: Nathan Malamud
# Date: 2024.11.14
#

# Analysis Packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(scales)
library(smatr)
library(GGally)

# Load data for traits
# REMINDER: Set Working Directory -> Source File Location
data <- read_csv("./data/traits.csv")

# Define Factor Levels (Treatment and Species)
data$species <- factor(data$species, levels=c("R. sativus", "B. officinalis", "H. vulgare"))

# Define a custom theme for all plots
custom_theme <- theme_classic() +  # Start with a minimal theme
  theme(
    text = element_text(family = "sans", size = 12),  # Set font family and size
    axis.title.x = element_text(size = 12),  # Customize x-axis title
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # Add margin to y-axis title
    axis.text = element_text(size = 10),  # Customize axis text
    legend.text = element_text(size = 11),  # Customize legend text
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth=0.2),  # Customize major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and style the plot title
    aspect.ratio = 1  # Fix the aspect ratio (1:1)
  )

# Calculate metrics - unit leaf rate, leaf area ratio, leaf weight ratio, plant RGR
data$GRT <- (data$dry_whole_g / 6 * 7)

# // ------------------------------------------------------------- //
# Load necessary packages
library(tidyverse)
library(ggplot2)
library(GGally)
library(smatr)  # For regression analysis

# Define variables for regression analysis
variables <- c("CHL", "LMA", "LDMC", "GRT")

data <- data %>%
  filter(treatment_mmol > 0 & treatment_mmol < 35)

# Apply log10 transformations to relevant variables
data <- data %>%
  mutate(
    log_CHL = log10(CHL+1),  # Adding 1 to avoid log10(0)
    log_LMA = log10(LMA+1),
    log_LDMC = log10(LDMC+1),
    log_GRT = log10(GRT+1)
  )

# Define custom colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Function to perform pairwise regressions for each species using GLMs with link functions
pairwise_species_glms <- function(data, variables) {
  results <- list()
  
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      x <- variables[i]
      y <- variables[j]
      
      # Store the variable names and models for each species
      models <- data %>%
        group_by(species) %>%
        do(model = glm(reformulate(x, y), family = gaussian(link = "log"), data = .))
      
      results[[length(results) + 1]] <- list(
        x = x,
        y = y,
        models = models
      )
    }
  }
  
  return(results)
}

# Perform species-specific GLMs
glm_results <- pairwise_species_glms(data, c("CHL", "LMA", "LDMC"))

# Function to generate plots for each GLM regression with log scales
generate_species_glm_plots <- function(results, data) {
  plots <- list()
  
  for (res in results) {
    x <- res$x
    y <- res$y
    models <- res$models
    
    # Define custom labels with units (no log transformation in the labels)
    x_label <- switch(
      x,
      "CHL" = "CHL (μg/cm²)",         # Chlorophyll content in micrograms per square centimeter
      "LDMC" = "LDMC (mg/g)",         # Leaf dry matter content in milligrams per gram
      "LMA" = "LMA (g/m²)" 
    )
    
    y_label <- switch(
      y,
      "CHL" = "CHL (μg/cm²)",         # Chlorophyll content in micrograms per square centimeter
      "LDMC" = "LDMC (mg/g)",         # Leaf dry matter content in milligrams per gram
      "LMA" = "LMA (g/m²)" 
    )
    
    # Create scatter plot with separate GLM regression lines for each species
    plot <- ggplot(data, aes_string(x = x, y = y, color = "species", shape = "species")) +
      geom_point(alpha=0.55) +
      geom_smooth(method = "glm", method.args = list(family = gaussian(link = "log")), se = FALSE) +  # Fit GLM regression lines
      scale_x_log10() +  # Apply log scale to x-axis
      scale_y_log10() +  # Apply log scale to y-axis
      scale_color_manual(values = josef_colors, name=NULL) +  # Use josef_colors for species
      scale_shape_manual(values = c(16, 17, 18), name=NULL) +  # Use different shapes for species
      labs(x = x_label, y = y_label) +
      custom_theme  # Apply custom theme
    
    plots[[length(plots) + 1]] <- plot
  }
  
  return(plots)
}

# Generate updated GLM regression plots
glm_regression_plots <- generate_species_glm_plots(glm_results, data)

# Display the updated plots
ggarrange(plotlist = glm_regression_plots, labels = c("a", "b", "c"),
          nrow=1, common.legend = T, legend="bottom")
