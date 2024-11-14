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

# TODO: move prospect_rtm_output to data
# Load prospect measurements from spec curves
prospect <- read_csv("./data/molecular_content.csv") %>%
  select(-c(sampleID, species, treatment_mmol))
data <- merge(data, prospect, by = "barcodeID")

# Define Factor Levels (Treatment and Species)
data$species <- factor(data$species, levels=c("R. sativus", "B. officinalis", "H. vulgare"))

# Define a custom theme for all plots
custom_theme <- theme_minimal(base_family = "sans") + 
  theme(
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25, color = "grey80"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    strip.text = element_text(hjust = 0, size = 10, face = "italic"),
    legend.position = "none"
  )

# // ------------------------------------------------------------- //
# Load necessary packages
library(tidyverse)
library(ggplot2)
library(GGally)
library(smatr)  # For regression analysis

# Define variables for regression analysis
variables <- c("CHL", "LMA", "LDMC", "dry_whole_g")

# Apply log10 transformations to relevant variables
data <- data %>%
  mutate(
    log_CHL = log10(CHL + 1),  # Adding 1 to avoid log10(0)
    log_LMA = log10(LMA + 1),
    log_LDMC = log10(LDMC + 1),
    log_dry_whole_g = log10(dry_whole_g + 1)
  )

# Define custom colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Function to perform pairwise regressions for each species
pairwise_species_regressions <- function(data, variables) {
  results <- list()
  
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      x <- variables[i]
      y <- variables[j]
      
      # Store the variable names and models for each species
      models <- data %>%
        group_by(species) %>%
        do(model = lm(reformulate(x, y), data = .))
      
      results[[length(results) + 1]] <- list(
        x = x,
        y = y,
        models = models
      )
    }
  }
  
  return(results)
}

# Perform species-specific regressions
regression_results <- pairwise_species_regressions(data, c("log_CHL", "log_LMA", "log_LDMC", "log_dry_whole_g"))

# Function to generate plots for each regression
generate_species_regression_plots <- function(results, data) {
  plots <- list()
  
  for (res in results) {
    x <- res$x
    y <- res$y
    models <- res$models
    
    # Create scatter plot with separate regression lines for each species
    plot <- ggplot(data, aes_string(x = x, y = y, color = "species")) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +  # Separate lines for each species
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 1:1 line
      scale_color_manual(values = josef_colors) +  # Use josef_colors for species
      labs(x = x, y = y) +
      custom_theme  # Apply custom theme
    
    plots[[length(plots) + 1]] <- plot
  }
  
  return(plots)
}

# Generate regression plots
regression_plots <- generate_species_regression_plots(regression_results, data)

# Display the plots
ggarrange(plotlist = regression_plots)
