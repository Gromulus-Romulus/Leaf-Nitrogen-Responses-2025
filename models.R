##' Compare multiple models of nitrogen saturation
##' for LMA, LDMC, CHL, and PhiPSII
##' 
##' @author [Nathan Malamud]
##' @date [2025-01-22]
##' 

# Libraries ----
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(scales)
library(smatr)
library(factoextra)

# Import Data ----
# REMINDER: Set Working Directory -> Source File Location
# Define factor levels as species
traits <-  read_csv("./data/traits.csv"),
traits$species <- factor(traits$species,
                         levels=c("R. sativus", "B. officinalis", "H. vulgare"))

# Styles ----
# Calculate rate of growth
growth_period_days <- 6 * 7 # 6 week experiment
traits$GRT <- (traits$dry_whole_g / growth_period_days)

# Filter by metrics of interest only
traits <- traits %>%
  select(barcodeID, species, treatment_mmol,
         LDMC, LMA, CHL, Phi_PS2, GRT)

# Define a custom theme for all plots
custom_theme <- theme_classic() +  # Base theme
  theme(
    # Text and font styling
    text = element_text(family = "sans", size = 12),
    axis.text = element_text(size = 10),  
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    
    # Axis labels
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # Margin for y-axis title
    
    # Panel and grid styling
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.1),  
    panel.grid.minor = element_blank(),  # No minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # White background
    
    # Aspect ratio
    aspect.ratio = 1  # 1:1 ratio
  )

# Custom colors advised by J. Garen
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Define units for variables
# TODO: format with Latex expressions
label_units <- c(
  "CHL" = "CHL (ug / cm²)",
  "LDMC" = "LDMC (mg / g)",
  "LMA" = "LMA (g / m²)",
  "Phi_PS2" = "PSII Fraction",
  "treatment_mmol" = "N (mM)"
)

# Model validation specifications ----
#   For each trait, regress onto nitrogen treatment.
#   Produce species x method table containing AIC values.
#   Produce PDF of model fits.

# TODO: models_to_try <- 

# Select model types to assess.
# TODO

## LMA model comparison ----
# TODO

## LDMC model comparison ----
# TODO

## CHL model comparison ----
# TODO

## PSII model comparison ----
# TODO

# Format results ----
# TODO