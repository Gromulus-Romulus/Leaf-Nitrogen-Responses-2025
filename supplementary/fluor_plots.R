##' Objective: Plot differences in FV/FM / gsw
##' across low and high treatments
##'
##' Adapted from code written for 7th
##' annual plant functional traits course.
##' 
##' TODO: Do a wavelength-by-wavelength T test
##' and add p values to each chart
##'   Reference: https://journals.ashs.org/hortsci/view/journals/hortsci/53/5/article-p669.xml
##' 
##' @author [Nicole Bison, Nathan Malamud]
##' @date [2025.01.06]
##' 
##' Source:
##'   https://github.com/MichaletzLab/pftc7_spectroscopy
##'   

library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

# - - - - - 
# Read in input files and remove "date" column
# as it is irrelevant for now.
traits <- read.csv("../data/traits.csv")

custom_theme <- theme_classic() +  # Start with a minimal theme
  theme(
    text = element_text(family = "sans", size = 12, face="bold"),  # Set font family and size
    axis.title.x = element_text(size = 12),  # Customize x-axis title
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # Add margin to y-axis title
    axis.text = element_text(size = 10),  # Customize axis text
    legend.text = element_text(size = 11),  # Customize legend text
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.1),  # Customize major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background
    plot.title = element_blank(), # No title
    aspect.ratio = 1  # Fix the aspect ratio (1:1)
  )

# Define colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# - - - - -
# Helper function for making spectral plots with four treatment levels
# Create boxplots of Phi_PSII and gsw across treatment faceted by species
phipsII_plot <- ggplot(traits, aes(x = as.factor(treatment_mmol), y = Phi_PS2)) +
  geom_boxplot(aes(alpha = 0.2)) +
  geom_jitter(shape=1) +
  facet_wrap(~species) +
  labs(x = "Treatment", y = "Phi_PS2") +
  custom_theme + theme(legend.position = "none")

print(phipsII_plot)
  
gsw_plot <- ggplot(traits, aes(x = as.factor(treatment_mmol), y = gsw)) +
  geom_boxplot(aes(alpha = 0.2)) +
  geom_jitter(shape=1) +
  facet_wrap(~species) +
  labs(x = "Treatment", y = "gsw") +
  custom_theme + theme(legend.position = "none")

print(gsw_plot)

phys_plot <- ggarrange(phipsII_plot, gsw_plot, nrow = 2, labels = "AUTO")
  
# Save spectral Plot
ggsave(filename = "../figures/phys_plot.pdf", 
       plot = phys_plot, device = "pdf", width = 10, height = 6)
