##' Objective: Plot differences in spectral signatures
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
##' @date [2024.10.02]
##' 
##' Source:
##'   https://github.com/MichaletzLab/pftc7_spectroscopy
##'   

library(dplyr)
library(tidyverse)
library(RColorBrewer)

# - - - - - 
# Read in input files and remove "date" column
# as it is irrelevant for now.
spectra_long <- read.csv("./spec_data_long.csv") |>
  select(-date)

spectra_wide <- read.csv("./spec_data_wide.csv") |>
  select(-date)

# Remove "X" characters from beginning of wavelengths
spec_names <- colnames(spectra_wide)[5:ncol(spectra_wide)]
colnames(spectra_wide)[5:ncol(spectra_wide)] <- gsub("X", "", spec_names)

custom_theme <- theme_classic() +  # Start with a minimal theme
  theme(
    text = element_text(family = "sans", size = 12),  # Set font family and size
    axis.title.x = element_text(size = 12),  # Customize x-axis title
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # Add margin to y-axis title
    axis.text = element_text(size = 10),  # Customize axis text
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.text = element_text(size = 11),  # Customize legend text
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.1),  # Customize major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background
    plot.title = element_blank(), # No title
    aspect.ratio = 1  # Fix the aspect ratio (1:1)
  )


# - - - - -
# Average measurements per-species and treatment
per_species_treatment <- spectra_long %>% group_by(species, treatment_mmol, wavelength) %>% 
  summarize(mean_r = mean(reflectance),
            sd_r = sd(reflectance),
            n = n(),
            se_r = sd_r / sqrt(n)) %>% ungroup()

# Relable treatment levels with units
per_species_treatment <- per_species_treatment %>%
  mutate(treatment_mmol = case_when(
    treatment_mmol == 0 ~ "0 mM",
    treatment_mmol == 5 ~ "5 mM", 
    treatment_mmol == 10 ~ "10 mM",
    treatment_mmol == 15 ~ "15 mM",
    treatment_mmol == 20 ~ "20 mM",
    treatment_mmol == 25 ~ "25 mM",
    treatment_mmol == 30 ~ "30 mM",
    treatment_mmol == 35 ~ "35 mM"
  ))

# Set levels so it's in order: 0, 5, ... 25, ... 30, 35
per_species_treatment$treatment_mmol <- factor(per_species_treatment$treatment_mmol, levels = c("0 mM", "5 mM", "10 mM", "15 mM", "20 mM", "25 mM", "30 mM", "35 mM"))

# Assign wavelength ranges
# Don't include below < 400 and above 2400
# due to increased noise at range edges
# Ranges set using SVC i-series Field Spectroscopy guide
FULL = seq(400, 2400, by = .1)
VIS = seq(400, 700, by = .1)
NIR = seq(700, 1000, by = .1)
SWIR = seq(1000, 2400, by = .1)

spectra_long <- spectra_long %>%
  mutate(wave_type = case_when(
    wavelength %in% VIS ~ "Visible",
    wavelength %in% NIR ~ "NIR",
    wavelength %in% SWIR ~ "SWIR"
  ))

# - - - - -
# Helper function for making spectral plots with four treatment levels
spectra_plot <- ggplot(per_species_treatment, aes(x = wavelength, color = "black", fill = "grey")) +
    geom_ribbon(aes(ymin = 100 * (mean_r - 2 * se_r), ymax = 100 * (mean_r + 2 * se_r)), 
                alpha = 0.3, color = "grey") +
    geom_line(aes(y = 100 * mean_r), linewidth = 0.75, color = "red") +  # Thicker lines for better visibility
    ylab("Average % Reflectance") + 
    xlab("Wavelength (nm)") +
    facet_grid(species~treatment_mmol) + custom_theme +
  theme(legend.position = "None")
  
# Save spectral Plot
ggsave(filename = "../figures/spectral_plot.pdf", 
       plot = spectra_plot, device = "pdf", width = 11, height = 8.5)

# Show in plot viewer
print(spectra_plot)

