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

# - - - - - 
# Categorize treatments into "LO" and "HI" based on set threshold for mmol
# Create a new column for aggregated treatment levels
spectra_long <- spectra_long %>%
  mutate(treatment_level = case_when(
    treatment_mmol <= 5 ~ "0 - 5 mmol",
    treatment_mmol <= 15 ~ "10 - 15 mmol",
    treatment_mmol <= 25 ~ "20 - 25 mmol",
    treatment_mmol <= 35 ~ "30 - 35 mmol",
    TRUE ~ "Other"  # Optional: catch any values above 35
  ))

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


# - - - - -
# Average measurements per-species and treatment
per_species_treatment <- spectra_long %>% group_by(species, treatment_level, wavelength) %>% 
  summarize(mean_r = mean(reflectance),
            sd_r = sd(reflectance),
            n = n(),
            se_r = sd_r / sqrt(n)) %>% ungroup()

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
spectra_plot <- ggplot(per_species_treatment, aes(x = wavelength, color = treatment_level, fill = treatment_level)) +
    geom_ribbon(aes(ymin = 100 * (mean_r - 2 * se_r), ymax = 100 * (mean_r + 2 * se_r)), 
                alpha = 0.3, color = NA) +
    geom_line(aes(y = 100 * mean_r), linewidth = 0.75) +  # Thicker lines for better visibility
    ylab("% Reflectance") + 
    xlab("Wavelength (nm)") +
    scale_fill_manual(name = "Treatment Level", 
                      values = c("0 - 5 mmol" = "#feb24c", 
                                 "10 - 15 mmol" = "#bda48e", 
                                 "20 - 25 mmol" = "#c9a9a9", 
                                 "30 - 35 mmol" = "#f20000")) +
    scale_color_manual(name = "Treatment Level", 
                       values = c("0 - 5 mmol" = "#feb24c", 
                                  "10 - 15 mmol" = "#bda48e", 
                                  "20 - 25 mmol" = "#c9a9a9", 
                                  "30 - 35 mmol" = "#f20000")) +
    facet_wrap(~species) + custom_theme
  
# Save spectral Plot
ggsave(filename = "../figures/spectral_plot.pdf", 
       plot = spectra_plot, device = "pdf", width = 10, height = 6)

# Show in plot viewer
print(spectra_plot)

