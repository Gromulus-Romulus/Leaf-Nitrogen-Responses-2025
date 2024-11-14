##' Create small multiples of different structural
##' and biochemical traits responding to nitrogen.
##' 
##' @author: Nathan D. Malamud
##' @date: 2021-11-11
##' 

library(tidyverse)
library(ggplot2)
library(ggpubr)

# Load data for traits
# REMINDER: Set Working Directory -> Source File Location
data <- read_csv("./data/traits.csv")

# TODO: move prospect_rtm_output to data
# Load prospect measurements from spec curves
prospect <- read_csv("./data/molecular_content.csv") %>%
  select(-c(sampleID, species, treatment_mmol))
data <- merge(data, prospect, by = "barcodeID")

# Define Factor Levels (Treatment and Species)
data$treatment_mmol <- data$treatment_mmol
data$species <- factor(data$species, levels=c("R. sativus", "B. officinalis", "H. vulgare"))
data$LDMC <- as.numeric(data$LDMC)

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

# Assign factor levels
data$species <- factor(data$species, levels = c("R. sativus", "B. officinalis", "H. vulgare"))

# Ensure numeric values (TODO: quality check all of this)
data$treatment_mmol <- as.numeric(data$treatment_mmol) 
data$Qamb <- as.double(data$Qamb)

# Define colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Define average growth rate for each species, assume 6 week interval
data$shoot_growth <- data$dry_whole_g / (6 * 7) # gram per day
data$CHL_cm2 <- data$CHL
data$CHL_g <- data$CHL_cm2 / (data$LMA) # TODO: check your units here
data$ETR <- (data$ETR / data$Qamb) * 50.0 # Based on PhiPSII measurements

# Filter out traits of interest
data <- data %>% select(barcodeID, species, treatment_mmol,
  LDMC, # Leaf dry matter content
  LMA, # Leaf mass per area
  CHL, # Chlorophyll pigment content
  N, # Effective number of mesophyll layers
  #ETR, # Expected rate of electron transport
  dry_whole_g # Dry weight of whole plant
)

# TODO: decide whether log axis transformation is appropriate
# And adjust labels as such
# Log-transformed data
data_log10 <- data %>%
  mutate(
    LMA = log10(LMA),
    LDMC = log10(LDMC),
    N = log10(N),
    #ETR = log10(ETR),
    CHL = log10(CHL),
    dry_whole_g = log10(dry_whole_g)
  )

# Create long form of data
data_long <- data %>% pivot_longer(
  cols = -c(barcodeID, species, treatment_mmol), names_to = "trait",
  values_to = "value")

# Updated Y-axis labels with units
label_units <- c(
  CHL = "CHL (μg/cm²)",
  dry_whole_g = "AB Dry Weight (g)",
  LDMC = "LDMC (mg/g)",
  LMA = "LMA (g/m²)",
  N = "N Layers"
)

# Create small multiples of boxplots with updated labels
ggplot(data_long, aes(x = as.factor(treatment_mmol), y = value, fill = species)) +
  geom_boxplot() + 
  facet_grid(trait ~ species, scales = "free", labeller = labeller(trait = label_units)) +
  scale_fill_manual(values = josef_colors) +
  custom_theme +
  labs(x = "Ammoniacal Nitrogen (mM)", y = NULL)

