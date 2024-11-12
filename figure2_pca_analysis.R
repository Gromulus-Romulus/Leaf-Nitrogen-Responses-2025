##' View response of plants to nitrate treatments
##' in the trait PCA space
##' 
##' @author: Nathan D. Malamud
##' @date: 2021-11-11
##' 

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(factoextra)

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

# Create a new column for aggregated treatment levels
data <- data %>%
  mutate(treatment_level = case_when(
    treatment_mmol <= 15 ~ "0 - 15 mmol",
    treatment_mmol <= 35 ~ "20 - 35 mmol"
  ))

# TODO: decide whether log axis transformation is appropriate
# And adjust labels as such
# Log-transformed data
data_log10 <- data %>%
  mutate(
    LMA = log10(LMA),
    LDMC = log10(LDMC),
    N = log10(N),
    ETR = log10(ETR),
    CHL = log10(CHL),
    dry_whole_g = log10(dry_whole_g)
  )

# - - - - - 
# TODO: Example definition of pca_data (species across treatment levels)
pca_data <- data %>%
  select(species, treatment_level, CHL, N, LMA, LDMC, dry_whole_g)

# Structural trait PCA
# Split data by species
species_list <- split(pca_data, pca_data$species)

# Initialize empty lists to store PCA results and biplots
pca_results <- list()
pca_biplots <- list()

# Loop over each species in the species list
for (species_name in names(species_list)) {
  # Subset the data for the current species
  data_subset <- species_list[[species_name]]
  
  # Fit the PCA model
  pca <- prcomp(data_subset %>% select(CHL, N, LMA, LDMC, dry_whole_g), center = TRUE, scale. = TRUE)
  
  # Store the PCA result
  pca_results[[species_name]] <- pca
  
  # Generate the biplot and store it
  pca_biplots[[species_name]] <- fviz_pca_biplot(pca, 
                             geom.ind = "point", # Use points for individuals
                             col.ind = data_subset$treatment_level, # Color by treatment level
                             addEllipses = TRUE, # Add concentration ellipses
                             legend.title = "Treatment Level",
                             repel = TRUE, # Avoid label overlap
                             title = paste("PCA Biplot -", species_name))
}

# Arrange the biplots using ggarrange with 3 columns and 1 row
ggarrange(plotlist = pca_biplots, nrow=1, common.legend = TRUE)

# - - - - - - 
# TODO: Split data by treatment level
# Fit one PCA model for the entire dataset
# Define your custom colors for species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Structural trait PCA
# Split data by treatment level
treatment_list <- split(pca_data, pca_data$treatment_level)

# Initialize empty lists to store PCA results and biplots
pca_results <- list()
pca_biplots <- list()

# Loop over each treatment level in the treatment list
for (treatment_level in names(treatment_list)) {
  # Subset the data for the current treatment level
  data_subset <- treatment_list[[treatment_level]]
  
  # Fit the PCA model
  pca <- prcomp(data_subset %>% select(CHL, N, LMA, LDMC, dry_whole_g), center = TRUE, scale. = TRUE)
  
  # Store the PCA result
  pca_results[[treatment_level]] <- pca
  
  # Generate the biplot and store it
  pca_biplots[[treatment_level]] <- fviz_pca_biplot(
    pca, 
    geom.ind = "point", # Use points for individuals
    col.ind = data_subset$species, # Color by species
    addEllipses = TRUE, # Add concentration ellipses
    palette = josef_colors, # Use custom colors for species
    legend.title = "Species",
    repel = TRUE, # Avoid label overlap
    title = paste("PCA Biplot - Treatment Level", treatment_level)
  )
}

# Arrange the biplots using ggarrange with 1 row and a common legend
ggarrange(plotlist = pca_biplots, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
