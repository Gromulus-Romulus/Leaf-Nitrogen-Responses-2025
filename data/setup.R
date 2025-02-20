## ------------------------------------------------------ ##
## Objective: Prepare raw traits data for analysis.
##
## Sections:
##   * Read in leaf mass Excel sheet
##   * Run ImageJ software on leaf scans
##   * Calculate LMA and LDMC from weight and area measurements
##   * Load fluorometry data from LI-600 and rename columns
##   * Write traits and fluorometry to CSV file
##
## Output: R script output will be fed into plsr.R for model 
##         fitting and cross-validation.
##
## Notes:
## - Dependencies: Requires ImageJ to be downloaded.
## - May need devtools for R package installations.
##
## @author Nathan Malamud
## @date   2024.09.20
## ------------------------------------------------------ ##

# ------------------------------------------------------ #
# LIBRARIES
library(tidyverse)
library(readxl)
library(reshape2)
library(LeafArea)

# ------------------------------------------------------ #
# Read in leaf mass data
# REMINDER: Set Working Directory -> Source File Location
mass_data <- read_excel("./leaf_mass_data.xlsx", sheet="trial2024")

# Define Factor Levels (treatment, species, barcode)
mass_data$treatment_mmol <- as.factor(mass_data$treatment_mmol)
mass_data$barcodeID <- as.factor(mass_data$barcodeID)
mass_data$LDMC <- as.numeric(mass_data$LDMC)
mass_data$species <- factor(mass_data$species, levels = c("radish", "borage", "barley"))

# ------------------------------------------------------ #
# RUN IMAGEJ SOFTWARE - ESTIMATE LEAF AREA
# This initiates a pop-up window and may take a while to run.
# Afterwards, merge values for image scans with the dataframe.
scans <- run.ij(set.directory = './scans/raw', path.imagej = '/Applications/ImageJ.app')
names(scans) <- c("barcodeID", "area_cm2")
#mass_data <- merge(mass_data, scans, by = "barcodeID")

# Allow missing values so that we can do data imputation
mass_data <- merge(mass_data, scans, by = "barcodeID", all.x = TRUE)

# ------------------------------------------------------ #
# CALCULATION OF LMA + EWT // 
# Follow formula reccomendations from Perez-Harguindeguy et al. (2016)
# Calculate LMA from weight and area measurements (units = g/m2)
mass_data$LMA <- ((mass_data$dry_leaf_g) /
                    (mass_data$area_cm2)) * (100^2)

# Calculate EWT from weight and area measurements (units = g/m2)
mass_data$EWT <- ((mass_data$wet_leaf_g - mass_data$dry_leaf_g) /
                    (mass_data$area_cm2)) * (100^2)

# ------------------------------------------------------ #
# LOAD FLUOROMETRY DATA
# Li-600 uses Fs and Fm' to calculate quantum yield of fluorescence
# Merge porometer measurements with data afterwards
#   Source: https://www.licor.com/env/products/LI-600/
fluor <- read_csv("./fluor/fluor_data.csv") %>%
  select(c(unique_id, gsw, gtw, gbw, PhiPS2))
names(fluor) <- c("barcodeID", "gsw", "gtw", "gbw", "Phi_PS2")

# ------------------------------------------------------ #
# LOAD PROSPECT-PRO OUTPUT
# Load PROSPECT-PRO output from run_prospect.R
prospect <- read_csv("./prospect.csv") %>%
  select(-c(sampleID, species, treatment_mmol))

# ------------------------------------------------------ #
# WRITE TRAITS DATA TO CSV
# Also merge measured traits with LI-COR fluorometry data
traits <- mass_data |>
  select("barcodeID",
         "sampleID",
         "treatment_mmol",
         "species",
         "dry_whole_g", "dry_leaf_g", "area_cm2",
         "LDMC", "LMA", "EWT")

# Remove duplicate values
traits <- traits %>%
  distinct(barcodeID, .keep_all = TRUE)

# Create a lookup table for species to Latin names
latin.names <- c("borage" = "B. officinalis",
                 "barley" = "H. vulgare",
                 "radish" = "R. sativus")

# Rename species column using the lookup table
traits <- traits %>%
  mutate(species = recode(species,
                          "borage" = latin.names["borage"],
                          "barley" = latin.names["barley"],
                          "radish" = latin.names["radish"])) %>%
  group_by(species)

# Merge with prospect-PRO and fluorometry values - allow NA for missing cells
traits <- merge(traits, prospect, by = "barcodeID", all.x = T)
traits <- merge(traits, fluor, by = "barcodeID", all.x = T)

# Merge duplicate rows
traits <- traits %>%
  distinct(barcodeID, .keep_all = TRUE)

# Save final traits data
write_csv(traits, file = ("./traits.csv"))
