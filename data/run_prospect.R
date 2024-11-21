## ------------------------------------------------------ ##
## Objective: Infer Chlorophyll (CHL), Carotenoid (CAR), 
##            and Anthocyanin (ANT) content using the 
##            empirically-validated method of PROSPECT-PRO 
##            inversion.
##
## Reference:
##   https://jbferet.gitlab.io/prospect/index.html
##
## @author Nathan Malamud
## @date   2024.09.27
##
## Chat-GPT4 Log:
##   https://chatgpt.com/share/66f7258c-28a4-8012-84b4-4f00c74ff8f6
## ------------------------------------------------------ ##

# ------------------------------------------------------ #
# LIBRARIES
library(tidyverse)
library(prospect)
library(dplyr)

# ------------------------------------------------------ #
# Load spectral data
spec_wide <- read_csv("./spec_data_wide.csv")

# Reorder columns
spec_wide <- spec_wide %>%
  select(barcodeID, sampleID, species, treatment_mmol, everything())

# Filter wavelengths of interest
FULL <- seq(500.0, 2400.0, .01)
NIR <- seq(700.0, 1100.0, .01)
SWIR <- seq(1100.0, 2400.0, .01)

# Remove duplicates
spec_wide <- spec_wide %>%
  distinct(barcodeID, .keep_all = TRUE) %>%
  distinct(sampleID, .keep_all = TRUE)

# Pivot dataframe for lambda, reflectance, and transmittance
spec_long <- spec_wide |>
  pivot_longer(
    cols = c(6:ncol(spec_wide)), 
    names_to = "lambda", 
    values_to = "reflectance"
  ) |>
  mutate(
    lambda = as.numeric(lambda),         # Convert lambda to numeric
    reflectance = as.numeric(reflectance) # Convert reflectance to numeric
  ) %>%
  filter(lambda %in% FULL)

# Aggregate lambda values by rounding and calculate mean reflectance
spec_long$lambda <- round(spec_long$lambda, 0)
spec_long_agg <- aggregate(
  reflectance ~ barcodeID + species + treatment_mmol + lambda,
  data = spec_long,
  FUN = mean,
  na.rm = TRUE
)

# ------------------------------------------------------ #
# Estimate parameters using PROSPECT-PRO
# Define set of parameters to estimate
Parms2Estimate <- c('N', 'CHL', 'CAR', 'ANT')

# Create an empty list to store results for each leaf
results_list <- list()

# Get unique barcodeIDs
unique_barcodes <- unique(spec_long_agg$barcodeID)

# Loop through each barcodeID
for (barcode in unique_barcodes) {
  
  # Subset data for current barcode
  leaf_spec <- spec_long_agg %>% filter(barcodeID == barcode)
  
  # Adjust spectral domain
  SubData <- FitSpectralData(
    lambda = leaf_spec$lambda,
    Refl = leaf_spec$reflectance, 
    Tran = NULL
  )
  
  # Invert PROSPECT to estimate parameters
  OutPROSPECTPRO <- Invert_PROSPECT(
    SpecPROSPECT = SubData$SpecPROSPECT, 
    Refl = SubData$Refl, 
    Tran = NULL,
    Parms2Estimate = Parms2Estimate, 
    PROSPECT_version = 'PRO'
  )
  
  # Extract estimated parameters (ug / cm2) and store with metadata
  estimated_params <- data.frame(
    barcodeID = barcode,
    N = OutPROSPECTPRO$N,
    CHL = OutPROSPECTPRO$CHL,
    CAR = OutPROSPECTPRO$CAR,
    ANT = OutPROSPECTPRO$ANT
  )
  
  # Append to results list
  results_list[[as.character(barcode)]] <- estimated_params
}

# Combine results into a single dataframe
final_results <- do.call(rbind, results_list)

# Merge with original data for metadata
final_results <- merge(
  spec_wide %>%
    select(barcodeID, sampleID, species, treatment_mmol),
  final_results,
  by = "barcodeID"
)

# Write results to CSV
write.csv(final_results, file = "./prospect.csv", row.names = FALSE)
