##' Compare multiple models of nitrogen saturation
##' for plant functional traits using delta_AIC.
##' 
##' @author [Nathan Malamud]
##' @date [2025-01-27]
##' 

# Libraries ----
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(scales)
library(smatr)
library(RColorBrewer)
library(patchwork)
library(gridExtra)
library(cowplot)

# Styles ----
# Display first rows of the dataset
# Define consistent font size
base_font_size <- 10

# Define a custom theme for all plots
custom_theme <- theme_classic() +  # Base theme
  theme(
    # Text and font styling
    text = element_text(family = "sans", size = base_font_size),
    axis.text = element_text(size = 7),  
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "None",
    
    # No axis titles
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),

    # Title - left justified, 12 point font, bold
    plot.title = element_text(size = base_font_size, hjust = 0.5),
    
    # Panel and grid styling
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.1),  
    panel.grid.minor = element_blank(),  # No minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # White background
  )

# Define units for variables
label_units <- c(
  "LDMC" = "LDMC",
  "LMA" = "LMA",
  "area_cm2" = "Leaf Area",
  #"Phi_PS2" = "PSII",
  "CHL" = "CHL",
  "GRT" = "GRT",
  "treatment_mmol" = "N (mM)"
)

# Import Data ----
# REMINDER: Set Working Directory -> Source File Location
# Define factor levels as species
traits <-  read_csv("./data/traits.csv")

traits$species <- factor(traits$species,
                         levels=c("R. sativus",
                                  "B. officinalis",
                                  "H. vulgare"))

# Calculate rate of growth
growth_period_days <- 6 * 7 # 6 week experiment
traits$GRT <- (traits$dry_whole_g / growth_period_days)

# Filter by metrics of interest only
traits <- traits %>%
  select(barcodeID, species, treatment_mmol,
         LDMC, LMA, area_cm2, CHL, Phi_PS2, GRT)

# Model validation specifications ----
# Define helper function for log-model fitting ----
#   - Filters dataframe for one trait and one species
#   - Fits multiple transformation of variables
#   - Compares AIC values
fit_lm_models <- function(df, var, species_filter) {
  
  # Filter for the selected species
  df <- df %>%
    filter(species == species_filter) %>%
    select(species, treatment_mmol, all_of(var))
  
  # Rename variable for compatibility in models
  colnames(df)[colnames(df) == var] <- "var"
  
  # Define the list of model formulas
  model_formulas <- list(
    var ~ treatment_mmol,
    var ~ treatment_mmol**2,
    var ~ treatment_mmol**3,
    var ~ log(treatment_mmol + 0.00001),
    var ~ exp(treatment_mmol),
    var ~ sqrt(treatment_mmol)
  )
  
  # Fit all models and store them in a list
  models <- lapply(model_formulas, function(formula) lm(formula, data = df, method="qr"))
  
  # Function to extract model statistics
  model_stats <- function(model) {
    data.frame(
      AIC_val = AIC(model),
      R2 = summary(model)$r.squared,
      P_Value = summary(model)$coefficients[2, 4]
    )
  }
  
  # Extract statistics for all models
  model_statistics <- lapply(models, model_stats)
  
  # Combine statistics into a single data frame
  stats_df <- do.call(rbind, model_statistics)
  
  # Add model names for clarity
  stats_df$model <- paste0("Model ", seq_along(models))
  
  # Combine results into a single dataframe with proper labels
  model_labels <- c(
    "X ~ N", "X ~ N**2", "X ~ N**3", "X ~ log(N)", "X ~ exp(N)", "X ~ sqrt(N)"
  )
  
  # Add model labels
  stats_df$Model <- factor(model_labels, levels = model_labels)
  
  # Reorder columns for better readability
  stats_df <- stats_df[, c("Model", "AIC_val", "R2", "P_Value")]
  
  stats_df$Model <- factor(stats_df$Model,
                           levels = model_labels)
  
  return(stats_df)
}


# Initialize empty nested structure
model_results <- list() 

# Loop fits regression models to traits across species
for (species in levels(traits$species)) {
  # Fit models to separate species for traits of interest
  
  # Fit linear and non-linear models for leaf LMA
  lma <- fit_lm_models(traits, "LMA", species_filter = species)
  
  # Fit linear and non-linear models for LDMC
  ldmc <- fit_lm_models(traits, "LDMC", species_filter = species)
  
  # # Fit linear and non-linear models for leaf area
  acm2 <- fit_lm_models(traits, "area_cm2", species_filter = species)
  # 
  # # Fit linear and non-linear models for chlorophyll content
  chl <- fit_lm_models(traits, "CHL", species_filter = species)
  # 
  # # Fit linear and non-linear models for growth rate
  grt <- fit_lm_models(traits, "GRT", species_filter = species)
  
  # Combine results as entry in nested structure
  model_results[[species]] <- list(
    LMA = lma,
    LDMC = ldmc,
    area_cm2 = acm2,
    CHL = chl,
    GRT = grt
  )
}

# Save model results as CSV ----

# Define traits of interest
traits_of_interest <- c("LMA", "LDMC", "area_cm2", "CHL", "GRT")

# Initialize an empty data frame with the desired column names
model_frame <- data.frame(
  trait = character(),
  species = character(),
  model = character(),
  AIC_val = numeric(),
  R2 = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

# Collect model results for each species and trait
for (species in levels(traits$species)) {
  for (trait in traits_of_interest) {
    # Extract model statistics for this species and trait
    model_stats <- model_results[[species]][[trait]]
    
    # Create a temporary data frame for this species-trait combination
    temp_frame <- data.frame(
      trait = trait,
      species = species,
      model = model_stats$Model,
      AIC_val = model_stats$AIC_val,
      R2 = model_stats$R2,
      p_val = model_stats$P_Value,
      stringsAsFactors = FALSE
    )
    
    # Combine with the main data frame
    model_frame <- rbind(model_frame, temp_frame)
  }
}

# Order models in dataframe by length of formula
model_frame$model <- factor(model_frame$model, levels = unique(model_frame$model))

# Save the model results to a CSV file
write.csv(model_frame, "aic_results.csv", row.names = FALSE)

# View the result
head(model_frame)
