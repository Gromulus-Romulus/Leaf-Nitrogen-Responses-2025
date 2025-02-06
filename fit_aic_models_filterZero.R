##' Compare pairwise combinations of variable transforms
##' for trait X and nitrogen treatment N.
##' 
##' @author [Nathan D. Malamud]
##' @date [2025-01-29]
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

# Filter out treatment_level == 0
# Avoid numerical issues
#traits <- traits %>%
#  filter(treatment_mmol > 0)
# traits <- traits %>%
#  mutate(treatment_mmol = treatment_mmol + 1e-6)

# Calculate rate of growth
growth_period_days <- 6 * 7 # 6 week experiment
traits$GRT <- (traits$dry_whole_g / growth_period_days)

# Filter by metrics of interest only
traits <- traits %>%
  select(barcodeID, species, treatment_mmol,
         LDMC, LMA, area_cm2, CHL, Phi_PS2, GRT)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# AIC model fitting - extract results ----

# Approach: define a script for each of these tasks.
source("fit_model_function.R")

traits_of_interest <- c("LMA", "LDMC", "CHL", "area_cm2", "GRT")
species_list <- unique(traits$species)

# Feed data to scripts and extract results
aic_model_stats <- data.frame(
  Trait = character(),
  Model = character(),
  AIC = numeric(),
  Singular = logical(),
  R2_marginal = numeric(),
  R2_conditional = numeric(),
  Slope = numeric(),
  P_Value = numeric()
)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
## No transformation, X ~ N ----
lme_models <- lapply(traits_of_interest,
                    function(trait) fit_mixed_models(traits, trait, "treatment_mmol",
                                                     species_list,
                                                     plot_visreg = TRUE,
                                                     plot_dir = "./visreg_plots/lme")) |> bind_rows()

lme_models$Formula <- "trait ~ nitrogen"
aic_model_stats <- rbind(aic_model_stats, lme_models)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
## Log transform one variable, X ~ log(N) ----
## Log transform one variable, log(X) ~ N ----
## Log transform both variables log(X) ~ log(N) ----
#
#   Note: Mengel and Kirkby 1987 suggest that the relationship
#   between growth and nutrient uptake is generally logarithmic.
#
# Apply log transformation to all traits of interest,
# adding a small value to avoid log(0) numerical errors.
#
log_traits <- traits %>%
  mutate(LMA = log(LMA),
         LDMC = log(LDMC),
         area_cm2 = log(area_cm2),
         CHL = log(CHL),
         GRT = log(GRT),
         log_treatment_mmol = log(treatment_mmol)
         )

traits <- traits %>%
  mutate(log_treatment_mmol = log(treatment_mmol))

# Fit all potential combinations of axis transforms.
# Account for convergence issues help('isSingular') by recording
# whether different models converged.
log_trait_models <- lapply(traits_of_interest,
                           function(trait)
                             fit_mixed_models(log_traits, trait,
                                              "treatment_mmol",
                                              species_list,
                                              plot_visreg = TRUE,
                                              plot_dir = "./visreg_plots/log_trait")) |> bind_rows()

log_trait_models$Formula <- "log(trait) ~ nitrogen"

log_nitrogen_models <- lapply(traits_of_interest,
                              function(trait)
                                fit_mixed_models(traits, trait, "log_treatment_mmol",
                                                 species_list,
                                                 plot_visreg = TRUE,
                                                 plot_dir = "./visreg_plots/log_nitrogen"
                                                 )) |> bind_rows()

log_nitrogen_models$Formula <- "trait ~ log(nitrogen)"

log_both_models <- lapply(traits_of_interest,
                          function(trait)
                            fit_mixed_models(log_traits, trait,
                                             "log_treatment_mmol", species_list,
                                             plot_visreg = TRUE,
                                             plot_dir = "./visreg_plots/log_both"
                                             )) |> bind_rows()

log_both_models$Formula <- "log(trait) ~ log(nitrogen)"

# Add to aic_model_stats dataframe
aic_model_stats <- rbind(aic_model_stats, log_trait_models, log_nitrogen_models, log_both_models)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
## Square root transformation X ~ sqrt(N) ----
sqrt_traits <- traits %>%
  mutate(LMA = sqrt(LMA),
         LDMC = sqrt(LDMC),
         area_cm2 = sqrt(area_cm2),
         CHL = sqrt(CHL),
         GRT = sqrt(GRT),
         sqrt_treatment_mmol = sqrt(treatment_mmol)
         )

traits <- traits %>%
  mutate(sqrt_treatment_mmol = sqrt(treatment_mmol))

# Do the same things as log models (last section), now with square roots.
sqrt_trait_models <- lapply(traits_of_interest,
                            function(trait)
                              fit_mixed_models(sqrt_traits, trait,
                                               "treatment_mmol",
                                               species_list,
                                               plot_visreg = TRUE,
                                               plot_dir = "./visreg_plots/sqrt_trait"
                                               )) |> bind_rows()
sqrt_trait_models$Formula <- "sqrt(trait) ~ nitrogen"

sqrt_nitrogen_models <- lapply(traits_of_interest,
                               function(trait)
                                 fit_mixed_models(traits, trait, "sqrt_treatment_mmol",
                                                  species_list,
                                                  plot_visreg=TRUE,
                                                  plot_dir = "./visreg_plots/sqrt_nitrogen"
                                                  )) |> bind_rows()
sqrt_nitrogen_models$Formula <- "trait ~ sqrt(nitrogen)"

sqrt_both_models <- lapply(traits_of_interest,
                           function(trait)
                             fit_mixed_models(sqrt_traits, trait,
                                              "sqrt_treatment_mmol",
                                              species_list,
                                              plot_visreg=TRUE,
                                              plot_dir = "./visreg_plots/sqrt_both"
                                              )) |> bind_rows()
sqrt_both_models$Formula <- "sqrt(trait) ~ sqrt(nitrogen)"

# Add to aic_model_stats dataframe
aic_model_stats <- rbind(aic_model_stats, sqrt_trait_models, sqrt_nitrogen_models, sqrt_both_models)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
## Reciprocal transformation 1/X ~ 1/N. Add small value to avoid 1/0 errors. ----
recip_traits <- traits %>%
  mutate(LMA = 1 / (LMA),
         LDMC = 1 / (LDMC),
         area_cm2 = 1 / (area_cm2),
         CHL = 1 / (CHL),
         GRT = 1 / (GRT),
         recip_treatment_mmol = 1 / (treatment_mmol)
         )

traits <- traits %>%
  mutate(recip_treatment_mmol = 1 / (treatment_mmol))

recip_nitrogen_models <- lapply(traits_of_interest,
                                function(trait)
                                  fit_mixed_models(recip_traits, trait, "recip_treatment_mmol",
                                                   species_list,
                                                   plot_visreg = TRUE,
                                                   plot_dir = "./visreg_plots/recip_nitrogen")) |> bind_rows()

recip_nitrogen_models$Formula <- "1/trait ~ 1/nitrogen"

# Add to aic_model_stats dataframe
aic_model_stats <- rbind(aic_model_stats, recip_nitrogen_models)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
## Reciprocal sqrt transformation 1/sqrt(X) ~ 1/sqrt(N). ----
# Add small value to avoid 1/0 errors.
recip_sqrt_traits <- traits %>%
  mutate(LMA = 1 / sqrt(LMA),
         LDMC = 1 / sqrt(LDMC),
         area_cm2 = 1 / sqrt(area_cm2),
         CHL = 1 / sqrt(CHL),
         GRT = 1 / sqrt(GRT),
         recip_sqrt_treatment_mmol = 1 / sqrt(treatment_mmol)
         )

traits <- traits %>%
  mutate(recip_sqrt_treatment_mmol = 1 / sqrt(treatment_mmol))

recip_sqrt_nitrogen_models <- lapply(traits_of_interest,
                                     function(trait)
                                       fit_mixed_models(recip_sqrt_traits, trait, "recip_sqrt_treatment_mmol",
                                                        species_list,
                                                        plot_visreg = TRUE,
                                                        plot_dir = "./visreg_plots/recip_sqrt_nitrogen"
                                                        )) |> bind_rows()

recip_sqrt_nitrogen_models$Formula <- "1/sqrt(trait) ~ 1/sqrt(nitrogen)"

# Add to aic_model_stats dataframe
aic_model_stats <- rbind(aic_model_stats, recip_sqrt_nitrogen_models)

# Add quadratic and cubic terms to the traits dataframe
traits <- traits %>%
  mutate(
    treatment_mmol_sq = treatment_mmol^2,       # Quadratic term
    treatment_mmol_cu = treatment_mmol^3        # Cubic term
  )

# Quadratic polynomial with respect to nitrogen X^2
quad_nitrogen_models <- lapply(traits_of_interest, function(trait) {
  fit_mixed_models(traits, trait, "treatment_mmol_sq",
                   species_list,
                   plot_visreg = TRUE,
                   plot_dir = "./visreg_plots/quad_nitrogen")
}) |> bind_rows()
quad_nitrogen_models$Formula <- "trait ~ nitrogen^2"

# Cubic polynomial with respect to nitrogen X^3
cubic_nitrogen_models <- lapply(traits_of_interest, function(trait) {
  fit_mixed_models(traits, trait, "treatment_mmol_cu",
                   species_list,
                   plot_visreg = TRUE,
                   plot_dir = "./visreg_plots/cubic_nitrogen")
}) |> bind_rows()
cubic_nitrogen_models$Formula <- "trait ~ nitrogen^3"

# Add to aic_model_stats dataframe
aic_model_stats <- rbind(aic_model_stats, quad_nitrogen_models, cubic_nitrogen_models)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Export AIC metrics and model fitting data to CSV ----
# - Reorganize columns for ease of readability.
aic_model_stats <- aic_model_stats %>% 
  select(Trait, Formula, Model, AIC, Singular, R2_marginal, R2_conditional, P_Value)
write_csv(aic_model_stats, "./aic_model_stats.csv")

# Based on R2 metrics, find which model formulas work best for each trait
aic_model_stats_filt <- aic_model_stats |>
  # Filter for non-singular models
  filter(!Singular) |> select(-Singular)

# Formula by trait table.
# The model that works best for ALL traits is trait ~ sqrt(nitrogen)
# Findings from this indicate that most working models have untransformed y values
aic_model_stats_filt |>
  select(Formula, Model) |>
  table()

aic_model_stats_filt.2 <- aic_model_stats_filt |>
  select(Trait, Model, Formula, AIC, R2_marginal, R2_conditional, P_Value) |>
  arrange(R2_conditional) |> filter(R2_conditional > 0.30) |>
  write_csv("./aic_model_stats_filt.csv")

# There is only 1 working model for LMA :
# - trait ~ log(Nitrogen)
table(aic_model_stats_filt.2$Trait)

# There are 5 working models for trait ~ log(Nitrogen)
# There are only 3 working models for trait ~ 1/Nitrogen
# These traits are... LMA, GRT, CHL, area_cm2, and LDMC
table(aic_model_stats_filt.2 %>% select(Model, Formula))

# This appears to be our best candidate - log transformation of nitrogen + linear regression
# AIC metrics indicate Random intercept works better than Random Intercept + Slope as well
aic_model_stats_filt.3 <- aic_model_stats_filt.2 %>%
  filter(Formula == "trait ~ log(nitrogen)", Model == "Random Intercept")

# Export the final model statistics for further analysis
write_csv(aic_model_stats_filt.3, "./aic_logNmodel_stats.csv")
