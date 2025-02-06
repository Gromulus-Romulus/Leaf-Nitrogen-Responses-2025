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
  "grt_g_d" = "grt_g_d",
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
traits$grt_g_d <- (traits$dry_whole_g / growth_period_days)

# Filter by metrics of interest only
traits <- traits %>%
  select(barcodeID, species, treatment_mmol,
         LDMC, LMA, area_cm2, CHL, Phi_PS2, grt_g_d)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# AIC model fitting - extract results ----

# Approach: define a script for each of these tasks.
source("fit_model_function.R")

traits_of_interest <- c("LMA", "LDMC", "CHL", "grt_g_d")
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
# Extract features from traits dataset ----
# by selecting for growth rate, LDMC, LMA, CHL and applying transforms
features <- traits %>% select(grt_g_d, species, LDMC, LMA, CHL)

# Add columns for log transformed LDMC, LMA, CHL
features <- features %>%
  mutate(log_LDMC = log(LDMC),
         log_LMA = log(LMA),
         log_CHL = log(CHL))

# Add columns for square root transformed LDMC, LMA, CHL
features <- features %>%
  mutate(sqrt_LDMC = sqrt(LDMC),
         sqrt_LMA = sqrt(LMA),
         sqrt_CHL = sqrt(CHL))

# Add columns for inverse transformed LDMC, LMA, CHL
features <- features %>%
  mutate(inv_LDMC = 1 / LDMC,
         inv_LMA = 1 / LMA,
         inv_CHL = 1 / CHL)

# Add columns for inverse square root transformed
features <- features %>%
  mutate(inv_sqrt_LDMC = 1 / sqrt(LDMC),
         inv_sqrt_LMA = 1 / sqrt(LMA),
         inv_sqrt_CHL = 1 / sqrt(CHL))

# Add columns for squared transformed LDMC, LMA, CHL
features <- features %>%
  mutate(sq_LDMC = LDMC^2,
         sq_LMA = LMA^2,
         sq_CHL = CHL^2)

# Add columns for cubed transformed LDMC, LMA, CHL
features <- features %>%
  mutate(cub_LDMC = LDMC^3,
         cub_LMA = LMA^3,
         cub_CHL = CHL^3)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Fit linear model for each trait + transform combo and growth rate
linear_mod_grt <- lm(grt_g_d ~ LDMC + LMA + CHL, data = features)
log_mod_grt <- lm(grt_g_d ~ log_LDMC + log_LMA + log_CHL, data = features)
sqrt_mod_grt <- lm(grt_g_d ~ sqrt_LDMC + sqrt_LMA + sqrt_CHL, data = features)
inv_mod_grt <- lm(grt_g_d ~ inv_LDMC + inv_LMA + inv_CHL, data = features)
inv_sqrt_mod_grt <- lm(grt_g_d ~ inv_sqrt_LDMC + inv_sqrt_LMA + inv_sqrt_CHL, data = features)
sq_mod_grt <- lm(grt_g_d ~ sq_LDMC + sq_LMA + sq_CHL, data = features)
cub_mod_grt <- lm(grt_g_d ~ cub_LDMC + cub_LMA + cub_CHL, data = features)

models <- list(linear_mod_grt, log_mod_grt, sqrt_mod_grt, inv_mod_grt, inv_sqrt_mod_grt, sq_mod_grt, cub_mod_grt)
names(models) <- c("linear", "log", "sqrt", "inv", "inv_sqrt", "sq", "cub")

# Define function to extract model statistics
model_stats_f <- function(mod) {
  summ <- summary(mod)
  
  data.frame(
    AIC = AIC(mod),
    R2 = summ$r.squared,  # Fixed effects only
    P_Value = anova(mod)$`Pr(>F)`[1]  # Extracting p-value for first predictor
  )
}

stats_df <- do.call(rbind, Map(function(mod, name) {
  df <- model_stats_f(mod)
  df$Model <- name
  df
}, models, names(models)))

# Print results
print(stats_df)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Load required packages
library(lme4)
library(MuMIn)

# Fit mixed-effects models with different trait transformations (species as random intercept)
linear_mod_grt <- lmer(grt_g_d ~ LDMC + LMA + CHL + (1 | species), data = features)  # No transformation
log_mod_grt <- lmer(grt_g_d ~ log_LDMC + log_LMA + log_CHL + (1 | species), data = features)  # Log-transformed traits
sqrt_mod_grt <- lmer(grt_g_d ~ sqrt_LDMC + sqrt_LMA + sqrt_CHL + (1 | species), data = features)  # Square-rooted traits
inv_mod_grt <- lmer(grt_g_d ~ inv_LDMC + inv_LMA + inv_CHL + (1 | species), data = features)  # Inverse-transformed traits
inv_sqrt_mod_grt <- lmer(grt_g_d ~ inv_sqrt_LDMC + inv_sqrt_LMA + inv_sqrt_CHL + (1 | species), data = features)  # Inverse sqrt
sq_mod_grt <- lmer(grt_g_d ~ sq_LDMC + sq_LMA + sq_CHL + (1 | species), data = features)  # Squared traits
cub_mod_grt <- lmer(grt_g_d ~ cub_LDMC + cub_LMA + cub_CHL + (1 | species), data = features)  # Cubed traits

# Store models in a named list
mixed_models <- list(
  "Linear" = linear_mod_grt,
  "Log" = log_mod_grt,
  "Sqrt" = sqrt_mod_grt,
  "Inverse" = inv_mod_grt,
  "Inverse Sqrt" = inv_sqrt_mod_grt,
  "Squared" = sq_mod_grt,
  "Cubed" = cub_mod_grt
)

# Function to extract model statistics
mixed_model_stats_f <- function(mod, name) {
  r2_vals <- r.squaredGLMM(mod)  # Get Marginal and Conditional R²
  
  data.frame(
    AIC = AIC(mod),
    Marginal_R2 = r2_vals[1, "R2m"],  # Fixed effects R²
    Conditional_R2 = r2_vals[1, "R2c"]  # Full model R² (fixed + random effects)
  )
}

# Sort results by AIC value
mixed_stats_df <- mixed_stats_df[order(mixed_stats_df$AIC), ]

# Print sorted results
print(mixed_stats_df %>% select(-Model))
