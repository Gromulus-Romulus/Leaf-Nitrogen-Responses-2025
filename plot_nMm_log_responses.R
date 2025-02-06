##' Plot reciprocal transformeed
##' responses with random intercepts
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

# Model fitting ----
library(lme4)
library(MuMIn)

# Import data
traits <- read_csv("./data/traits.csv")

# Ensure species is a factor
traits$species <- factor(traits$species, 
                         levels = c("R. sativus",
                                    "B. officinalis",
                                    "H. vulgare"))

# Compute growth rate per day
growth_period_days <- 6 * 7  # 6-week experiment
traits$GRT <- traits$dry_whole_g /
  growth_period_days

# Filter only by traits of interest
traits <- traits %>% select(
  species, LDMC, LMA, CHL, GRT, treatment_mmol
  )

# Styles ----
# Display first rows of the dataset
# Define consistent font size
base_font_size <- 12

# Define a custom theme for all plots
custom_theme <- theme_classic() +  # Base theme
  theme(
    # Text and font styling
    text = element_text(family = "sans", size = base_font_size),
    axis.text = element_text(size = 9),  
    legend.position = "bottom",
    
    # Title - left justified, 12 point font, bold
    plot.title = element_text(size = base_font_size, hjust = 0.5),
    
    # Aspect ratio 1:1
    aspect.ratio = 1,
    
    # Panel and grid styling
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.1),  
    panel.grid.minor = element_blank(),  # No minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # White background
  )

# Custom colors advised by J. Garen
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Define units for variables with LaTeX formatting
label_units <- list(
  "LDMC" = expression("LDMC"~"("*mg~g^-1*")"),
  "LMA" = expression("LMA"~"("*g~m^-2*")"),
  "CHL" = expression("CHL"~"("*mu*g~cm^-2*")"),
  "treatment_mmol" = "Nitrogen (mM)"  # Treatment nitrogen label remains unchanged
)

# Import Data ----
# REMINDER: Set Working Directory -> Source File Location
# Define traits of interest
traits_of_interest <- c("LMA", "LDMC", "CHL")

# Apply specific transformations for each trait
traits <- traits %>%
   mutate(
     LMA_trans = I(LMA),
     LDMC_trans = I(LDMC),
     CHL_trans = I(CHL)
)

# Create an empty list to store plots
plots <- list()

# Iterate over traits and fit models accordingly
for (trait in traits_of_interest) {
  
  # Determine the trait name and formula
  formula <- as.formula(paste(trait, "~ treatment_mmol + (1 | species)"))
  
  # Fit the mixed-effects model using REML = FALSE for AIC comparison
  mod <- lmer(formula, data = traits, REML = FALSE)
  
  # Get model predictions
  traits$predicted <- predict(mod)
  
  # Create ggplot for the untransformed trait with log10-scaled axes
  p <- ggplot(traits, aes(x = treatment_mmol, y = !!sym(trait), color = species)) +
    geom_point(aes(alpha = 0.5, shape = species), size = 2) +
    geom_line(aes(y = predicted), linewidth = 1) +  # Model predictions
    scale_y_log10() +  # Log10 scale for y-axis
    #scale_x_continuous(trans = "log10") +  # Log10 scale for x-axis (if needed)
    labs(
      x = label_units[["treatment_mmol"]],  # Nitrogen in mM
      y = label_units[[trait]],  # Use formatted LaTeX labels for the trait
      color = "Species",
      shape = "Species"
    ) +
    scale_color_manual(values = josef_colors) +
    custom_theme +
    guides(
      alpha = "none",  # Remove alpha from the legend
      color = guide_legend(override.aes = list(alpha = 1, size = 4)),
      shape = guide_legend(override.aes = list(alpha = 1, size = 4))
    )
  
  # Store plots
  plots[[trait]] <- p
}

# Arrange all plots in a grid with a common legend at the bottom
final_plot <- ggarrange(
  plotlist = plots,
  labels = c("a", "b", "c"),
  nrow = 1,
  ncol = 3,
  common.legend = TRUE,
  legend = "bottom"  # Position the legend at the bottom
)

# Display the final plot
print(final_plot)

##  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Chat GPT4 used for assistance to print p values for mixed effects models
# across species, R2 values, and p values for simple linear models (lm, separate species)

# Create an empty list to store p-values
p_values <- list()

# Iterate over traits and fit models accordingly
for (trait in traits_of_interest) {
  
  # Determine the transformed trait name and formula
  transformed_trait <- paste0(trait, "_trans")
  formula <- as.formula(paste(transformed_trait, "~ treatment_mmol + (1 | species)"))
  
  # Fit the mixed-effects model using REML = FALSE for hypothesis testing
  mod <- lmer(formula, data = traits, REML = FALSE)
  
  # Extract p-value for treatment effect
  p_value <- anova(mod)$`Pr(>F)`[1]
  
  # Store in the list with trait name
  p_values[[trait]] <- p_value
}

# Print p-values to R console
print(p_values)

##  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# Load required packages
library(lme4)
library(MuMIn)

# Create empty lists to store R² values
r2_values <- list()

# Iterate over traits and fit models accordingly
for (trait in traits_of_interest) {
  
  # Determine the transformed trait name and formula
  transformed_trait <- paste0(trait, "_trans")
  formula <- as.formula(paste(transformed_trait, "~ treatment_mmol + (1 | species)"))
  
  # Fit the mixed-effects model using REML = FALSE for hypothesis testing
  mod <- lmer(formula, data = traits, REML = FALSE)
  
  # Extract R² values using MuMIn package
  r2_result <- r.squaredGLMM(mod)
  r2_marginal <- r2_result[1]  # Marginal R² (fixed effects only)
  r2_conditional <- r2_result[2]  # Conditional R² (fixed + random effects)
  
  # Store in list with trait name
  r2_values[[trait]] <- list(
    "Marginal R²" = r2_marginal,
    "Conditional R²" = r2_conditional
  )
}

# Print R² values to R console
print(r2_values)

##  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# Create an empty list to store p-values for each species and trait
species_p_values <- list()

# Iterate over traits
for (trait in traits_of_interest) {
  
  # Subset the data for the current transformed trait
  transformed_trait <- paste0(trait, "_trans")
  
  # Initialize a nested list for the current trait
  species_p_values[[trait]] <- list()
  
  # Iterate over each species
  for (sp in unique(traits$species)) {
    
    # Subset data for the current species
    species_data <- traits %>% filter(species == sp)
    
    # Fit a simple linear model for the current species and trait
    formula <- as.formula(paste(transformed_trait, "~ treatment_mmol"))
    lm_model <- lm(formula, data = species_data)
    
    # Extract the p-value for the slope (treatment_mmol)
    p_value <- summary(lm_model)$coefficients["treatment_mmol", "Pr(>|t|)"]
    
    # Store the p-value in the nested list
    species_p_values[[trait]][[as.character(sp)]] <- p_value
  }
}

# Print p-values for each species and trait
print(species_p_values)
