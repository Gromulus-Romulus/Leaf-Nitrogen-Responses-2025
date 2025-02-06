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
  "LDMC" = expression(frac(1, LDMC) ~ "(g mg"^{-1} * ")"),  # LDMC as fraction
  "LMA" = expression(frac(1, LMA) ~ "(g m"^{-2} * ")"),  # LMA as fraction
  "CHL" = expression(frac(1, CHL) ~ "(" * mu * "g cm"^{-2} * ")"),  # CHL with micrograms
  "treatment_mmol" = expression(frac(1, N) ~ "(mM"^{-1} *")")  # Treatment as fraction
)

# Import Data ----
# REMINDER: Set Working Directory -> Source File Location
# Define factor levels as species
traits <-  read_csv("./data/traits.csv")

traits$species <- factor(traits$species,
                         levels=c("R. sativus",
                                  "B. officinalis",
                                  "H. vulgare"))

# Filter by metrics of interest only
traits <- traits %>%
  select(barcodeID, species, treatment_mmol,
         LDMC, LMA, CHL)

# Plot visual comparison of log-transformed y axis (or log models) ----

# As of 2025-01-29
# Exploratory model and consultation with literature indicate
# reciprocal fitting with random intercepts (not slopes) is appropriate.

# ------ Attempt fitting with numerical buffer -----
traits_of_interest <- c("LMA", "LDMC", "CHL")

# Create an empty list to store plots
plots <- list()

# Apply reciprocal transform to treatment
# Set guideline to 1% of smallest non-zero value

N_BUFF = 100 

traits <- traits %>%
  mutate(treatment_reciprocal = 1 / (treatment_mmol + N_BUFF))  # Avoid divide by zero

for (trait in traits_of_interest) {
  # Apply reciprocal transform to the response variable
  traits <- traits %>%
    mutate(!!sym(trait) := 1 / (!!sym(trait) + 1e-6))  # Avoid divide by zero
  
  # Fit the mixed-effects model (linear reciprocal response)
  formula <- as.formula(paste(trait, "~ treatment_reciprocal + (1 | species)"))
  mod <- lmer(formula, data = traits)
  
  # Get model predictions
  traits$predicted <- predict(mod)
  
  # Extract R² using MuMIn
  r2_values <- r.squaredGLMM(mod)
  r2_marginal <- round(r2_values[1], 3)  # Marginal R² (fixed effects)
  r2_conditional <- round(r2_values[2], 3)  # Conditional R² (fixed + random)
  model_AIC <- AIC(mod)  # Extract AIC
  
  # Extract model slope
  slope <- round(fixef(mod)["treatment_reciprocal"], 3)  # Extract slope coefficient
  
  # Create annotation text
  annotation_text <- paste0("Slope = ", slope, "\n",
                            "Marginal R² = ", r2_marginal, "\n",
                            "Conditional R² = ", r2_conditional, "\n",
                            "AIC = ", model_AIC)
  
  # Define placement coordinates dynamically
  x_position <- max(traits$treatment_reciprocal) * 0.9  # Near the left edge
  y_position <- max(traits[[trait]]) * 0.95  # Adjust for top placement
  
  p <- ggplot(traits, aes(x = treatment_reciprocal, y = !!sym(trait), color = species)) +
    geom_point(aes(alpha = 0.5, shape = species), size = 2) +
    geom_line(aes(y = predicted), linewidth = 1) +  # Model predictions
    scale_y_continuous() +
    labs(
      x = label_units[["treatment_mmol"]],  # Reciprocal nitrogen
      y = label_units[[trait]],  # Use formatted LaTeX labels
      color = "Species"
    ) +
    scale_color_manual(values = josef_colors) +
    custom_theme +
    annotate("text", x = x_position, y = y_position,
             label = annotation_text, hjust = 1, vjust = 1, size = 3.5)
  
  # Store plots
  plots[[trait]] <- p
}

# Arrange all plots in a grid
buffer_plot <- ggarrange(plotlist = plots,
          labels = c("a", "b", "c"),
          nrow = 1, ncol = 3)

# ------ Attempt fitting with filtering 0s -----
traits_of_interest <- c("LMA", "LDMC", "CHL")

# Create an empty list to store plots
plots <- list()

# Apply filter zeros
traits <- traits %>%
  filter(treatment_mmol > 0)  # Remove zero values 

for (trait in traits_of_interest) {
  # Apply reciprocal transform to the response variable
  traits <- traits %>%
    mutate(!!sym(trait) := 1 / (!!sym(trait) + 1e-6))  # Avoid divide by zero
  
  # Fit the mixed-effects model (linear reciprocal response)
  formula <- as.formula(paste(trait, "~ treatment_reciprocal + (1 | species)"))
  mod <- lmer(formula, data = traits)
  
  # Get model predictions
  traits$predicted <- predict(mod)
  
  # Extract R² using MuMIn
  r2_values <- r.squaredGLMM(mod)
  r2_marginal <- round(r2_values[1], 3)  # Marginal R² (fixed effects)
  r2_conditional <- round(r2_values[2], 3)  # Conditional R² (fixed + random)
  model_AIC <- AIC(mod)  # Extract AIC
  
  # Extract model slope
  slope <- round(fixef(mod)["treatment_reciprocal"], 3)  # Extract slope coefficient
  
  # Create annotation text
  annotation_text <- paste0("Slope = ", slope, "\n",
                            "Marginal R² = ", r2_marginal, "\n",
                            "Conditional R² = ", r2_conditional, "\n",
                            "AIC = ", model_AIC)
  
  # Define placement coordinates dynamically
  x_position <- max(traits$treatment_reciprocal) * 0.75  # Near the left edge
  y_position <- max(traits[[trait]]) * 0.95  # Adjust for top placement
  
  p <- ggplot(traits, aes(x = treatment_reciprocal, y = !!sym(trait), color = species)) +
    geom_point(aes(alpha = 0.5, shape = species), size = 2) +
    geom_line(aes(y = predicted), linewidth = 1) +  # Model predictions
    scale_y_continuous() +
    labs(
      x = label_units[["treatment_mmol"]],  # Reciprocal nitrogen
      y = label_units[[trait]],  # Use formatted LaTeX labels
      color = "Species"
    ) +
    scale_color_manual(values = josef_colors) +
    custom_theme +
    annotate("text", x = x_position, y = y_position,
             label = annotation_text, hjust = 1, vjust = 1, size = 3.5)
  
  # Store plots
  plots[[trait]] <- p
}

# Arrange all plots in a grid
filter_plot <- ggarrange(plotlist = plots,
                         labels = c("d", "e", "f"),
                         nrow = 1, ncol = 3)

# Combine plots
combined_plot <- plot_grid(buffer_plot, filter_plot, nrow = 2)
print(combined_plot)
