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
                         levels = c("R. sativus", "B. officinalis", "H. vulgare"))

# Compute growth rate per day
growth_period_days <- 6 * 7  # 6-week experiment
traits$GRT <- traits$dry_whole_g / growth_period_days

# Filter only by traits of interest
traits <- traits %>% select(species, LDMC, LMA, CHL, GRT, treatment_mmol)

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
  "LDMC" = expression(1 / LDMC ~ "(g mg"^{-1} * ")"),  # LDMC as fraction
  "LMA" = expression(1 / LMA ~ "(g m"^{-2} * ")"),  # LMA as fraction
  "CHL" = expression(1 / sqrt(CHL) ~ "(" * mu * "g"^{-1/2} * " cm)"),  # CHL with sqrt transformation in both value and units
  "treatment_mmol" = "Nitrogen (mM)" # Treatment nitrogen
)

# Import Data ----
# REMINDER: Set Working Directory -> Source File Location
# Define traits of interest
traits_of_interest <- c("LMA", "LDMC", "CHL")

# Apply specific transformations for each trait
traits <- traits %>%
  mutate(
    LMA_trans = 1 / LMA,
    LDMC_trans = 1 / LDMC,
    CHL_trans = 1 / sqrt(CHL)
  )

# Create an empty list to store plots
plots <- list()

# Iterate over traits and fit models accordingly
for (trait in traits_of_interest) {
  
  # Determine the transformed trait name and formula
  if (trait %in% c("LMA", "LDMC")) {
    transformed_trait <- paste0(trait, "_trans")
    formula <- as.formula(paste(transformed_trait, "~ treatment_mmol + (1 | species)"))
  } else if (trait == "CHL") {
    transformed_trait <- "CHL_trans"
    formula <- as.formula(paste(transformed_trait, "~ treatment_mmol + (treatment_mmol | species)"))
  }
  
  # Fit the mixed-effects model using REML = FALSE for AIC comparison
  mod <- lmer(formula, data = traits, REML = FALSE)
  
  # Get model predictions
  traits$predicted <- predict(mod)
  
  # # Extract R² using MuMIn
  # r2_values <- r.squaredGLMM(mod)
  # r2_marginal <- round(r2_values[1], 3)  # Marginal R² (fixed effects)
  # r2_conditional <- round(r2_values[2], 3)  # Conditional R² (fixed + random)
  # p_value <- anova(mod)$`Pr(>F)`[1]  # P-value
  # #model_AIC <- AIC(mod)  # Extract AIC
  # 
  # # Create annotation text
  # annotation_text <- paste0("P value: ", round(p_value, digits = 3))
  #                           #"AIC = ", round(model_AIC, digits = 2))
  
  # # Define placement coordinates dynamically
  # x_position <- max(traits$treatment_mmol) * 0.95  # Near the left edge
  # y_position <- max(traits[[paste0(trait, "_trans")]]) * 0.995  # Adjust for top placement
  
  p <- ggplot(traits, aes(x = treatment_mmol, y = !!sym(paste0(trait, "_trans")), color = species)) +
    geom_point(aes(alpha = 0.5, shape = species), size = 2) +
    geom_line(aes(y = predicted), linewidth = 1) +  # Model predictions
    scale_y_continuous() +
    labs(
      x = label_units[["treatment_mmol"]],  # Reciprocal nitrogen
      y = label_units[[trait]],  # Use formatted LaTeX labels
      color = "Species",
      shape = "Species"
    ) +
    scale_color_manual(values = josef_colors) +
    custom_theme +
    # annotate("text", x = x_position, y = y_position,
    #          label = annotation_text, hjust = 1, vjust = 1, size = 3.5) +
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

print(final_plot)
