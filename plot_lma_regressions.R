# Visualize regression models

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
    text = element_text(family = "sans", size = 12),
    axis.text = element_text(size = 10),  
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    
    # Axis labels
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # Margin for y-axis title
    
    # Panel and grid styling
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.1),  
    panel.grid.minor = element_blank(),  # No minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # White background
    
    # Aspect ratio
    aspect.ratio = 1  # 1:1 ratio
  )

# Custom colors advised by J. Garen
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Define units for variables
# TODO: format with Latex expressions
label_units <- c(
  "CHL" = "CHL (ug / cm²)",
  "LDMC" = "LDMC (mg / g)",
  "LMA" = "LMA (g / m²)",
  "area_cm2" = "Area (cm²)",
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

# Visualize regression models ----
traits_of_interest <- c("LMA", "LDMC", "area_cm2", "CHL", "GRT")

# Open pdf for visual comparison
pdf("trait_geom_smooths.pdf", width = 11, height = 8.5)

# Loop through each trait of interest
for (trait in traits_of_interest) {
  
  # Extract data for the current trait
  trait_data <- traits %>% select(species, treatment_mmol, all_of(trait))
  
  # Define the list of model formulas
  model_formulas <- list(
    "linear" = as.formula("y ~ x"),
    "quadratic" = as.formula("y ~ x + poly(x, 2)"),
    "cubic" = as.formula("y ~ poly(x, 3)"),
    "log" = as.formula("y ~ log(x + 0.00001)"),
    "exp" = as.formula("y ~ exp(x)"),
    "sqrt" = as.formula("y ~ sqrt(x)")
  )
  
  # Store plots
  plots <- list()
  
  # Loop through each model
  for (model_name in names(model_formulas)) {
    
    # Create plot for the current model
    p <- ggplot(trait_data, aes(x = treatment_mmol, y = .data[[trait]], color = species)) +
      geom_point(alpha = 0.6) +
      stat_smooth(method = "lm", formula = model_formulas[[model_name]], se = FALSE) +
      scale_color_manual(values = josef_colors, name = "Species") +
      scale_shape_manual(values = c(16, 17, 18), name = "Species") +
      labs(
        title = paste(model_name),
        x = "Treatment (mmol)",
        y = trait
      ) +
      guides(color = guide_legend(override.aes = list(shape = c(16, 17, 18)))) + custom_theme
    
    # Add plot to the list
    plots[[model_name]] <- p
  }
  
  # Dynamically determine layout for combined plots
  ncol <- ceiling(sqrt(length(model_formulas)))
  nrow <- ceiling(length(model_formulas) / ncol)
  
  # Combine all model plots using ggarrange
  final_plot <- ggarrange(plotlist = plots,
                          ncol = ncol, nrow = nrow,
                          common.legend = TRUE)
  
  # Print to viewer
  print(final_plot)
}

# Close pdf
dev.off()

