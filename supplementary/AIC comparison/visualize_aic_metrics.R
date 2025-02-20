##' Create lollipop plots of delta_AIC for
##' model selection.
##'
##' @author [Nathan D Malamud]
##' @date [2025-01-28]
##' 
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
aic_results <- read.csv("aic_results.csv")

# Visualizing model comparison ----
# Define traits of interest
traits_of_interest <- c("LMA", "LDMC", "area_cm2", "CHL", "GRT")

# Initialize an empty list to store plots
plots <- list(
  "R. sativus" = list(),
  "B. officinalis" = list(),
  "H. vulgare" = list()
)

model_species <- unique(aic_results$species)

for (species_name in model_species) {
  for (trait_name in traits_of_interest) {
    
    #  Extract data for trait and species
    aic_compare <- aic_results %>% filter(species == species_name,
                                          trait == trait_name)
    
    # Calculate weighted AIC values from raw AICs, normalize to sum
    aic_compare$deltaAIC = aic_compare$AIC_val - min(aic_compare$AIC_val) 
    aic_compare$deltaAIC_norm = aic_compare$AIC_val / sum(aic_compare$AIC_val)
    
    # Create lollipop chart for this species and trait
    lollipop_plot <- ggplot(aic_compare, aes(x = model, y = deltaAIC)) +
      # Add the "head" and "stick" for the lollipop
      geom_segment(aes(x = model, xend = model, y = 0, yend = deltaAIC), color = "gray") +
      geom_point(size=2.0, color = "black") +
      custom_theme +
      ggtitle(paste(label_units[[trait_name]])) # Add title for each plot
    
    # Add the plot to the list
    plots[[species_name]] <- append(plots[[species_name]], list(lollipop_plot))
  }
}

# Combine plots for each species into a single row with a common legend
p1 <- ggarrange(plotlist = plots[["R. sativus"]], nrow = 1)
p2 <- ggarrange(plotlist = plots[["B. officinalis"]], nrow = 1)
p3 <- ggarrange(plotlist = plots[["H. vulgare"]], nrow = 1)

# Combine all rows into a single plot
# Add species labels to the side
# Combine rows into a single plot
final_plot <- ggarrange(
  p1, p2, p3,
  nrow = 3,
  labels = NULL # Remove default labels
)

# Add species labels to the side manually
annotated_plot <- ggarrange(
  ggarrange(
    text_grob("RADISH (n = 39)", rot = 90, size = 12, face="bold"),
    text_grob("BORAGE (n = 37)", rot = 90, size = 12, face="bold"),
    text_grob("BARLEY (n = 21)", rot = 90, size = 12, face="bold"),
    nrow = 3, ncol = 1, widths = c(1) # Create a column for species labels
  ),
  final_plot,
  ncol = 2, widths = c(0.1, 1) # Adjust widths for label and plot
)

# Display the annotated plot
print(annotated_plot)
