##' Create small multiples of different structural
##' and biochemical traits responding to nitrogen.
##' 
##' @author: Nathan D. Malamud
##' @date: 2021-11-11
##' 

library(tidyverse)
library(ggplot2)
library(ggpubr)

# Load data for traits
# REMINDER: Set Working Directory -> Source File Location
data <- read_csv("./data/traits.csv")

# TODO: move prospect_rtm_output to data
# Load prospect measurements from spec curves
prospect <- read_csv("./data/molecular_content.csv") %>%
  select(-c(sampleID, species, treatment_mmol))
data <- merge(data, prospect, by = "barcodeID")

# Define Factor Levels (Treatment and Species)
data$treatment_mmol <- data$treatment_mmol
data$species <- factor(data$species, levels=c("R. sativus", "B. officinalis", "H. vulgare"))
data$LDMC <- as.numeric(data$LDMC)

# Define a custom theme for all plots
custom_theme <- theme_minimal(base_family = "sans") + 
  theme(
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25, color = "grey80"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    strip.text = element_text(hjust = 0, size = 10, face = "italic"),
    legend.position = "none"
  )

# Assign factor levels
data$species <- factor(data$species, levels = c("R. sativus", "B. officinalis", "H. vulgare"))

# Ensure numeric values (TODO: quality check all of this)
data$treatment_mmol <- as.numeric(data$treatment_mmol) 
data$Qamb <- as.double(data$Qamb)

# Define colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Define average growth rate for each species, assume 6 week interval
data$shoot_growth <- data$dry_whole_g / (6 * 7) # gram per day
data$CHL_cm2 <- data$CHL
data$CHL_g <- data$CHL_cm2 / (data$LMA) # TODO: check your units here
data$ETR <- (data$ETR / data$Qamb) * 50.0 # Based on PhiPSII measurements

# Filter out traits of interest
data <- data %>% select(barcodeID, species, treatment_mmol,
  LDMC, # Leaf dry matter content
  LMA, # Leaf mass per area
  CHL, # Chlorophyll pigment content
  Phi_PS2,
  dry_whole_g # Dry weight of whole plant
)

# TODO: decide whether log axis transformation is appropriate
# And adjust labels as such
# Log-transformed data
# data_log10 <- data %>%
#   mutate(
#     LMA = log10(LMA),
#     LDMC = log10(LDMC),
#     N = log10(N),
#     #ETR = log10(ETR),
#     CHL = log10(CHL),
#     dry_whole_g = log10(dry_whole_g)
#   )

# Define the custom theme with a fixed aspect ratio
custom_theme <- theme_classic() +  # Start with a minimal theme
  theme(
    text = element_text(family = "sans", size = 12),  # Set font family and size
    axis.title.x = element_text(size = 12),  # Customize x-axis title
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # Add margin to y-axis title
    axis.text = element_text(size = 10),  # Customize axis text
    legend.text = element_text(size = 11),  # Customize legend text
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth=0.2),  # Customize major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and style the plot title
    aspect.ratio = 1  # Fix the aspect ratio (1:1)
  )

# Loop through each variable and create plots
# Define variables of interest
variables_of_interest <- c("CHL", "LDMC", "LMA")

# Define y-axis labels with units
label_units <- c(
  CHL = "CHL (μg/cm²)",         # Chlorophyll content in micrograms per square centimeter
  LDMC = "LDMC (mg/g)",         # Leaf dry matter content in milligrams per gram
  LMA = "LMA (g/m²)"            # Leaf mass per area in grams per square meter
)

# Initialize an empty list to store plots
plot_list <- list()

# Loop through each variable and create plots
for (variable in variables_of_interest) {
  # Filter data for the current variable
  data_filtered <- data %>%
    select(barcodeID, species, treatment_mmol, !!sym(variable)) %>%
    rename(value = !!sym(variable))
  
  # Create the plot with points and LOESS curves
  plot <- ggplot(data_filtered, aes(x = treatment_mmol, y = value, color = species, shape = species)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE, linetype = "solid") +  # LOESS curve in black
    scale_color_manual(values = josef_colors, name=NULL) +
    scale_shape_manual(values = c(16, 17, 18), name=NULL) +  # Set shapes for species
    labs(x = "Nitrogen (mM)", y = label_units[[variable]]) +
    custom_theme
  
  # Add the plot to the list
  plot_list[[variable]] <- plot
}

# Arrange the plots using ggarrange and assign to combined_plot
combined_plot <- ggarrange(
  plotlist = plot_list,
  nrow = 1,
  labels = c("A", "B", "C"),  # Add labels to each plot
  common.legend = TRUE,
  legend = "bottom"
)

# Display the combined plot
print(combined_plot)
