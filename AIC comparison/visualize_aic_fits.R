##' Import stats_csv data and visualize summed AIC
##'
##' @author [Nathan D. Malamud]
##' @date [2025-01-30]
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
    text = element_text(family = "sans", size = 12),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),  
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    
    # Axis labels
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
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

traits_of_interest <- c("LMA", "LDMC", "CHL")

# Load stats csv file
stats_csv <- read_csv("./aic_model_stats.csv") %>%
  filter(Trait %in% traits_of_interest) %>%
  filter(Model %in% c("Random Intercept", "Random Slope"))
  #filter(Singular == FALSE)

# Sum AIC values by model formula
sum_stats_csv <- stats_csv %>%
  select(Trait, Model, Formula, AIC) %>%
  group_by(Trait, Formula) %>%
  summarise(Sum_AIC = sum(AIC))

# Use with ball-and-stick to visualize model fits for each trait and formula
# - y axis: sum_AIC, x axis: formula, facet by: trait
# Ball-and-stick plot: AIC on x-axis, Formula on y-axis
ggplot(sum_stats_csv, aes(x = Sum_AIC, y = Formula)) +
  geom_segment(aes(xend = Sum_AIC, yend = Formula, x = min(Sum_AIC)), color="grey", size = 0.5) +  # Stick
  geom_point(color = "black", size = 2) +  # Ball
  facet_grid(Trait~., scales = "free_y") +
  labs(
    title = "Summed AIC Values for Different Model Formulas",
    x = "Summed AIC",
    y = "Model Formula"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_manual(values = josef_colors) +
  custom_theme

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Now, let's do the same thing compared to R2 fit
avgR2_stats_csv <- stats_csv %>%
  select(Trait, Model, Formula, R2_conditional) %>%
  group_by(Trait, Formula) %>%
  summarise(mean_R2 = mean(R2_conditional))

# Use with ball-and-stick to visualize model fits for each trait and formula
# - y axis: sum_AIC, x axis: formula, facet by: trait
# Ball-and-stick plot: AIC on x-axis, Formula on y-axis
ggplot(avgR2_stats_csv, aes(x = mean_R2, y = Formula)) +
  geom_segment(aes(xend = mean_R2, yend = Formula, x = min(mean_R2)), color="grey", size = 0.5) +  # Stick
  geom_point(color = "black", size = 2) +  # Ball
  facet_grid(Trait ~ ., scales = "free_y") +
  labs(
    title = "Mean R2 Values for Different Model Formulas",
    x = "Mean R2",
    y = "Model Formula"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_manual(values = josef_colors) +
  custom_theme

