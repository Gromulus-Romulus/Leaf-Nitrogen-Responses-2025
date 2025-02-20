# Compare pairwise combinations of variable transforms
# for trait X and nitrogen treatment N.

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


# Import Data ----
# REMINDER: Set Working Directory -> Source File Location
# Define factor levels as species
traits <-  read_csv("./data/traits.csv")

traits$species <- factor(traits$species,
                         levels=c("R. sativus",
                                  "B. officinalis",
                                  "H. vulgare"))


# Filter and apply transforms ----
traits_of_interest <- c("LMA", "LDMC", "area_cm2", "CHL", "GRT")

# Calculate rate of growth
growth_period_days <- 6 * 7 # 6 week experiment
traits$GRT <- (traits$dry_whole_g / growth_period_days)

traits <- traits %>% select(species,
                            treatment_mmol,
                            all_of(traits_of_interest))

# Numerical buffer to avoid errors for log and inverse functions
NUM_BUFFER = 0.000001

# Define transform "link" functions
link_formulas <- list(
  "identity" = y ~ x,
  "log" = y ~ log(x + NUM_BUFFER),
  "sqrt" = y ~ sqrt(x),
  #"logit" = y ~ log(x / (1 - x)), Not working right now
  "inverse" = y ~ 1 / (x + NUM_BUFFER),
  "inverse_squared" = y ~ 1 / (x + NUM_BUFFER)^2
)

label_units <- c(
  "CHL" = "ug per cm2",
  "LDMC" = "mg per g",
  "LMA" = "g per m2",
  "GRT" = "g per day",
  "area_cm2" = "cm2",
  "treatment_mmol" = "mM"
)

# Apply link formulas to each trait of interest ----
linked_traits <- list()

for (trait in traits_of_interest) {
  for (link_name in names(link_formulas)) {
    link_formula <- link_formulas[[link_name]]
    
    # Apply the transformation
    transformed_data <- traits %>%
      select(all_of(trait)) %>%
      rename(x = all_of(trait)) %>%
      mutate(y = eval(link_formula[[3]])) %>%
      rename(!!trait := y) %>%
      mutate(trait = trait, transformation = link_name)
    
    # Store in list
    linked_traits[[paste(trait, link_name, sep = "_")]] <- transformed_data
  }
}

# Combine all transformations into one data frame
transformed_traits <- bind_rows(linked_traits)

# Add units for each trait
transformed_traits <- transformed_traits %>%
  mutate(units = label_units[trait])

# Melt dataframe ----
melted_traits <- transformed_traits %>%
  select(species, trait, treatment_mmol, transform, units, value)

# Produce box plots of each transformed trait ---

# TODO: open pdf
species_list <- unique(melted_traits$species)

trait <- "LMA"

for (trait in traits_of_interest) {
  for (species in species_list) {
    
    # Filter data for the current species and trait
    plot_data <- melted_traits %>%
      filter(species == species, trait == trait)
    
    # Generate box plot
    box_plot <- ggplot(plot_data, aes(x = as.factor(treatment_mmol), y = value, fill = species)) +
      geom_boxplot() +
      facet_grid(transform ~ species, scales = "free") +
      labs(
        title = paste(trait, ":", label_units[[trait]]),
        x = "Nitrogen Treatment (mmol)", fill = "Species"
      ) +
      scale_fill_manual(values = josef_colors) +
      custom_theme
    
    # Print to save in PDF
    print(box_plot)
  }
}

# TODO: close pdf

