##' Compare multiple models of nitrogen saturation
##' for plant functional traits
##' 
##' @author [Nathan Malamud]
##' @date [2025-01-23]
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

# Import Data ----
# REMINDER: Set Working Directory -> Source File Location
# Define factor levels as species
traits <-  read_csv("./data/traits.csv")
traits$species <- factor(traits$species,
                         levels=c("R. sativus",
                                  "B. officinalis",
                                  "H. vulgare"))

# Styles ----
# Calculate rate of growth
growth_period_days <- 6 * 7 # 6 week experiment
traits$GRT <- (traits$dry_whole_g / growth_period_days)

# Filter by metrics of interest only
traits <- traits %>%
  select(barcodeID, species, treatment_mmol,
         LDMC, LMA, area_cm2, CHL, Phi_PS2, GRT)

# Define consistent font size
base_font_size <- 12

# Define a custom theme for all plots
custom_theme <- theme_classic() +  # Base theme
  theme(
    # Text and font styling
    text = element_text(family = "sans", size = base_font_size),
    axis.text = element_text(size = base_font_size),  
    legend.text = element_text(size = base_font_size),
    
    # Title - left justified, 12 point font, bold
    plot.title = element_text(size = base_font_size, hjust = 0.5, face = "bold"),
    
    # Axis labels
    axis.title.x = element_text(size = base_font_size),
    axis.title.y = element_text(size = base_font_size, margin = margin(r = 10)),  # Margin for y-axis title
    
    # Panel and grid styling
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.1),  
    panel.grid.minor = element_blank(),  # No minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # White background
    
    # Aspect ratio
    aspect.ratio = 1  # 1:1 ratio
  )

# Custom colors advised by J. Garen
josef_colors <- c("R. sativus" = "#299680",
                  "B. officinalis" = "#7570b2",
                  "H. vulgare" = "#ca621c")

# Define units for variables
# TODO: format with Latex expressions
label_units <- c(
  "LDMC" = "LDMC (mg / g)",
  "LMA" = "LMA (g / m²)",
  "area_cm2" = "Leaf Area (cm²)",
  "Phi_PS2" = "PSII Fraction",
  "CHL" = "CHL (ug / cm²)",
  "GRT" = "Growth Rate (g / day)",
  "treatment_mmol" = "N (mM)"
)

# Model validation specifications ----
#   For each trait, regress onto nitrogen treatment.
#   Produce species x method table containing AIC values.
#   Produce PDF of model fits.

# Helper function
#   - Selects variable X and nitrogen treatment N
#   - Fits X ~ N, log(X) ~ N, X ~ log(N), log(X) ~ log(N)
#   - Returns (df, plot)
#   - df: data frame with model fits
#   - plot: ggplot2 object with model fits
fit_models <- function(df, var, species_filter) {
  library(dplyr)
  
  # Filter for the selected species
  df <- df %>%
    filter(species == species_filter) %>%
    select(species, treatment_mmol, all_of(var))
  
  # Rename variable for compatibility in models
  colnames(df)[colnames(df) == var] <- "var"
  
  # Fit models
  fit1 <- lm(var ~ treatment_mmol, data = df)
  fit2 <- lm(log(var + 0.01) ~ treatment_mmol, data = df)
  fit3 <- lm(var ~ log(treatment_mmol + 0.01), data = df)
  fit4 <- lm(log(var + 0.01) ~ log(treatment_mmol + 0.01), data = df)
  
  # Extract model statistics
  model_stats <- function(model) {
    data.frame(
      AIC = AIC(model),
      R2 = summary(model)$r.squared,
      Adj_R2 = summary(model)$adj.r.squared,
      P_Value = summary(model)$coefficients[2, 4]
    )
  }
  
  stats1 <- model_stats(fit1)
  stats2 <- model_stats(fit2)
  stats3 <- model_stats(fit3)
  stats4 <- model_stats(fit4)
  
  # Combine results into a single dataframe
  model_df <- data.frame(
    Model = c("X ~ N", "log(X) ~ N", "X ~ log(N)", "log(X) ~ log(N)"),
    stats1, stats2, stats3, stats4
  )
  
  model_df <- rbind(
    cbind(Model = "X ~ N", stats1),
    cbind(Model = "log(X) ~ N", stats2),
    cbind(Model = "X ~ log(N)", stats3),
    cbind(Model = "log(X) ~ log(N)", stats4)
  )
  
  return(model_df)
}

## Fit model comparisons ----
# TODO: iterate through species of interes

# Initialize empty nested structure
model_results <- list() 

for (species in levels(traits$species)) {
  # Fit models to separate species for traits of interest
  lma <- fit_models(traits, "LMA", species_filter = species)
  ldmc <- fit_models(traits, "LDMC", species_filter = species)
  acm2 <- fit_models(traits, "area_cm2", species_filter = species)
  chl <- fit_models(traits, "CHL", species_filter = species)
  psii <- fit_models(traits, "Phi_PS2", species_filter = species)
  grt <- fit_models(traits, "GRT", species_filter = species)
  
  # Combine results as entry in nested structure
  model_results[[species]] <- list(
    LMA = lma,
    LDMC = ldmc,
    area_cm2 = acm2,
    CHL = chl,
    Phi_PS2 = psii,
    GRT = grt
  )
}

## Visualizing model comparison ----

# Define traits of interest
traits_of_interest <- c("LMA", "LDMC", "area_cm2", "CHL", "Phi_PS2", "GRT")

# Open PDF for output
pdf("model_fit_comparison.pdf", width = 11, height = 7) # Landscape orientation

# Iterate through each trait
for (trait in traits_of_interest) {
  
  # Subset data
  data <- traits %>% select(species, treatment_mmol, all_of(trait))
  
  # Plot the data for the current trait
  plot <- ggplot(data, aes_string(x = "treatment_mmol", y = trait)) +
    geom_point(size = 2, shape = 1) +
    # Linear model: X ~ N
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = "X ~ N")) +
    # Linear model: log(X) ~ N
    geom_smooth(method = "lm", formula = exp(log(y + 0.01)) ~ x, se = FALSE, aes(color = "log(X) ~ N")) +
    # Linear model: X ~ log(N)
    geom_smooth(method = "lm", formula = y ~ log(x + 0.01), se = FALSE, aes(color = "X ~ log(N)")) +
    # Linear model: log(X) ~ log(N)
    geom_smooth(method = "lm", formula = exp(log(y + 0.01)) ~ log(x + 0.01), se = FALSE, aes(color = "log(X) ~ log(N)")) +
    facet_wrap(~species) +
    labs(
      x = label_units[["treatment_mmol"]],
      y = label_units[[trait]],
      title = paste("Model Fits for", trait, "by Species"),
      color = "Model"
    ) +
    scale_color_brewer(palette = "Set1") + 
    custom_theme +  
    theme(legend.position = "bottom")
  
  # Define species order and initialize a list to store grobs
  species_list <- c("R. sativus", "B. officinalis", "H. vulgare")
  library(dplyr)
  library(gridExtra)
  
  # Initialize a data frame to store combined results
  combined_table <- data.frame()
  
  # Iterate through each species to extract model fits
  for (species in species_list) {
    # Extract model_df for the current species and specified trait
    model_df <- model_results[[species]][[trait]]
    
    # Create a summary table for AIC and p-values
    species_table <- model_df %>%
      mutate(
        AIC = round(AIC, 1),
        P_Value = signif(P_Value, 2)
      ) %>%
      select(Model, AIC, P_Value) %>%
      rename(
        Model_Formula = Model,
        !!paste(species, "AIC", sep = "_") := AIC,
        !!paste(species, "P_Value", sep = "_") := P_Value
      )
    
    # Combine the tables
    if (nrow(combined_table) == 0) {
      combined_table <- species_table
    } else {
      combined_table <- full_join(combined_table, species_table, by = "Model_Formula")
    }
  }
  
  # Print the combined table
  print(combined_table)
  
  # Convert the combined table to a grob for visualization
  final_combined <- tableGrob(
    combined_table,
    theme = ttheme_minimal(base_size = 10)
  )
  
  # Combine the plot and the table
  final_combined_layout <- plot / 
    wrap_elements(final_combined) + 
    plot_layout(heights = c(3, 1)) +  # Adjust height ratios
    plot_annotation(
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  # Display the combined layout
  print(final_combined_layout)
}

# Close and save the PDF
dev.off()
