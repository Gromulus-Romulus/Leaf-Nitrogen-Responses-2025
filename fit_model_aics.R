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
  "Phi_PS2" = "PSII",
  "CHL" = "CHL",
  "GRT" = "GRT",
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

# Model validation specifications ----
# Define helper function for log-model fitting ----
#   - Filters dataframe for one trait and one species
#   - Fits multiple transformation of variables
#   - Compares AIC values
fit_log_models <- function(df, var, species_filter) {
  
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
  fit5 <- lm(sqrt(var) ~ treatment_mmol, data = df)
  fit6 <- lm(var ~ sqrt(treatment_mmol), data = df)
  fit7 <- lm(sqrt(var) ~ sqrt(treatment_mmol), data=df)
  
  # Extract model statistics (mini-function)
  model_stats <- function(model) {
    data.frame(
      AIC = AIC(model),
      R2 = summary(model)$r.squared,
      Adj_R2 = summary(model)$adj.r.squared,
      P_Value = summary(model)$coefficients[2, 4]
    )
  }
  
  # ... get stats for each model
  stats1 <- model_stats(fit1)
  stats2 <- model_stats(fit2)
  stats3 <- model_stats(fit3)
  stats4 <- model_stats(fit4)
  stats5 <- model_stats(fit5)
  stats6 <- model_stats(fit6)
  stats7 <- model_stats(fit7)
  
  # Combine results into a single dataframe
  model_df <- rbind(
    cbind(Model = "X ~ N", stats1),
    cbind(Model = "log(X) ~ N", stats2),
    cbind(Model = "X ~ log(N)", stats3),
    cbind(Model = "log(X) ~ log(N)", stats4),
    cbind(Model = "sqrt(X) ~ N", stats5),
    cbind(Model = "X ~ sqrt(N)", stats6),
    cbind(Model = "sqrt(X) ~ sqrt(N)", stats7)
  )
  
  model_df$Model <- factor(model_df$Model,
                           levels = c(
                             "X ~ N", "log(X) ~ N", "X ~ log(N)", "log(X) ~ log(N)",
                             "X ~ sqrt(N)", "sqrt(X) ~ N", "sqrt(X) ~ sqrt(N)"))
  
  return(model_df)
}

# Initialize empty nested structure
model_results <- list() 

# Loop fits regression models to traits across species
for (species in levels(traits$species)) {
  # Fit models to separate species for traits of interest
  lma <- fit_log_models(traits, "LMA", species_filter = species)
  ldmc <- fit_log_models(traits, "LDMC", species_filter = species)
  acm2 <- fit_log_models(traits, "area_cm2", species_filter = species)
  chl <- fit_log_models(traits, "CHL", species_filter = species)
  psii <- fit_log_models(traits, "Phi_PS2", species_filter = species)
  grt <- fit_log_models(traits, "GRT", species_filter = species)
  
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

# Visualizing model comparison ----

# Define traits of interest
traits_of_interest <- c("LMA", "LDMC", "area_cm2", "CHL", "Phi_PS2", "GRT")

# Initialize an empty list to store plots
plots <- list(
  "R. sativus" = list(),
  "B. officinalis" = list(),
  "H. vulgare" = list()
)

for (species in levels(traits$species)) {
  for (trait in traits_of_interest) {
    
    # Extract model stats from nested dataframe
    model_stats <- model_results[[species]][[trait]]
    
    # Obtain AIC value for models
    aic_compare <- model_stats %>% select(Model, AIC)
    
    # Create lollipop chart for this species and trait
    lollipop_plot <- ggplot(aic_compare, aes(x = Model, y = AIC)) +
      # Add the "head" and "stick" for the lollipop
      geom_segment(aes(x = Model, xend = Model, y = 0, yend = AIC), color = "gray") +
      geom_point(size=2.0, color = "black") +
      custom_theme +
      ggtitle(paste(label_units[[trait]])) # Add title for each plot
    
    # Add the plot to the list
    plots[[species]] <- append(plots[[species]], list(lollipop_plot))
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

# Save model results as CSV ----
# Initialize an empty data frame with the desired column names
model_frame <- data.frame(
  trait = character(),
  species = character(),
  model = character(),
  AIC = numeric(),
  R2 = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

# Collect model results for each species and trait
for (species in levels(traits$species)) {
  for (trait in traits_of_interest) {
    # Extract model statistics for this species and trait
    model_stats <- model_results[[species]][[trait]]
    
    # Create a temporary data frame for this species-trait combination
    temp_frame <- data.frame(
      trait = trait,
      species = species,
      model = model_stats$Model,
      AIC = model_stats$AIC,
      R2 = model_stats$R2,
      p_val = model_stats$P_Value,
      stringsAsFactors = FALSE
    )
    
    # Combine with the main data frame
    model_frame <- rbind(model_frame, temp_frame)
  }
}

# Save the model results to a CSV file
write.csv(model_frame, "aic_results.csv", row.names = FALSE)

# View the result
head(model_frame)
