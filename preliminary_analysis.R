# Author: Nathan Malamud
# Date: 2025.01.15

#Analysis Packages----
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(reshape2)
library(scales)
library(smatr)
library(GGally)
library(factoextra) # PCA analysis (TODO switch to RDA)

#Load data for traits----
# REMINDER: Set Working Directory -> Source File Location
traits <- tryCatch(
  read_csv("./data/traits.csv"),
  error = function(e) stop("Error loading traits.csv: ", e)
)

# Define Factor Levels (Treatment and Species)
traits$species <- factor(traits$species,
                         levels=c("R. sativus", "B. officinalis", "H. vulgare")
                         )

##Calculate rate of growth----
growth_period_days <- 6 * 7 # 6 week experiment
traits$GRT <- (traits$dry_whole_g / growth_period_days)

##Filter by traits of interest only----
traits <- traits %>%
  select(barcodeID, species, treatment_mmol,
         LDMC, LMA, CHL, Phi_PS2, GRT)

#Styling and aesthetics----
# Define a custom theme for all plots

##GGplot settings---
custom_theme <- theme_classic() +  # Start with a minimal theme
  theme(
    text = element_text(family = "sans", size = 12),  # Set font family and size
    axis.title.x = element_text(size = 12),  # Customize x-axis title
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # Add margin to y-axis title
    axis.text = element_text(size = 10),  # Customize axis text
    legend.text = element_text(size = 11),  # Customize legend text
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth=0.1),  # Customize major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and style the plot title
    aspect.ratio = 1  # Fix the aspect ratio (1:1)
  )

##Custom Colors----
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Units - # TODO: format with Latex expressions
label_units <- c(
   "CHL" = "CHL (ug / cm²)",         # Chlorophyll content in micrograms per square centimeter
   "LDMC" = "LDMC (mg / g)",         # Leaf dry matter content in milligrams per gram
   "LMA" = "LMA (g / m²)",
   "Phi_PS2" = "PSII Fraction",
   "treatment_mmol" = "N (mM)"
)

#Treatment Responses-----
# 1) Show response of LDMC, LMA, CHL, and Phi_PS2 to treatment_mmol

##LMA saturation----
lma_Nmm <- smatr::sma(LMA ~ treatment_mmol*species, data=traits, method="OLS")
summary(lma_Nmm)

lma_Nmm_plot <- ggplot(traits, aes(y = LMA, x = treatment_mmol)) +
  geom_point(size=2, alpha=0.5, aes(color=species, shape=species)) +
  stat_ma_line(aes(color = species), method="OLS", se=F) +
  scale_color_manual(values=josef_colors, name="Species") +
  scale_shape_manual(values=c(16, 17, 18), name="Species") +
  labs(x = label_units[["treatment_mmol"]], y = label_units[["LMA"]]) +
  ylim(0, NA) +  # Set minimum x-axis limit to 0, allow maximum to be auto-scaled
  theme(
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    aspect.ratio=1,
    legend.position="bottom"
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

##LDMC saturation----
ldmc_Nmm <- smatr::sma(LDMC ~ treatment_mmol*species, data = traits, method="OLS")
summary(ldmc_Nmm)

ldmc_Nmm_plot <- ggplot(traits, aes(y = LDMC, x = treatment_mmol)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +
  stat_ma_line(aes(color = species), method = "OLS", se = F) +
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  labs(x = label_units[["treatment_mmol"]], y = label_units[["LDMC"]]) +
  theme(
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    aspect.ratio=1,
    legend.position="bottom"
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

##CHL saturation----
chl_Nmm <- smatr::sma(CHL ~ treatment_mmol * species, data = traits, method="OLS")
summary(chl_Nmm)

chl_Nmm_plot <- ggplot(traits, aes(x = treatment_mmol, y = CHL)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +
  stat_ma_line(aes(color = species), method = "OLS", se = FALSE) +
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  labs(x = label_units[["treatment_mmol"]], y = label_units[["CHL"]]) +
  theme(
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    aspect.ratio=1,
    legend.position="bottom"
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

##Phi_PS2 saturation----
phi_ps2_Nmm <- smatr::sma(Phi_PS2 ~ treatment_mmol * species, data = traits, method="OLS")
summary(phi_ps2_Nmm)

phi_ps2_Nmm_plot <- ggplot(traits, aes(y = Phi_PS2, x = treatment_mmol)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +
  stat_ma_line(aes(color = species), method = "OLS", se = FALSE) +
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  labs(x = label_units[["treatment_mmol"]], y = label_units[["Phi_PS2"]]) +
  theme(
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    aspect.ratio=1,
    legend.position="bottom"
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

##Arrange figures with ggarrange and save pdf----
# Update each plot to move the legend to the bottom
lma_Nmm_plot <- lma_Nmm_plot + theme(legend.position = "bottom")
ldmc_Nmm_plot <- ldmc_Nmm_plot + theme(legend.position = "bottom")
chl_Nmm_plot <- chl_Nmm_plot + theme(legend.position = "bottom")
phi_ps2_Nmm_plot <- phi_ps2_Nmm_plot + theme(legend.position = "bottom")

# Combine plots into one row with legends at the bottom
figure2_saturation_plot <- ggarrange(
  lma_Nmm_plot, ldmc_Nmm_plot, chl_Nmm_plot, phi_ps2_Nmm_plot,
  ncol = 2, nrow = 2,  # Arrange in a single row
  labels = c("a", "b", "c", "d"),  # Add subplot labels
  common.legend = TRUE,  # Combine legends into one
  legend = "bottom"  # Place the combined legend at the bottom
)

# Print to console
print(figure2_saturation_plot)

# Save to figures directory
# TODO: way to include R2 in figures
ggsave(
  filename = "./figures/figure2_saturation_plots.tiff",
  plot = figure2_saturation_plot,
  width = 5, height = 5  # Adjust height to accommodate the legend
)

#Pairwise Regressions----
# 2) Show pairwise responses of LDMC, LMA, CHL, and Phi_PS2 using glms (log-log sma)

## LMA vs LDMC----
# Compare trait-trait scaling relationships using SMA (standard major axis)
# regression models.
# Milos: OLS vs SMA, sum of squared residuals is changing. R2 is the same, trendlines different.
#   Syntax: * does slope, + does slope and elevation
lma_ldmc_sma <- smatr::sma(LMA ~ LDMC * species, log = "XY",
                           method = "SMA",
                           data = traits)
summary(lma_ldmc_sma, method = "SMA")

lma_ldmc_plot <- ggplot(traits, aes(y = LDMC, x = LMA)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +  # Adjusted point size and opacity
  stat_ma_line(aes(color = species), method = "SMA", se=F) +  # SMA regression line
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  scale_y_log10() +  # Log scale for Y-axis
  scale_x_log10() +  # Log scale for X-axis
  labs(x = label_units[["LMA"]], y = label_units[["LDMC"]]) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    aspect.ratio = 1,
    legend.position = "bottom"  # Legend position at the bottom
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

## LDMC vs CHL----
ldmc_chl_sma <- smatr::sma(LDMC ~ CHL * species, log = "XY",
                           method = "SMA",
                           data = traits)
summary(ldmc_chl_sma, method = "SMA")

ldmc_chl_plot <- ggplot(traits, aes(y = LDMC, x = CHL)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +  # Adjusted point size and opacity
  stat_ma_line(aes(color = species), method = "SMA", se=F) +  # SMA regression line
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  scale_y_log10() +  # Log scale for Y-axis
  scale_x_log10() +  # Log scale for X-axis
  labs(x = label_units[["CHL"]], y = label_units[["LDMC"]]) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    aspect.ratio = 1,
    legend.position = "bottom"  # Legend position at the bottom
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

## LMA vs CHL----
lma_chl_sma <- smatr::sma(LMA ~ CHL * species, log = "XY",
                           method = "SMA",
                           data = traits)
summary(lma_chl_sma)

lma_chl_plot <- ggplot(traits, aes(y = LMA, x = CHL)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +  # Adjusted point size and opacity
  stat_ma_line(aes(color = species), method = "SMA", se=F) +  # SMA regression line
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  scale_y_log10() +  # Log scale for Y-axis
  scale_x_log10() +  # Log scale for X-axis
  labs(x = label_units[["CHL"]], y = label_units[["LMA"]]) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    aspect.ratio = 1,
    legend.position = "bottom"  # Legend position at the bottom
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

## Phi_PS2 vs LDMC----
phi_ps2_ldmc_sma <- smatr::sma(Phi_PS2 ~ LDMC * species, log="XY",
                               method = "SMA", data = traits)
summary(phi_ps2_ldmc_sma, method = "SMA")

phi_ps2_ldmc_plot <- ggplot(traits, aes(x = LDMC, y = Phi_PS2)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +
  stat_ma_line(aes(color = species), method = "SMA", se=F) +
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  labs(x = label_units[["LDMC"]], y = label_units[["Phi_PS2"]]) +  # Use proper units or descriptions
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    aspect.ratio = 1,
    legend.position = "bottom"
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

## Phi_PS2 vs LMA----
phi_ps2_lma_sma <- smatr::sma(Phi_PS2 ~ LMA * species, log="XY",
                              method = "SMA", data = traits)
summary(phi_ps2_lma_sma)

phi_ps2_lma_plot <- ggplot(traits, aes(x = LMA, y = Phi_PS2)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +
  stat_ma_line(aes(color = species), method = "SMA", se=F) +
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  labs(x = label_units[["LMA"]], y = label_units[["Phi_PS2"]]) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    aspect.ratio = 1,
    legend.position = "bottom"
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme


## Phi_PS2 vs CHL----
phi_ps2_chl_sma <- smatr::sma(Phi_PS2 ~ CHL * species, method = "SMA", data = traits)
summary(phi_ps2_chl_sma)

phi_ps2_chl_plot <- ggplot(traits, aes(x = CHL, y = Phi_PS2)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +
  stat_ma_line(aes(color = species), method = "SMA", se=F) +
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  labs(x = label_units[["CHL"]], y = label_units[["Phi_PS2"]]) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    aspect.ratio = 1,
    legend.position = "bottom"
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

## apply ggarrange and save as figure 3. ----
# Combine all plots
figure3_regression_plot <- ggarrange(
  lma_ldmc_plot, ldmc_chl_plot, lma_chl_plot,
  phi_ps2_lma_plot, phi_ps2_ldmc_plot, phi_ps2_chl_plot,
  ncol = 3, nrow = 2,  # 2 rows, 3 columns
  labels = c("a", "b", "c", "d", "e", "f"),  # Subplot labels
  common.legend = TRUE,  # Combine legends
  legend = "bottom"  # Legend at the bottom
)

# Print to console
print(figure3_regression_plot)

# Save the figure
ggsave(
  filename = "./figures/figure3_regression_plots.tiff",
  plot = figure3_regression_plot,
  width = 10, height = 7.5  # Adjust width and height for layout
)
 
#PCA Analysis---
# 3) Show pairwise regression plot of variables 
# Fit PCA model
# PCA Analysis ---
# Select relevant variables and scale data
# PCA Analysis ---
# Select relevant variables and scale data
pca_data <- traits %>% select(species, LMA, LDMC, CHL, Phi_PS2, GRT)

# Fit PCA model
pca <- prcomp(pca_data %>% select(-species), center = TRUE, scale. = TRUE)

# Extract explained variation
explained_variation <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)  # Calculate % variance explained by each PC

# PCA Biplot
pca_plot <- fviz_pca_biplot(
  pca,
  geom.ind = "point",  # Plot individuals as points
  col.ind = pca_data$species,  # Color individuals by species
  palette = josef_colors,  # Custom species colors
  addEllipses = TRUE,  # Add confidence ellipses for each group
  ellipse.level = 0.95,  # Set confidence level for ellipses
  legend.title = "Species",  # Legend title
  repel = TRUE  # Avoid overlap of labels
) +
  labs(
    title = NULL,
    x = paste0("PC1 (", explained_variation[1], "% variance)"),
    y = paste0("PC2 (", explained_variation[2], "% variance)")
  ) +
  custom_theme +
  theme(
    legend.position = "bottom",  # Position legend below the plot
    aspect.ratio = 1  # Equal aspect ratio for x and y axes
  )

# Print the PCA plot
print(pca_plot)

# Save the PCA plot as Figure 4
ggsave(
  filename = "./figures/figure4_pca_plot.tiff",
  plot = pca_plot,
  width = 5, height = 5  # Adjust dimensions to fit biplot
)
