##' Produce preliminary figures for thesis proposal.
##'
##' @author [Nathan Malamud]
##' @date [2025-01-]

# Libraries ----
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(scales)
library(smatr)
library(factoextra)

# Import Data ----
# REMINDER: Set Working Directory -> Source File Location
# Define factor levels as species
traits <-  read_csv("./data/traits.csv")
traits$species <- factor(traits$species,
                         levels=c("R. sativus", "B. officinalis", "H. vulgare"))

# Calculate rate of growth
growth_period_days <- 6 * 7 # 6 week experiment
traits$GRT <- (traits$dry_whole_g / growth_period_days)

# Filter by metrics of interest only
traits <- traits %>%
  select(barcodeID, species, treatment_mmol,
         LDMC, LMA, CHL, area_cm2, GRT)

# Log-Log Trait Regressions ----
# 2) Show pairwise responses of LDMC, LMA, CHL, and area_cm2 using glms (log-log sma)

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

## area_cm2 vs LDMC----
leaf_area_ldmc_sma <- smatr::sma(area_cm2 ~ LDMC * species, log="XY",
                               method = "SMA", data = traits)
summary(leaf_area_ldmc_sma, method = "SMA")

leaf_area_ldmc_plot <- ggplot(traits, aes(x = LDMC, y = area_cm2)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +
  stat_ma_line(aes(color = species), method = "SMA", se=F) +
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  labs(x = label_units[["LDMC"]], y = label_units[["area_cm2"]]) +  # Use proper units or descriptions
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

## area_cm2 vs LMA----
leaf_area_lma_sma <- smatr::sma(area_cm2 ~ LMA * species, log="XY",
                              method = "SMA", data = traits)
summary(leaf_area_lma_sma)

leaf_area_lma_plot <- ggplot(traits, aes(x = LMA, y = area_cm2)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +
  stat_ma_line(aes(color = species), method = "SMA", se=F) +
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  labs(x = label_units[["LMA"]], y = label_units[["area_cm2"]]) +
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


## area_cm2 vs CHL----
leaf_area_chl_sma <- smatr::sma(area_cm2 ~ CHL * species, method = "SMA", data = traits)
summary(leaf_area_chl_sma)

leaf_area_chl_plot <- ggplot(traits, aes(x = CHL, y = area_cm2)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +
  stat_ma_line(aes(color = species), method = "SMA", se=F) +
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  labs(x = label_units[["CHL"]], y = label_units[["area_cm2"]]) +
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

# Apply ggarrange and save as figure 3.
# Combine all plots
figure3_regression_plot <- ggarrange(
  lma_ldmc_plot, ldmc_chl_plot, lma_chl_plot,
  leaf_area_lma_plot, leaf_area_ldmc_plot, leaf_area_chl_plot,
  ncol = 3, nrow = 2,  # 2 rows, 3 columns
  labels = c("a", "b", "c", "d", "e", "f"),  # Subplot labels
  common.legend = TRUE,  # Combine legends
  legend = "bottom"  # Legend at the bottom
)

# Print to console
print(figure3_regression_plot)

# Save the figure
ggsave(
  filename = "./figures/prelim/figure3_regression_plots.png",
  plot = figure3_regression_plot,
  width = 10, height = 7.5  # Adjust width and height for layout
)