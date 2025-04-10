##' Run PCA analysis of traits.
##' View differences across crops and nitrogen treatment.
##' 
##' @author [Nathan D. Malamud]
##' @date [2025-01-27] 
##'

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

# Define traits of interest here
traits_of_interest <- c("LDMC", "LMA", "CHL", "GRT")

# Filter by metrics of interest only
traits <- traits %>%
  select(barcodeID, species, treatment_mmol,
         all_of(traits_of_interest))

# Styling and aesthetics ----
# Define a custom theme for all plots
custom_theme <- theme_classic() +  # Base theme
  theme(
    # Text and font styling
    text = element_text(family = "sans", size = 12),
    axis.text = element_text(size = 10),  
    legend.text = element_text(size = 11),
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
  "area_cm2" = "Leaf Area (cm²)",
  "treatment_mmol" = "N (mM)"
)

# PCA Analysis ----
# Select relevant variables and scale data
pca_data <- traits %>% select(species, treatment_mmol, LMA, LDMC, CHL, GRT)
pca_data$treatment_mmol <- as.factor(pca_data$treatment_mmol)

pca <- prcomp(pca_data %>% select(-c(species, treatment_mmol)), center = TRUE, scale. = TRUE)

# Calculate % variance explained by each PC
explained_variation <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1) 

## Ordination biplot, groups by species ----
pca_plot_species <- fviz_pca_biplot(
  pca,
  geom.ind = "point",  # Plot individuals as points
  col.ind = pca_data$species,  # Color individuals by species
  palette = josef_colors,  # Custom species colors
  addEllipses = FALSE,  # Add confidence ellipses for each group
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
    legend.position = "bottom",
    aspect.ratio = 1
  )

# Print the PCA plot
print(pca_plot_species)

## Ordination groups by treatment ----
# Define custom colors for specific nitrogen levels
custom_colors1 <- c(
  "0"  = "#D2B48C",  # Darker beige
  "5"  = "#E6A96B",  # Light brownish-orange
  "10" = "#FDBE85",  # Soft orange
  "15" = "#FDAE61",  # Bright orange
  "20" = "#F46D43",  # Deep orange
  "25" = "#E31A1C",  # Red
  "30" = "#BD0026",  # Dark red
  "35" = "#800026"   # Deep maroon
)

# Highlight key treatments using custom colors for specific nitrogen levels
custom_colors2 <- c(
  "0"  = "#FDBE85",  # Darker beige
  "5"  = "grey",  # Light brownish-orange
  "10" = "grey",  # Soft orange
  "15" = "#FDAE61",  # Bright orange
  "20" = "grey",  # Deep orange
  "25" = "grey",  # Red
  "30" = "#BD0026",  # Dark red
  "35" = "grey"   # Deep maroon
)


# Define distinct shapes (must match the number of nitrogen levels)
custom_shapes <- c(16, 17, 15, 3, 8, 10, 7, 4)  # 8 values

pca_plot_treatment_ellipse <- fviz_pca_biplot(
  pca,
  geom.ind = "point",  # Plot individuals as points
  col.ind = factor(pca_data$treatment_mmol),  # Ensure discrete coloring
  shape.ind = factor(pca_data$treatment_mmol),  # Ensure discrete shapes
  addEllipses = TRUE,  # Add confidence ellipses for each group
  ellipse.level = 0.95,  # Set confidence level for ellipses
  ellipse.alpha = 0.0,  # Set fill transparency to 0 for ellipses
  legend.title = "Nitrogen addition (mM)",  # Legend title
  repel = TRUE  # Avoid overlap of labels
) +
  labs(
    title = NULL,
    x = paste0("PC1 (", explained_variation[1], "% variance)"),
    y = paste0("PC2 (", explained_variation[2], "% variance)")
  ) +
  custom_theme +
  scale_color_manual(values = custom_colors2) +  # Assign custom colors
  scale_shape_manual(values = custom_shapes) +  # Assign distinct shapes
  theme(
    legend.position = "bottom",
    aspect.ratio = 1
  )

pca_plot_treatment_no_ellipse <- fviz_pca_biplot(
  pca,
  geom.ind = "point",  # Plot individuals as points
  col.ind = factor(pca_data$treatment_mmol),  # Ensure discrete coloring
  shape.ind = factor(pca_data$treatment_mmol),  # Ensure discrete shapes
  addEllipses = FALSE,  # No confidence ellipses for each group
  ellipse.level = 0.95,  # Set confidence level for ellipses
  ellipse.alpha = 0.0,  # Set fill transparency to 0 for ellipses
  legend.title = "Nitrogen addition (mM)",  # Legend title
  repel = TRUE  # Avoid overlap of labels
) +
  labs(
    title = NULL,
    x = paste0("PC1 (", explained_variation[1], "% variance)"),
    y = paste0("PC2 (", explained_variation[2], "% variance)")
  ) +
  custom_theme +
  scale_color_manual(values = custom_colors1) +  # Assign custom colors
  scale_shape_manual(values = custom_shapes) +  # Assign distinct shapes
  theme(
    legend.position = "bottom",
    aspect.ratio = 1
  )

# Full plot with ellipses
pca_plot_treatment <- ggarrange(
  pca_plot_treatment_ellipse,
  pca_plot_treatment_no_ellipse,
  nrow = 1, labels = c("a", "b")
)

# Print the PCA plot
print(pca_plot_treatment)

# Create a scree plot for eigenvalues
pca_eigenplot <- fviz_eig(pca, addlabels = TRUE, geom = c("line", "point")) +
  custom_theme +
  geom_segment(
    aes(x = seq_along(pca$sdev), xend = seq_along(pca$sdev), 
        y = 0, yend = pca$sdev^2 / sum(pca$sdev^2) * 100),
    linetype = "dashed", color = "gray"
  ) +
  geom_point(
    aes(x = seq_along(pca$sdev), 
        y = pca$sdev^2 / sum(pca$sdev^2) * 100), color = "black"
  ) +
  theme(
    aspect.ratio = 1
  ) +
  labs(
    y = "% Variance",
    title = NULL
  )

# Full plot with ellipses
pca_plot_treatment <- ggarrange(
  pca_plot_treatment_ellipse,
  pca_plot_treatment_no_ellipse,
  nrow = 1, labels = c("", "")
)

print(pca_plot_treatment)

print(pca_eigenplot)

# Save the PCA plot as Figure 4
# TODO: investigated ignored "override.aes" warnings
ggsave(
  filename = "./figures/treatment_pca_plot.png",
  plot = pca_plot_treatment,
  width = 10, height = 5,
  bg = "white"
)

# Save as svg
library(svglite)
ggsave("treatment_pca_plot.svg", pca_plot_treatment,
       width = 16, height = 6, bg="white")

# Save scree plot
ggsave(
  filename = "./figures/pca_eigenplot.png",
  plot = pca_eigenplot,
  width = 5, height = 5,
  bg = "white"
)
