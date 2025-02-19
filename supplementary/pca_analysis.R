##' Objective: Plot differences in spectral signatures
##' across low and high treatments
##'
##' Adapted from code written for 7th
##' annual plant functional traits course.
##' 
##' TODO: Do a wavelength-by-wavelength T test
##' and add p values to each chart
##'   Reference: https://journals.ashs.org/hortsci/view/journals/hortsci/53/5/article-p669.xml
##' 
##' @author [Nicole Bison, Nathan Malamud]
##' @date [2024.10.02]
##' 
##' Source:
##'   https://github.com/MichaletzLab/pftc7_spectroscopy
##'   

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(factoextra)

# Load data for traits
data <- read_csv("./data/traits.csv")

# Define a custom theme for all plots
custom_theme <- theme_classic() +
  theme(
    text = element_text(family = "sans", size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 14, hjust = 0.5),
    aspect.ratio = 1
  )

# Define custom colors for species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Select relevant columns and preprocess the data
pca_data <- data %>%
  select(species, CHL, LMA, LDMC, Phi_PS2, dry_whole_g, treatment_mmol) %>%
  drop_na() %>%
  mutate(
    GRT = dry_whole_g / 6 * 7,  # Calculate growth rate (example metric)
    treatment_level = case_when(
      treatment_mmol <= 15 ~ "0 - 15 mM",
      treatment_mmol <= 35 ~ "20 - 35 mM",
      TRUE ~ "Other"
    )
  ) %>%
  filter(treatment_level != "Other")  # Exclude "Other" bin if needed

# Filter out treatment_mmol
pca_data <- pca_data %>% select(-c(treatment_mmol, dry_whole_g))

# List to store plots
plots <- list()

# Loop through treatment bins to create PCA plots
for (bin in unique(pca_data$treatment_level)) {
  bin_data <- pca_data %>% filter(treatment_level == bin)
  
  # Fit PCA model
  pca <- prcomp(bin_data %>% select(-species, -treatment_level), center = TRUE, scale. = TRUE)
  
  # Generate PCA plot for each bin
  pca_plot <- fviz_pca_biplot(
    pca,
    geom.ind = "point",
    col.ind = bin_data$species,
    palette = josef_colors,
    addEllipses = TRUE,  # Ellipses are only calculated if enough points are available
    legend.title = "Species",
    repel = TRUE
  ) +
    custom_theme +
    ggtitle(paste("PCA for", bin))
  
  # Add plot to list
  plots[[bin]] <- pca_plot
}

# Arrange plots in a grid (2 plots)
arranged_plots <- ggarrange(
  plots[["0 - 15 mM"]],
  plots[["20 - 35 mM"]],
  ncol = 2, nrow = 1, labels=c("a", "b"),
  common.legend = TRUE, legend = "bottom"
)

# Display the arranged plots in Rstudio
print(arranged_plots)

# Save the combined plot as a PDF
ggsave(filename = "./figures/pca_plot_treatments.pdf", 
       plot = arranged_plots, device = "pdf", width = 10, height = 5)
