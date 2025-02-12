## Visualize Final SEM Model 1.1 (Automated Version)
## Author: Nathan D. Malamud
## Date: February 3rd, 2025

# Load necessary libraries ----
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(dplyr)

# Load the final SEM model ----
mod1.1 <- sem1.1  # Ensure this is defined earlier in your workflow
mod <- sem1.1

# Extract coefficients from the model ----
sem_summary <- summary(mod)

# Save coefficients for later use ----
coef_values <- list(
  "log10_LMA ~ treatment_mmol" = -0.0023,
  "log10_LDMC ~ treatment_mmol" = -0.0020,
  "log10_CHL ~ treatment_mmol" = 0.0061,
  "log10_grt_g_d ~ log10_LMA" = 2.4784,
  "log10_grt_g_d ~ log10_LDMC" = -2.4999,
  "log10_grt_g_d ~ log10_CHL" = 1.3530,
  "log10_grt_g_d ~ treatment_mmol" = 0.0119,
  "log10_LMA ~~ log10_LDMC" = 0.7344,
  "log10_LMA ~~ log10_CHL" = -0.0069,
  "log10_LDMC ~~ log10_CHL" = 0.0412
)

# Function to determine line properties ----
get_edge_style <- function(estimate) {
  color <- ifelse(estimate < 0, "red", "black")
  style <- ifelse(abs(estimate) < 0.1, "dashed", "solid")
  penwidth <- max(1, abs(estimate) * 2)
  return(paste0("color=", color, ", style=", style, ", penwidth=", penwidth))
}

# Print coefficients to the R console ----
print(sem_summary$coefficients)

# Define node positions ----
nodes <- paste0(
  '"treatment_mmol" [label="Soil Nitrogen Content"];
  "log10_LMA" [label="LMA", pos="1 2!"];
  "log10_LDMC" [label="LDMC", pos="2 2!"];
  "log10_CHL" [label="Chlorophyll Content", pos="3 2!"];
  "log10_grt_g_d" [label="Growth Rate"];' 
)

# Define edges ----
edges <- paste0(
  '"treatment_mmol" -> "log10_LMA" [label="', coef_values[["log10_LMA ~ treatment_mmol"]], '", ', get_edge_style(coef_values[["log10_LMA ~ treatment_mmol"]]), '];
  "treatment_mmol" -> "log10_LDMC" [label="', coef_values[["log10_LDMC ~ treatment_mmol"]], '", ', get_edge_style(coef_values[["log10_LDMC ~ treatment_mmol"]]), '];
  "treatment_mmol" -> "log10_CHL" [label="', coef_values[["log10_CHL ~ treatment_mmol"]], '", ', get_edge_style(coef_values[["log10_CHL ~ treatment_mmol"]]), '];
  "log10_LMA" -> "log10_grt_g_d" [label="', coef_values[["log10_grt_g_d ~ log10_LMA"]], '", ', get_edge_style(coef_values[["log10_grt_g_d ~ log10_LMA"]]), '];
  "log10_LDMC" -> "log10_grt_g_d" [label="', coef_values[["log10_grt_g_d ~ log10_LDMC"]], '", ', get_edge_style(coef_values[["log10_grt_g_d ~ log10_LDMC"]]), '];
  "log10_CHL" -> "log10_grt_g_d" [label="', coef_values[["log10_grt_g_d ~ log10_CHL"]], '", ', get_edge_style(coef_values[["log10_grt_g_d ~ log10_CHL"]]), '];
  "treatment_mmol" -> "log10_grt_g_d" [label="', coef_values[["log10_grt_g_d ~ treatment_mmol"]], '", ', get_edge_style(coef_values[["log10_grt_g_d ~ treatment_mmol"]]), '];
  "log10_LMA" -> "log10_LDMC" [dir=both, label="', coef_values[["log10_LMA ~~ log10_LDMC"]], '", color=gray, penwidth=1];
  "log10_LMA" -> "log10_CHL" [dir=both, label="', coef_values[["log10_LMA ~~ log10_CHL"]], '", color=gray, penwidth=1, style=dashed];
  "log10_LDMC" -> "log10_CHL" [dir=both, label="', coef_values[["log10_LDMC ~~ log10_CHL"]], '", color=gray, penwidth=1, style=dashed];'
)

# Manually construct path diagram using DiagrammeR ----
graph_code <- paste0(
  "digraph SEM {\n",
  "  rankdir=TB;\n",
  "  node [shape=box, style=filled, fillcolor=white, fontname=\"Helvetica\", fontsize=12];\n",
  nodes, "\n",
  edges, "\n",
  "}"
)

# Render the SEM diagram ----
grViz(graph_code)

# Export SEM Diagram as SVG and PNG ----
# Export SEM Diagram as SVG and PNG ----
svg_filename <- paste0("SEM_Model_", "sem1.1", ".svg")

svg_output <- export_svg(sem_diagram)
writeLines(svg_output, svg_filename)
rsvg_png(charToRaw(svg_output), png_filename)