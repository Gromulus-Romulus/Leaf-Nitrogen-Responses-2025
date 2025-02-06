## Visualize Final SEM Model
## Author: Nathan D. Malamud
## Date: February 3rd, 2025

# Load necessary libraries ----
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# Modeling prototype (hypothesized structure) ----

library(DiagrammeR)

grViz("
  digraph SEM {
    # Global node style
    node [shape = box, style = filled, fillcolor = white, fontname = 'Helvetica', fontsize = 12, width = 1.5, height = 0.6]

    # Define nodes
    Nitrogen [label = 'Nitrogen']
    LMA [label = 'LMA']
    LDMC [label = 'LDMC']
    CHL [label = 'Chlorophyll']
    GrowthRate [label = 'Growth Rate']

    # Enforce vertical stacking of traits
    { rank = same; LMA; LDMC; CHL; }

    # Enforce position relationships
    Nitrogen -> LMA [penwidth = 1, color = black]
    Nitrogen -> LDMC [penwidth = 1.5, color = black]
    Nitrogen -> CHL [penwidth = 2, color = black]
    LMA -> GrowthRate [penwidth = 1, color = black]
    LDMC -> GrowthRate [penwidth = 3, color = black]
    CHL -> GrowthRate [penwidth = 2.5, color = black]

    # Correlation paths
    LMA -> LDMC [dir = both, style = dashed, color = red, penwidth = 1.5]
    LMA -> CHL [dir = both, style = dashed, color = red, penwidth = 1.5]
    LDMC -> CHL [dir = both, style = dashed, color = red, penwidth = 1.5]

    # Optional dashed path
    Nitrogen -> GrowthRate [style = dashed, color = black, constraint = false, penwidth = 1.5]

    # Graph attributes
    graph [rankdir = TB, splines = true, overlap = false] # Use TB (top-bottom) layout
  }
")


DiagrammeRsvg::export_svg(grViz("
  digraph SEM {
    ...
  }
")) %>%
  charToRaw() %>%
  rsvg::rsvg_svg("Tight_SEM_Diagram.svg")

# Import model from SEM_analysis_v2 (by Sean M) ----
# Note: unlike traditional SEM fits (e.g. in lavaan),
# the piecewiseSEM library fits hypothesize causal relationships
# with mixed effects.
#
# We can not use the conventional semPlot package for
# visualiztion of the DAG (Direct Acylclic Graph)
# However, there exist alternative visualization libraries.
# mod <- sem1.3

# Summary of final SEM model ----
sem_summary <- summary(mod)

# Print the coefficients table for reference
print(sem_summary$coefficients)

# Create a summary data frame based on the SEM model coefficients
sem_summary_table <- data.frame(
  Path = c(
    "Nitrogen -> LMA",
    "Nitrogen -> LDMC",
    "Nitrogen -> CHL",
    "LMA -> GrowthRate",
    "LDMC -> GrowthRate",
    "CHL -> GrowthRate",
    "Nitrogen -> GrowthRate",
    "LMA ~~ LDMC" # Correlation path
  ),
  Estimate = c(-0.0023, -0.0020, 0.0061, 2.4784, -2.4999, 1.3530, 0.0119, 0.7344),
  Std.Error = c(0.0007, 0.0006, 0.0011, 1.1769, 0.6122, 0.2159, 0.0028, NA),
  DF = c(93.1650, 93.0059, 93.1536, 2.5084, 90.0092, 88.6182, 88.6549, 97.0000),
  Crit.Value = c(-3.4738, -3.3602, 5.4763, 2.1059, -4.0834, 6.2670, 4.2218, 10.4904),
  P.Value = c(0.0008, 0.0011, 0.0000, 0.1434, 0.0001, 0.0000, 0.0001, 0.0000),
  Std.Estimate = c(-0.3074, -0.2162, 0.4865, 0.4493, -0.5599, 0.4119, 0.2897, 0.7344),
  Significant = c("***", "**", "***", "", "***", "***", "***", "***") # Significance stars
)

# Print the summary table for review
print(sem_summary_table)

# Save the summary table as a CSV file (optional)
write.csv(sem_summary_table, "SEM_Summary_Table.csv", row.names = FALSE)

# Load necessary libraries
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# Generate SEM Path Diagram
sem_diagram <- grViz("
  digraph SEM {
    # Global node style with Helvetica font
    node [shape = box, style = filled, fillcolor = white, fontname = \"Helvetica\", fontsize = 12]

    # Define nodes with consistent terminology
    Nitrogen [label = \"Soil Nitrogen Content\"]
    LMA [label = \"LMA\"]
    LDMC [label = \"LDMC\"]
    CHL [label = \"Chlorophyll Content\"]
    GrowthRate [label = \"Growth Rate\"]

    # Define causal paths with widths scaled by effect size
    Nitrogen -> LMA [penwidth = 1, color = black]      // Weak effect
    Nitrogen -> LDMC [penwidth = 1.5, color = black]   // Moderate effect
    Nitrogen -> CHL [penwidth = 2, color = black]      // Strong effect
    Nitrogen -> GrowthRate [penwidth = 1.2, color = black] // Moderate effect
    LDMC -> GrowthRate [penwidth = 3, color = black]   // Strongest effect
    CHL -> GrowthRate [penwidth = 2.5, color = black]  // Strong effect

    # Optional dashed lines for non-significant paths
    LMA -> GrowthRate [style = dashed, penwidth = 1, color = gray] // Non-significant path

    # Define correlations (bidirectional paths)
    LMA -> LDMC [dir=both, penwidth = 1.8, color = grey]
  }
")

# Export SEM Diagram as SVG for further editing in Inkscape
svg_output <- export_svg(sem_diagram) %>%
  charToRaw() %>%
  rsvg_svg("SEM_Diagram.svg")

# Print confirmation message
print("SEM Path Diagram saved as SEM_Diagram.svg. You can now edit it in Inkscape.")
