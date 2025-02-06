##' Fit piecewiseSEM framework.
##' 14 paths plus variables ~ 7 observations per path.
##' 
##' @author [Nathan D Malamud, Sean T. Michaletz]
##' @date [2025-01-31]
##' 

# Libraries ----
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(scales)
library(smatr)
library(piecewiseSEM)

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

# Specify model structure ----
# Define paths based on covariance vs directed
# structure. Recall path diagram with dashed and solid lines.
# - corClasses : goes into model specification syntax
# TODO: look up syntax
model <- psem(
  lme(y1 ~ x1 + x2, random = ~1 | group, data = mydata),
  lme(y2 ~ x3, random = ~1 | group, data = mydata),
  corClasses = corAR1(form = ~1 | group
)  # Adds correlation structure) summary(model)

# Fit models iteratively -----
# TODO

## 1) First pass - Naive fit ----
# TODO

## 2) Second pass - one path at a time ----
# TODO

## 3) Third pass - one path at a time ----
# TODO

# ...........

# Final model ----
# TODO

# Shipley reccommends d-separation - explores "missing" paths
# TODO

# Visualize path diagram ----
# TODO
