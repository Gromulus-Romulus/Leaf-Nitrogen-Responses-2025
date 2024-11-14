# Mixed effects modeling
# @author: Nathan Malamud
# @date: 2024.11.12

library(tidyverse)
library(lme4)

# Load data for traits
# REMINDER: Set Working Directory -> Source File Location
data <- read_csv("./data/traits.csv")

# TODO: move prospect_rtm_output to data
# Load prospect measurements from spec curves
prospect <- read_csv("./data/molecular_content.csv") %>%
  select(-c(sampleID, species, treatment_mmol))
data <- merge(data, prospect, by = "barcodeID")

# Define Factor Levels (Treatment and Species)
data$treatment_mmol <- data$treatment_mmol
data$species <- as.factor(data$species)
levels(data$species) <- c("R. sativus", "B. officinalis", "H. vulgare")
data$LDMC <- as.numeric(data$LDMC)

data <- data %>%
  mutate(treatment_level = case_when(
    treatment_mmol <= 15 ~ "0 - 15 mmol",
    treatment_mmol <= 35 ~ "20 - 35 mmol",
  ))

# // -------------------------------------------------------------- //
# Hypothesis: growth response to nitrogen is explained by
# coordinated shifts between CHL content and leaf mass per area

# Lesson learned.... there's an important interaction effect between nitrogen adsorption and leaf investment
# p vals are greatly reduced by including N - a somewhat species-specific parameter
z1 <- lm(dry_whole_g ~ LMA + CHL + LDMC, data=data) # P val is quite low, R2 is reasonable, AIC is low

# This model fits ok
z2 <- lm(dry_whole_g ~ LMA * CHL + LDMC, data=data) # P val is quite low, R2 is reasonable, AIC is low:w

# The R2 values for these are better
# p vals are greatly reduced by including N - a somewhat species-specific parameter
z3 <- lmer(dry_whole_g ~ (1|species), data=data) # reduced
z4 <- lmer(dry_whole_g ~ LMA * CHL + (LMA|species) + (CHL|species), data=data) # full

# // -------------------------------------------------------------- //
library(lavaan)
library(semPlot)

sem_dat <- data %>% select(dry_whole_g, LMA, LDMC, CHL, N, treatment_mmol)

# Log 10 scale for SMA regression
sem_dat <- sem_dat %>%
  mutate(dry_whole_g = log10(dry_whole_g),
         LMA = log10(LMA),
         CHL = log10(CHL),
         LDMC = log10(LDMC))

names(sem_dat) = c("GRW", "LMA", "LDMC", "CHL", "MES", "N")

# Research questions:
# 1). To what extent do LMA, CHL, and LDMC explain nitrogen response (growth ~ N)
# 2). Do they explain variation to the same extent?

# // -------------------------------- //
# Here's a naive model
model1 <- '
  # Regression
  GRW ~ b1*LMA + b2*CHL + b3*LDMC
  LMA ~ N
  CHL ~ N
  LDMC ~ N

  # Variances
  LMA ~~ LMA
  CHL ~~ CHL
  CHL ~~ LMA
  LDMC ~~ LDMC
  LDMC ~~ LMA
'

# // -------------------------------- //
# // MODEL ESTIMATION
# No regression context - no intercept
model1.fit <- sem(model1, data=sem_dat, estimator='ML', meanstructure = F)

# // -------------------------------- //
# // MODEL EVALUATION + DIAGNOSTIC
summary(model1.fit,
        rsquare = T,
        standardized = T,
        fit.measures = T)

# TODO: model-implied covariance mismatch
resid(model1.fit)

# Plotting
semPaths(model1.fit,
         rotation = 2,
         layout="tree2",
         what = "std",
         posCol = "black",
         edge.width= 0.5,
         style = "Lisrel",
         fade = T,
         edge.label.position = 0.55)

# // -------------------------------- //
# Model specification
# Must be able to turn structure of hypothetical links
# into a specific grammar
model2 <- '
  # Latent variables
  ABS =~ N 
  ULR =~ b1*MES + b2*CHL
  
  GRW ~ ULR
  CHL ~ ABS
  MES ~ ABS + CHL

  # Variances
  CHL ~~ CHL
  MES ~~ MES
  CHL ~~ MES
'

# // -------------------------------- //
# // MODEL ESTIMATION
# No regression context - no intercept
model2.fit <- sem(model2, data=sem_dat, estimator='ML', meanstructure = F)

# // -------------------------------- //
# // MODEL EVALUATION + DIAGNOSTIC
summary(model2.fit,
        rsquare = T,
        standardized = T,
        fit.measures = T)

# TODO: model-implied covariance mismatch
resid(model2.fit)

# Plotting
semPaths(model2.fit,
         rotation = 2,
         layout="tree2",
         what = "std",
         posCol = "black",
         edge.width= 0.5,
         style = "Lisrel",
         fade = T,
         edge.label.position = 0.55)

# Test hypothesis b1 = b2 "Omega statistic"
#lavTestWald(model2.fit, constraints = "b1==b2")
  
