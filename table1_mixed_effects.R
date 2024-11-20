# Mixed effects modeling
# @author: Nathan Malamud
# @date: 2024.11.12

library(tidyverse)
library(lme4)
library(ggpubr)

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
data$species <- factor(data$species, levels=c("R. sativus", "B. officinalis", "H. vulgare"))
data$LDMC <- as.numeric(data$LDMC)

data <- data %>%
  mutate(treatment_level = case_when(
    treatment_mmol <= 5 ~ "0 - 5 mmol",
    treatment_mmol <= 15 ~ "10 - 15 mmol",
    treatment_mmol <= 25 ~ "20 - 25 mmol",
    treatment_mmol <= 35 ~ "30 - 35 mmol",
  ))

# Calculate metrics - unit leaf rate, leaf area ratio, leaf weight ratio, plant RGR
data$GRT <- (data$dry_whole_g / 6 * 7)

data <- data %>% filter(treatment_mmol > 0 & treatment_mmol < 35)

# Josef colors
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# facet wrap ULR with respect to nmmol across species
p1 <- ggplot(data, aes(x = log(CHL), y = log(GRT), color = species)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() + scale_color_manual(values = josef_colors)

# facet wrap LAR with respect to nmmol across species
p2 <- ggplot(data, aes(x = log(LMA), y = log(GRT), color = species)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() + scale_color_manual(values = josef_colors)

# facet wrap RGR with respect to nmmol across species
p3 <- ggplot(data, aes(x = log(CHL), y = log(LDMC), color = species)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() + scale_color_manual(values = josef_colors)

ggpubr::ggarrange(p1, p2, p3, ncol = 3, common.legend = T)

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

sem_dat <- data %>% select(GRT, LMA, LDMC, CHL, N, Phi_PS2, treatment_mmol, species)

# Log 10 scale for SMA regression
sem_dat <- sem_dat %>%
  mutate(GRT = log10(GRT),
         LMA = log10(LMA),
         CHL = log10(CHL),
         LDMC = log10(LDMC))

names(sem_dat) = c("GRT", "LMA", "LDM", "CHL", "MES", "N", "PS2", "SPC")

# Research questions:
# 1). To what extent do LMA, CHL, and LDMC explain nitrogen response (growth ~ N)
# 2). Do they explain variation to the same extent?

# // -------------------------------- //
# Here's a naive model
model1 <- '
  # Regression
  GRT ~ b1*LMA + b2*CHL + b3*LDM + b4*MES
  LMA ~ N
  CHL ~ N
  LDM ~ N
  MES ~ N

  # Co-variances
  LMA ~~ LMA
  CHL ~~ CHL
  CHL ~~ LMA
  CHL ~~ LDM
  LDM ~~ LDM
  LDM ~~ LMA
  MES ~~ CHL
  MES ~~ MES
  MES ~~ LMA
  MES ~~ LDM
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
  # Regression
  GRT ~ CHL + LDM + LMA
  LMA ~ N
  CHL ~ N
  LDM ~ N

  # Variances
  CHL ~~ CHL
  LDM ~~ LDM
  CHL ~~ LDM
  LMA ~~ LDM
  LMA ~~ CHL
  LMA ~~ LMA
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

# // -------------------------------- //
# Model specification
# Must be able to turn structure of hypothetical links
# into a specific grammar
model3 <- '
  # Regression
  GRT ~ LMA + LDM
  LMA ~ N
  LDM ~ N

  # Variances
  LMA ~~ LMA
  LDM ~~ LDM
  LMA ~~ LDM 
'

# // -------------------------------- //
# // MODEL ESTIMATION
# No regression context - no intercept
model3.fit <- sem(model3, data=sem_dat, estimator='ML', meanstructure = F)

# // -------------------------------- //
# // MODEL EVALUATION + DIAGNOSTIC
summary(model3.fit,
        rsquare = T,
        standardized = T,
        fit.measures = T)

# TODO: model-implied covariance mismatch
resid(model3.fit)

# Plotting
semPaths(model3.fit,
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

# // ------------------------------------//
# Introduction of hypothesized latent variables
# Leaf economics spectrum predicts tradeoff
# between structural and photosynthetic investment
model4 <- '
  # Latent variables,  structure and function
  STR =~ N
  FNC =~ CHL
  
  # Regression
  GRT ~ FNC + LMA + LDM
  CHL ~ N
  #MES ~ N + STR
  LMA ~ STR + LDM
  LDM ~ STR
  STR ~~ -1*FNC

  # Variances
  LMA ~~ LMA
  CHL ~~ CHL
  CHL ~~ LMA
  CHL ~~ LDM
  LDM ~~ LDM
  LDM ~~ LMA
  #MES ~~ CHL
  #MES ~~ MES
  #MES ~~ LMA
  #MES ~~ LDM 
'

# // -------------------------------- //
# // MODEL ESTIMATION
# No regression context - no intercept
model4.fit <- sem(model4, data=sem_dat, estimator='ML', meanstructure = F)

# // -------------------------------- //
# // MODEL EVALUATION + DIAGNOSTIC
summary(model4.fit,
        rsquare = T,
        standardized = T,
        fit.measures = T)

# TODO: model-implied covariance mismatch
resid(model4.fit)

# Plotting
semPaths(model4.fit,
         rotation = 2,
         layout="tree2",
         what = "std",
         posCol = "black",
         edge.width= 0.5,
         style = "Lisrel",
         fade = T,
         edge.label.position = 0.55)

# // ------------------------------------//

model5 <- '
  # Regression
  GRT ~ CHL + LMA + N
  LDM ~ CHL
  LMA ~ CHL + LDM
  CHL ~ N

  # Variances
  LMA ~~ LMA
  CHL ~~ CHL
  LDM ~~ LDM
  CHL ~~ LMA
  CHL ~~ LDM
  LDM ~~ LMA
'

# // -------------------------------- //
# // MODEL ESTIMATION
# No regression context - no intercept
model5.fit <- sem(model5, data=sem_dat, estimator='ML', meanstructure = F)

# // -------------------------------- //
# // MODEL EVALUATION + DIAGNOSTIC
summary(model5.fit,
        rsquare = T,
        standardized = T,
        fit.measures = T)

# TODO: model-implied covariance mismatch
resid(model5.fit)

# Plotting
semPaths(model5.fit,
         rotation = 2,
         layout="tree2",
         what = "std",
         posCol = "black",
         edge.width= 0.5,
         style = "Lisrel",
         fade = T,
         edge.label.position = 0.55)

# // ------------------------------------//
model6 <- '
  # Regression
  GRT ~ PS2 + LMA + N
  LDM ~ CHL
  LMA ~ CHL + LDM
  PS2 ~ CHL
  CHL ~ N

  # Variances
  LMA ~~ LMA
  CHL ~~ CHL
  LDM ~~ LDM
  CHL ~~ LMA
  CHL ~~ LDM
  LDM ~~ LMA
  PS2 ~~ PS2
  CHL ~~ PS2
'

# // -------------------------------- //
# // MODEL ESTIMATION
# No regression context - no intercept
model6.fit <- sem(model6, data=sem_dat, estimator='ML', meanstructure = F)

# // -------------------------------- //
# // MODEL EVALUATION + DIAGNOSTIC
summary(model6.fit,
        rsquare = T,
        standardized = T,
        fit.measures = T)

# TODO: model-implied covariance mismatch
resid(model6.fit)

# Plotting
semPaths(model6.fit,
         rotation = 2,
         layout="tree2",
         what = "std",
         posCol = "black",
         edge.width= 0.5,
         style = "Lisrel",
         fade = T,
         edge.label.position = 0.55)

# // ------------------------------------//

library(piecewiseSEM)

# Select the variables and transform them
sem_dat <- data %>% 
  select(GRT, LMA, LDMC, CHL, Phi_PS2, treatment_mmol) %>%
  mutate(
    GRT = log10(GRT),
    LMA = log10(LMA),
    CHL = log10(CHL),
    LDMC = log10(LDMC)
  )

# Rename columns
names(sem_dat) <- c("GRT", "LMA", "LDM", "CHL", "PS2", "N")

# Specify the psem model
modelList <- psem(
  lm(GRT ~ LMA + PS2 + N, data=sem_dat),
  lm(LMA ~ CHL + LDM, data=sem_dat),
  lm(PS2 ~ CHL, data=sem_dat),
  lm(LDM ~ CHL, data=sem_dat),
  lm(CHL ~ N, data=sem_dat),
  sem_dat
)

# Evaluate model fit
summary(modelList)

# Address conflict using conserve = T
summary(modelList, conserve = T)

# lavaan fitting and plotting
lavaan_model <- '
  # Regression equations
  GRT ~ LMA + PS2 + N
  LMA ~ CHL + LDM
  PS2 ~ CHL
  LDM ~ CHL
  CHL ~ N

  # Covariance structure
  LMA ~~ LDM
  LMA ~~ CHL
  LDM ~~ CHL
'

# Fit the model using lavaan
library(lavaan)
sem_fit <- sem(lavaan_model, data = sem_dat)

# Summarize the model
summary(sem_fit, fit.measures = TRUE)

# Check residuals
resid(sem_fit)

# Plot as a graph
semPaths(sem_fit,
         rotation = 2,
         layout="tree2",
         what = "std",
         posCol = "black",
         edge.width= 0.5,
         style = "Lisrel",
         fade = T,
         edge.label.position = 0.55)
