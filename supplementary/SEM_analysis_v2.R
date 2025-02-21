# SEM analysis
# Sean Michaletz (sean.michaletz@gmail.com), 2 Feb 2025
##' @author [Sean T. Michaletz, Nathan D. Malamud]
##' @title [Path diagram for the best-fit SEM model`]

# 0. Initialize ----

# Load data
df <- read.csv("./data/traits.csv",header=T) # Data from Nate

# Calculate aboveground growth rate
df$grt_g_d <- df$dry_whole_g / 6*7 #days in experiment

# Here, we use untransformed nitrogen, log(traits), and log(growth rate). This
# is based entirely on theory and precedent in the literature. For example, 
# log(traits) ~ nitrogen was used in Halbritter et al. in review and many others,
# log(trait) ~ log(trait) is based on leaf economics spectrum theory, and
# log(growth rate) ~ log(traits) is based on allometry and metabolic theory. 
# Model selection is used to identify AIC preference for random effects by 
# species.

# This approach doesn't always provide the best description of the relationships
# (particularly for saturation of CHL and growth rates across nitrogen), but:
# 1) It enables retention of the 0 mmol treatment, which needs to be removed for
#    most other transformations such as 1/y ~ 1/x;
# 2) It yields a simple, parsimonious SEM that is consistent with existing trait-
#    based scaling theory;
# 3) It still captures the general, zeroth-order relationships;
# 3) We must remember that the goal of SEM is to evaluate hypothesized causal 
#    relationships, not to predict y from x. Thus, SEMs use linear models that 
#    provide effect size and significance even though they may not be the best 
#    description of the data (in terms of model form).

# If we ultimately use this approach, the above text should be included in the 
# Methods section.


# 1. Trait-nitrogen relationships ----
### Model selection ----
# Compete models for log(trait) ~ treatment_mmol using random intercept, random
# slope, and random slope+intercept by species.

# Load required packages
library(lmerTest)    # for lmer() with p-values
library(MuMIn)       # for r.squaredGLMM()
library(dplyr)

# Define the traits and predictor (only treatment_mmol)
traits <- c("LMA", "LDMC", "CHL")
predictors <- c("treatment_mmol")
predictor_labels <- c("treatment_mmol")
# Define random-effects types.
random_effects_types <- c("intercept", "slope", "intercept+slope")

# Initialize a list to store results.
results_list <- list()

# Loop over each trait, the predictor, and each random-effects structure.
for (trait in traits) {
  for (i in seq_along(predictors)) {
    pred <- predictors[i]
    pred_label <- predictor_labels[i]
    
    for (re_type in random_effects_types) {
      # Define the random-effects part of the formula based on the type.
      if (re_type == "intercept") {
        re_str <- "(1 | species)"
      } else if (re_type == "slope") {
        # Random slope only (no random intercept)
        re_str <- paste0("(0 + ", pred, " | species)")
      } else if (re_type == "intercept+slope") {
        re_str <- paste0("(", pred, " | species)")
      }
      
      # Build the full model formula.
      # The fixed-effects part is log(trait) ~ treatment_mmol.
      formula_str <- paste0("log10(", trait, ") ~ ", pred, " + ", re_str)
      formula_obj <- as.formula(formula_str)
      
      # Fit the model. Use tryCatch() to avoid halting on error.
      model_fit <- tryCatch(
        lmer(formula_obj, data = df),
        error = function(e) NULL
      )
      
      # If the model failed to fit, record NAs and continue.
      if (is.null(model_fit)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          trait = trait,
          predictor = pred_label,
          random_effects = re_type,
          singular = NA,
          AIC = NA,
          p_value = NA,
          r2_marginal = NA,
          r2_conditional = NA,
          stringsAsFactors = FALSE
        )
        next
      }
      
      # Check whether the fit was singular.
      sing <- isSingular(model_fit, tol = 1e-4)
      
      # Get the AIC.
      aic_val <- AIC(model_fit)
      
      # Extract the p-value for the fixed effect treatment_mmol.
      summ <- summary(model_fit)
      p_val <- summ$coefficients[2, "Pr(>|t|)"]
      
      # Compute marginal and conditional R2 using MuMIn.
      r2_vals <- r.squaredGLMM(model_fit)
      r2_marg <- r2_vals[1]
      r2_cond <- r2_vals[2]
      
      # Store the results.
      results_list[[length(results_list) + 1]] <- data.frame(
        trait = trait,
        predictor = pred_label,
        random_effects = re_type,
        singular = sing,
        AIC = aic_val,
        p_value = p_val,
        r2_marginal = r2_marg,
        r2_conditional = r2_cond,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine all results into one data frame.
results_df <- do.call(rbind, results_list)
# Subset only non-singular fits
results_df <- subset(results_df, singular == FALSE)
print(results_df)

# AIC prefers random intercept for LDMC, LMA, and grt_g_d, and cannot 
# distinguish between models for CHL. So, use random intercepts for all in SEMs.


### Plots ----
# Make plots, fitting regressions with species as random intercept
library(ggplot2)
library(lme4)
library(dplyr)
library(rlang)  # for tidy evaluation with sym()

# Define a function to plot a given trait against treatment_mmol
plot_trait <- function(trait) {
  # Build the model formula with a random intercept for species:
  # log(trait) ~ treatment_mmol + (1 | species)
  model_formula <- as.formula(paste0("log10(", trait, ") ~ treatment_mmol + (1 | species)"))
  
  # Fit the mixed model
  model <- lmer(model_formula, data = df)
  
  # Create new data for prediction: a grid of treatment_mmol values for each species
  pred_data <- expand.grid(
    treatment_mmol = seq(min(df$treatment_mmol), max(df$treatment_mmol), length.out = 100),
    species = unique(df$species)
  )
  
  # Get predictions from the model.
  # The predictions are on the log scale and will be back-transformed.
  pred_data$predicted_log_trait <- predict(model, newdata = pred_data, re.form = NULL)
  
  # Create the ggplot:
  # - geom_point() plots the raw data points (with y on original scale).
  # - geom_line() overlays the predicted values (back-transformed using exp()).
  # - scale_y_log10() sets the y-axis to a logarithmic scale.
  p <- ggplot() +
    geom_point(data = df, aes(x = treatment_mmol, y = !!sym(trait), color = species)) +
    geom_line(data = pred_data, aes(x = treatment_mmol, y = 10^(predicted_log_trait), color = species)) +
    scale_y_log10() +
    labs(x = "Nitrogen treatment (mmol)",
         y = trait,
         title = paste(trait, "vs. Nitrogen Treatment"),
         color = "Species") +
    theme_minimal()
  
  return(p)
}

# Generate plots for each of the three traits.
p_LMA  <- plot_trait("LMA")
p_LDMC <- plot_trait("LDMC")
p_CHL  <- plot_trait("CHL")
p_grt_g_d  <- plot_trait("grt_g_d")

# Display the plots (e.g., in an interactive session)
print(p_LMA)
print(p_LDMC)
print(p_CHL)
print(p_grt_g_d)

# Chlorophyll and growth rate appear to saturate across the nitrogen gradient,
# which the transformation and models don't capture well. However, the models do
# capture the general increase in chlorophyll and growth rates with nitrogen, 
# so in that sense they are performing adequately. We should remember the goal 
# here isn't to predict y from x, but rather to test hypothesized causal 
# relationships, so these models can adequately do that.


# 2. Growth-nitrogen relationships ----

### Model selection ----
library(lme4)        # For lmer and singular fit checking
library(MuMIn)       # For r.squaredGLMM()
library(performance) # For alternative R^2 calculation

# Define candidate models as functions that accept a REML argument.
models <- list(
  RI = function(reml) {
    lmer(log10(grt_g_d) ~ treatment_mmol + (1 | species),
         data = df, REML = reml)
  },
  RS = function(reml) {
    lmer(log10(grt_g_d) ~ treatment_mmol + (0 + treatment_mmol | species),
         data = df, REML = reml)
  },
  RIS = function(reml) {
    lmer(log10(grt_g_d) ~ treatment_mmol + (1 + treatment_mmol | species),
         data = df, REML = reml)
  }
)

# Prepare an empty results data frame
results_df <- data.frame(
  Model          = character(),
  AIC            = numeric(),
  p_value        = numeric(),
  R2_marginal    = numeric(),
  R2_conditional = numeric(),
  Singular_Fit   = logical(),
  stringsAsFactors = FALSE
)

# Loop over each candidate model
for(model_name in names(models)){
  
  # This will track if any step emits a singular-fit warning
  singular_warning <- FALSE
  
  # 1) Fit model with ML (REML = FALSE) for AIC and p-values
  fit_ml <- withCallingHandlers(
    models[[model_name]](reml = FALSE),
    warning = function(w) {
      if(grepl("singular", w$message, ignore.case = TRUE)) {
        singular_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  
  # Extract AIC from the ML fit
  model_aic <- AIC(fit_ml)
  
  # Extract p-value for treatment_mmol from the fixed-effects summary
  model_summary <- summary(fit_ml)
  p_val <- NA
  if("treatment_mmol" %in% rownames(model_summary$coefficients)){
    p_val <- model_summary$coefficients["treatment_mmol", "Pr(>|t|)"]
  }
  
  # 2) Refit the model with REML = TRUE to compute R^2 values
  fit_reml <- withCallingHandlers(
    update(fit_ml, REML = TRUE),
    warning = function(w) {
      if(grepl("singular", w$message, ignore.case = TRUE)) {
        singular_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  
  # 3) Attempt R^2 with MuMIn::r.squaredGLMM
  r2_vals <- withCallingHandlers(
    r.squaredGLMM(fit_reml),
    warning = function(w) {
      if(grepl("singular", w$message, ignore.case = TRUE)) {
        singular_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  
  r2_marg <- r2_vals["R2m"]
  r2_cond <- r2_vals["R2c"]
  
  # If either R^2 value is NA, attempt with performance::r2
  if(is.na(r2_marg) || is.na(r2_cond)){
    r2_perf <- withCallingHandlers(
      r2(fit_reml),
      warning = function(w) {
        if(grepl("singular", w$message, ignore.case = TRUE)) {
          singular_warning <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    )
    r2_marg <- r2_perf$R2_marginal
    r2_cond <- r2_perf$R2_conditional
  }
  
  # Append results to the data frame
  results_df <- rbind(
    results_df,
    data.frame(
      Model            = model_name,
      AIC              = model_aic,
      p_value          = p_val,
      R2_marginal      = r2_marg,
      R2_conditional   = r2_cond,
      Singular_Fit     = singular_warning,
      stringsAsFactors = FALSE
    )
  )
}

# Subset only non-singular fits
results_df <- subset(results_df, Singular_Fit == FALSE)
# Print the summary data frame
print(results_df)
# AIC prefers random intercept


### Plot ----
p_grt <- ggplot(df, aes(x = treatment_mmol, y = grt_g_d, color = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_y_log10() +
  labs(x = "Nitrogen treatment (mmol)",
       y = "Growth rate (g/day)",
       title = "Growth rate vs. Nitrogen Treatment",
       color = "Species") +
  theme_bw()
print(p_grt)
# As for chlorophyll, we're not really capturing the saturation in growth rates
# across nitrogen. However, we are capturing a general increase in growth rates
# with nitrogen.


# 3. Trait-trait relationships ----
### SMA regressions ----
library(smatr)
sma1 <- sma(LDMC ~ LMA * species, data = df, log = "xy")
summary(sma1) # all significant
sma2 <- sma(CHL ~ LMA * species, data = df, log = "xy")
summary(sma2) # B. officinalis and H. vulgare not significant
sma3 <- sma(CHL ~ LDMC * species, data = df, log = "xy")
summary(sma3) # B. officinalis and H. vulgare not significant (R. sativus is marginally significant)

### Plots ----
# Make plots, including regression lines only when significant
library(ggplot2)
library(ggpmisc)   # Provides stat_ma_line()
library(dplyr)
library(gridExtra)
library(scales)    # For hue_pal()

# A function that creates a pairwise plot with points and, conditionally, SMA regression lines by species.
plot_pair <- function(xvar, yvar) {
  p <- ggplot(df, aes_string(x = xvar, y = yvar, color = "species")) +
    geom_point() +
    # Plot raw data with log10-transformed axes.
    scale_x_log10() +
    scale_y_log10() +
    labs(x = xvar, y = yvar) +
    theme_bw()
  
  # Add regression lines conditionally based on xvar and yvar.
  if(xvar == "LMA" && yvar == "CHL") {
    # For CHL ~ LMA, add regression lines only for species other than B. officinalis and H. vulgare.
    p <- p + stat_ma_line(
      data = subset(df, !(species %in% c("B. officinalis", "H. vulgare"))),
      mapping = aes_string(x = xvar, y = yvar, group = "species"),
      method = "SMA", se = FALSE, linewidth = 1
    )
  } else if(xvar == "LDMC" && yvar == "CHL") {
    # For CHL ~ LDMC, add regression lines only for species other than B. officinalis and H. vulgare.
    p <- p + stat_ma_line(
      data = subset(df, !(species %in% c("B. officinalis", "H. vulgare"))),
      mapping = aes_string(x = xvar, y = yvar, group = "species"),
      method = "SMA", se = FALSE, linewidth = 1
    )
  } else {
    # For other pairs (e.g., LMA vs. LDMC), add regression lines for all species.
    p <- p + stat_ma_line(aes_string(group = "species"), method = "SMA", se = FALSE, linewidth = 1)
  }
  
  return(p)
}

# Create the three pairwise plots:
# TODO: arbitrary plot pair? Only significant results?
p1 <- plot_pair("LMA", "LDMC")   # Regression lines for all species.
p2 <- plot_pair("LMA", "CHL")    # Regression lines only for species other than B. officinalis and H. vulgare.
p3 <- plot_pair("LDMC", "CHL")   # Regression lines only for species other than B. officinalis and H. vulgare.

# Arrange plots
gridExtra::grid.arrange(p1, p2, p3, ncol = 2)


# 4. Growth-trait relationships ----
### SMA regressions ----
library(smatr)
sma4 <- sma(grt_g_d ~ LMA * species, data = df, log = "xy")
summary(sma4) # B. officinalis and H. vulgare not significant
sma5 <- sma(grt_g_d ~ LDMC * species, data = df, log = "xy")
summary(sma5) # B. officinalis and H. vulgare not significant
sma6 <- sma(grt_g_d ~ CHL * species, data = df, log = "xy")
summary(sma6) # All significant

### Plots ----
# Make plots, including regression lines only when significant
# Load required packages
library(ggplot2)
library(ggpmisc)   # Provides stat_ma_line() for Model II SMA lines
library(dplyr)
library(gridExtra)
library(scales)    # For hue_pal()

# A function to create a plot for growth rate vs. a given trait.
# Both axes are log-scaled (without transforming the underlying data).
# A separate SMA regression line is fit for each species,
# but for grt_g_d ~ LMA and grt_g_d ~ LDMC, lines are omitted for 
# B. officinalis and H. vulgare.
plot_growth_vs_trait <- function(trait) {
  p <- ggplot(df, aes_string(x = trait, y = "grt_g_d", color = "species")) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = trait,
         y = "Growth rate (g/d)") +
    theme_bw()
  
  # For LMA and LDMC, omit regression lines for B. officinalis and H. vulgare;
  # for other traits (e.g., CHL), add regression lines for all species.
  if (trait %in% c("LMA", "LDMC")) {
    p <- p + stat_ma_line(
      data = subset(df, !(species %in% c("B. officinalis", "H. vulgare"))),
      mapping = aes_string(x = trait, y = "grt_g_d", group = "species"),
      method = "SMA",
      se = FALSE,
      linewidth = 1
    )
  } else {
    p <- p + stat_ma_line(aes_string(group = "species"),
                          method = "SMA",
                          se = FALSE,
                          linewidth = 1)
  }
  
  return(p)
}

# Create the three plots:
p_LMA  <- plot_growth_vs_trait("LMA")   # grt_g_d vs. LMA: regression lines omitted for B. officinalis and H. vulgare
p_LDMC <- plot_growth_vs_trait("LDMC")  # grt_g_d vs. LDMC: regression lines omitted for B. officinalis and H. vulgare
p_CHL  <- plot_growth_vs_trait("CHL")   # grt_g_d vs. CHL: regression lines for all species

# Arrange the plots in a matrix-style layout (e.g., 2 columns)
gridExtra::grid.arrange(p_LMA, p_LDMC, p_CHL, ncol = 2)


### Model selection ----
#### Simple bivariate models ----

# Load required packages
library(lmerTest)    # for lmer() with p-values
library(MuMIn)       # for r.squaredGLMM()
library(dplyr)

# The three traits to test.
traits <- c("LMA", "LDMC", "CHL")

# Three random-effects structures by species:
# 1) Random intercept
# 2) Random slope
# 3) Random intercept + slope
random_effects_types <- c("intercept", "slope", "intercept+slope")

# Initialize a list to store model-fitting results.
results_list <- list()

# Loop over each trait and each random-effects structure.
for (trait in traits) {
  for (re_type in random_effects_types) {
    
    # Build the random-effects string based on the type.
    if (re_type == "intercept") {
      # Random intercept only
      re_str <- "(1 | species)"
    } else if (re_type == "slope") {
      # Random slope only (no random intercept)
      # i.e. (0 + log10(trait) | species)
      re_str <- paste0("(0 + log10(", trait, ") | species)")
    } else if (re_type == "intercept+slope") {
      # Random intercept + random slope
      # i.e. (log10(trait) | species)
      re_str <- paste0("(log10(", trait, ") | species)")
    }
    
    # Build the full model formula:
    # log10(grt_g_d) ~ log10(trait) + random effects
    formula_str <- paste0("log10(grt_g_d) ~ log10(", trait, ") + ", re_str)
    formula_obj <- as.formula(formula_str)
    
    # Attempt to fit the model with tryCatch (to skip errors).
    model_fit <- tryCatch(
      {
        lmer(formula_obj, data = df)
      },
      error = function(e) {
        # If model fitting fails, return NULL
        NULL
      }
    )
    
    # If the model failed, store NA rows and continue.
    if (is.null(model_fit)) {
      results_list[[length(results_list) + 1]] <- data.frame(
        trait           = trait,
        random_effects  = re_type,
        singular        = NA,
        AIC             = NA,
        p_value         = NA,
        r2_marginal     = NA,
        r2_conditional  = NA,
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Check for singular fit.
    sing <- isSingular(model_fit, tol = 1e-4)
    
    # Extract AIC
    aic_val <- AIC(model_fit)
    
    # Extract p-value for the fixed effect "log10(trait)"
    # (should be the second row in summary, but safer to find by name)
    coefs <- summary(model_fit)$coefficients
    row_name <- paste0("log10(", trait, ")")
    if (row_name %in% rownames(coefs)) {
      p_val <- coefs[row_name, "Pr(>|t|)"]
    } else {
      p_val <- NA
    }
    
    # Compute marginal and conditional R^2 (MuMIn)
    r2_vals <- r.squaredGLMM(model_fit)
    r2_marg <- r2_vals[1]
    r2_cond <- r2_vals[2]
    
    # Store the results
    results_list[[length(results_list) + 1]] <- data.frame(
      trait           = trait,
      random_effects  = re_type,
      singular        = sing,
      AIC             = aic_val,
      p_value         = p_val,
      r2_marginal     = r2_marg,
      r2_conditional  = r2_cond,
      stringsAsFactors = FALSE
    )
  }
}

# Combine all results into one final data frame.
results_df <- do.call(rbind, results_list)
print(results_df)

# Subset only non-singular fits:
results_df <- subset(results_df, singular == FALSE)
# AIC prefers slope for LMA, slope or intercept for LDMC, and intercept for CHL.

#### Multivariate models ----
# Load required packages
library(lmerTest)   # For lmer() with p-values
library(MuMIn)      # For r.squaredGLMM()
library(dplyr)

# Define the three random-effect options
re_options <- c("intercept", "slope", "intercept+slope")

# Create a dataframe with all combinations for LMA, LDMC, CHL.
combo_df <- expand.grid(LMA_opt = re_options,
                        LDMC_opt = re_options,
                        CHL_opt = re_options,
                        stringsAsFactors = FALSE)

# Initialize a list to store results.
results_list <- list()

# Loop over each combination (27 total).
for(i in seq_len(nrow(combo_df))) {
  # Extract the options for each predictor
  opt_LMA   <- combo_df$LMA_opt[i]
  opt_LDMC  <- combo_df$LDMC_opt[i]
  opt_CHL   <- combo_df$CHL_opt[i]
  
  # Build the random effects terms.
  re_terms <- c()
  intercept_flag <- FALSE
  
  # For LMA:
  if(opt_LMA == "slope") {
    re_terms <- c(re_terms, "(0 + log10(LMA)|species)")
  } else if(opt_LMA == "intercept+slope") {
    re_terms <- c(re_terms, "(log10(LMA)|species)")
  } else if(opt_LMA == "intercept") {
    intercept_flag <- TRUE
  }
  
  # For LDMC:
  if(opt_LDMC == "slope") {
    re_terms <- c(re_terms, "(0 + log10(LDMC)|species)")
  } else if(opt_LDMC == "intercept+slope") {
    re_terms <- c(re_terms, "(log10(LDMC)|species)")
  } else if(opt_LDMC == "intercept") {
    intercept_flag <- TRUE
  }
  
  # For CHL:
  if(opt_CHL == "slope") {
    re_terms <- c(re_terms, "(0 + log10(CHL)|species)")
  } else if(opt_CHL == "intercept+slope") {
    re_terms <- c(re_terms, "(log10(CHL)|species)")
  } else if(opt_CHL == "intercept") {
    intercept_flag <- TRUE
  }
  
  # If any predictor was set to "intercept", add a single random intercept.
  if(intercept_flag) {
    re_terms <- c(re_terms, "(1|species)")
  }
  
  # Remove duplicate intercept terms.
  re_terms <- unique(re_terms)
  
  # Build the full random effects string and the full model formula.
  if(length(re_terms) > 0) {
    re_str <- paste(re_terms, collapse = " + ")
    full_formula_str <- paste0("log10(grt_g_d) ~ log10(LMA) + log10(LDMC) + log10(CHL) + ", re_str)
  } else {
    full_formula_str <- "log10(grt_g_d) ~ log10(LMA) + log10(LDMC) + log10(CHL)"
  }
  
  # Convert the string to a formula.
  model_formula <- as.formula(full_formula_str)
  
  # Fit the model, using suppressWarnings() to avoid issues with warnings.
  fit <- try(suppressWarnings(lmer(model_formula, data = df)), silent = TRUE)
  
  # If the model failed, record NA's.
  if(inherits(fit, "try-error")) {
    results_list[[length(results_list) + 1]] <- data.frame(
      LMA_opt = opt_LMA,
      LDMC_opt = opt_LDMC,
      CHL_opt = opt_CHL,
      re_formula = full_formula_str,
      AIC = NA,
      singular = NA,
      p_logLMA = NA,
      p_logLDMC = NA,
      p_logCHL = NA,
      r2_marginal = NA,
      r2_conditional = NA,
      stringsAsFactors = FALSE
    )
    next
  }
  
  # Check for singular fit.
  is_sing <- isSingular(fit, tol = 1e-4)
  
  # Extract AIC.
  aic_val <- AIC(fit)
  
  # Get summary and extract fixed-effect coefficients.
  summ <- summary(fit)
  coefs <- summ$coefficients
  p_logLMA  <- if("log10(LMA)" %in% rownames(coefs))  coefs["log10(LMA)", "Pr(>|t|)"]  else NA
  p_logLDMC <- if("log10(LDMC)" %in% rownames(coefs)) coefs["log10(LDMC)", "Pr(>|t|)"] else NA
  p_logCHL  <- if("log10(CHL)" %in% rownames(coefs))  coefs["log10(CHL)", "Pr(>|t|)"]  else NA
  
  # Compute marginal and conditional R2 using MuMIn.
  r2_vals <- try(r.squaredGLMM(fit), silent = TRUE)
  if(inherits(r2_vals, "try-error")) {
    r2_marg <- NA
    r2_cond <- NA
  } else {
    r2_marg <- r2_vals[1]
    r2_cond <- r2_vals[2]
  }
  
  # Store results in a data frame row.
  results_list[[length(results_list) + 1]] <- data.frame(
    LMA_opt = opt_LMA,
    LDMC_opt = opt_LDMC,
    CHL_opt = opt_CHL,
    re_formula = full_formula_str,
    AIC = aic_val,
    singular = is_sing,
    p_logLMA = p_logLMA,
    p_logLDMC = p_logLDMC,
    p_logCHL = p_logCHL,
    r2_marginal = r2_marg,
    r2_conditional = r2_cond,
    stringsAsFactors = FALSE
  )
}

# Combine all results into one dataframe.
results_df <- do.call(rbind, results_list)
print(results_df)
results_df <- subset(results_df, results_df$singular == "FALSE")

# AIC cannot distinguish between three models. Let's go with random intercepts
# for CHL and LDMC, and random slope for LMA, though alternative models could 
# also be used. NOTE: Importantly, these results are the same as for the simple
# bivariate relationships above, which is great!


# 5. SEM analysis ----

library(piecewiseSEM)

### SEM v0 (starting model) ----
# For direct paths from nitrogen to traits, specify 
# log-transformed traits and growth rate and random intercept by species (see above):
m1.1 <- lmer(log10(LMA)  ~ treatment_mmol + (1 | species), data = df)
m1.2 <- lmer(log10(LDMC)  ~ treatment_mmol + (1 | species), data = df)
m1.3 <- lmer(log10(CHL) ~ treatment_mmol  + (1 | species), data = df)

# Specify a linear model with the following random effects (see above):
# - random intercepts for treatment_mmol, LDMC, and CHL
# - random slope for LMA
m1.4 <- lmer(
  log10(grt_g_d) ~ log10(LMA) + log10(LDMC) + log10(CHL) + treatment_mmol
  + (1 | species)                     # random intercept for species
  + (0 + log10(LMA) | species),       # random slope by species for LMA
  data = df
)


# Combine models into a piecewise SEM
sem1.1 <- psem(
  m1.1,
  m1.2,
  m1.3,
  m1.4,
  log10(LMA) %~~% log10(LDMC), # Note that these are correlations that should have double-headed arrow in path diagram
  log10(LMA) %~~% log10(CHL),
  log10(LDMC) %~~% log10(CHL)
)

# View summary
summary(sem1.1)
# Conclusions:
# 1)  p = 1 indicates the model is saturated, so there are no additional 
#     independence claims to test and no mismatch to measure
# 2) There are 3 non-significant paths, which could be removed in stepwise fashion
# 3) AIC = -505.700


### SEM v1 ----

# Remove most non-significant relationship: log10(LMA) ~ log10(CHL)

# Combine models into a piecewise SEM
sem1.2 <- psem(
  m1.1,
  m1.2,
  m1.3,
  m1.4,
  log10(LMA) %~~% log10(LDMC),
  log10(LDMC) %~~% log10(CHL)
)

# View summary
summary(sem1.2)
# Conclusions:
# 1)  p = 1 indicates the model is saturated, so there are no additional 
#     independence claims to test and no mismatch to measure
# 2) There are 2 non-significant paths, which could be removed in stepwise fashion
# 3) AIC = -505.700
# 4) The warning message can be ignored; this is happening because chi-square
#    can be negative or zero when model is saturated.

### SEM v2 ----

# Remove most non-significant relationship: log10(LDMC) ~ log10(CHL)

# Combine models into a piecewise SEM
sem1.3 <- psem(
  m1.1,
  m1.2,
  m1.3,
  m1.4,
  log10(LMA) %~~% log10(LDMC)
)

# View summary
summary(sem1.3)
# Conclusions:
# 1) Fisher's C and p value suggest model fits well
# 2) The singular fit warning can be ignored. It doesn't actually apply here
#    because isSingular() returns FALSE and the VarCorr() output shows non-zero
#    random effect variances. piecewiseSEM sometimes triggers this warning even 
#    when the model isn't truly singular.
# 3) There is one non-significant path: grt_g_d ~ LMA, which could be removed in stepwise fashion
# 4) AIC = -505.700
# 5) The warning message can be ignored; this is happening because chi-square
#    can be negative or zero when model is saturated.

### SEM v3 ----

# Remove most non-significant relationship: log10(grt_g_d) ~ log10(LMA)

# Specify a linear model with the following random effects (see above):
# - random intercepts for treatment_mmol, LDMC, and CHL
m1.4alt <- lmer(
  log10(grt_g_d) ~ log10(LDMC) + log10(CHL) + treatment_mmol
  + (1 | species),       # random intercepts by species
  data = df
)

# Combine models into a piecewise SEM
sem1.4 <- psem(
  m1.1,
  m1.2,
  m1.3,
  m1.4alt,
  log10(LMA) %~~% log10(LDMC)
)

# View summary
summary(sem1.4)
# Conclusions:
# 1) p values suggest model fits well
# 2) Can ignore singular fit warning as isSingular() returns FALSE
# 3) D-sep indicates missing path from LMA to grt_g_d (which we removed in this 
#    step)
# 4) AIC = -490.128, which lower than sem19. This, combined with d-sep, 
#    indicates that path from LMA to grt_g_d should be re-added

### SEM v4 (final model) ----

# Add back the path from LMA to grt_g_d

sem1.3 <- psem(
  m1.1,
  m1.2,
  m1.3,
  m1.4,
  log10(LMA) %~~% log10(LDMC)
)

# View summary
summary(sem1.3)
# Conclusions:
# 1) Fisher's C, Chi-squared (saturated), and p values indicate model fits well
# 2) The singular fit warning can be ignored. It doesn't actually apply here
#    because isSingular() returns FALSE and the VarCorr() output shows non-zero
#    random effect variances. piecewiseSEM sometimes triggers this warning even 
#    when the model isn't truly singular.
# 3) AIC = -505.700
# 4) The warning message can be ignored; this is happening because chi-square
#    can be negative or zero when model is saturated.
# 6) There is one non-significant path: grt_g_d ~ LMA

# --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --
