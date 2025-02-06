##' Fit a linear mixed-effects model for each trait + transform combo and growth rate.
##' Total of 1,019 multivariate regression models cross-compared with AIC and R2.
##'
##' Written with ChatGPT-4o assistance
##' Expected run time is 30-45 minutes. Actually appears to take 5
##'
##' @author: [Nathan D. Malamud]
##' @date: [2025-01-31]

library(tidyverse)
library(lme4)
library(MuMIn)
library(performance)  # For R² and model metrics

# Import data
traits <- read_csv("./data/traits.csv")

# Ensure species is a factor
traits$species <- factor(traits$species, 
                         levels = c("R. sativus", "B. officinalis", "H. vulgare"))

# Compute growth rate per day
growth_period_days <- 6 * 7  # 6-week experiment
traits$grt_g_d <- traits$dry_whole_g / growth_period_days

# Filter only by traits of interest
traits <- traits %>% select(species, LDMC, LMA, CHL, grt_g_d)

# Define transformations
apply_transformation <- function(x, type) {
  switch(type,
         # Define just as label for identity function
         "just" = x,
         "log" = log(x),
         "sqrt" = sqrt(x),
         "inv" = 1 / x,
         "inv_sqrt" = 1 / sqrt(x),
         "sq" = x^2,
         "cub" = x^3)
}

# Define transformations
transformations <- c("just", "log", "sqrt", "inv", "inv_sqrt", "sq", "cub")

# Generate all combinations using expand.grid function from base R
# It's essential to do this in order to iterate through all possible models.
transformation_combinations <- expand.grid(
  LDMC = transformations,
  LMA = transformations,
  CHL = transformations
)

# Initialize results storage
all_model_results <- tibble()

# Loop through each transformation combination
for (i in 1:nrow(transformation_combinations)) {
  trans <- transformation_combinations[i, ]
  
  # Apply transformations correctly
  traits_transformed <- traits %>% mutate(
    LDMC_trans = apply_transformation(LDMC, as.character(trans$LDMC)),
    LMA_trans = apply_transformation(LMA, as.character(trans$LMA)),
    CHL_trans = apply_transformation(CHL, as.character(trans$CHL))
  )
  
  # Remove invalid values
  traits_transformed <- traits_transformed %>% filter_all(all_vars(is.finite(.)))
  
  # Skip if data is empty
  if (nrow(traits_transformed) == 0) next
  
  # Loop through growth rate transformations
  for (grt_transform in transformations) {
    traits_transformed <- traits_transformed %>% mutate(
      grt_trans = apply_transformation(grt_g_d, grt_transform)
    )
    
    # Skip if data is empty after transformation
    if (nrow(traits_transformed) == 0) next
    
    # Fit mixed-effects models
    mod_intercept <- lmer(grt_trans ~ LDMC_trans + LMA_trans + CHL_trans + (1 | species), data = traits_transformed)
    mod_slope <- lmer(grt_trans ~ LDMC_trans + LMA_trans + CHL_trans + (LDMC_trans + LMA_trans + CHL_trans | species), data = traits_transformed)
    mod_slope_intercept <- lmer(grt_trans ~ LDMC_trans + LMA_trans + CHL_trans + (1 + LDMC_trans + LMA_trans + CHL_trans | species), data = traits_transformed)
    
    # Compute model metrics
    aic1 <- AIC(mod_intercept)
    aic2 <- AIC(mod_slope)
    aic3 <- AIC(mod_slope_intercept)
    
    # Also BIC
    bic1 <- BIC(mod_intercept)
    bic2 <- BIC(mod_slope)
    bic3 <- BIC(mod_slope_intercept)
    
    # And AICc
    aicc1 <- AICc(mod_intercept)
    aicc2 <- AICc(mod_slope)
    aicc3 <- AICc(mod_slope_intercept)
    
    # 1st is marginal R2, 2nd is conditional R2
    r2_1 <- r.squaredGLMM(mod_intercept)[2]
    r2_2 <- r.squaredGLMM(mod_slope)[2]
    r2_3 <- r.squaredGLMM(mod_slope_intercept)[2]
    
    # Store results
    AICmodel_results <- tibble(
      # Check for singular fit
      LDMC_trans = trans$LDMC, LMA_trans = trans$LMA, CHL_trans = trans$CHL,
      grt_trans = grt_transform,
      Model = c("Random Intercept", "Random Slope", "Random Slope + Intercept"),
      AIC = c(aic1, aic2, aic3),
      BIC = c(bic1, bic2, bic3),
      AICc = c(aicc1, aicc2, aicc3),
      R2 = c(r2_1, r2_2, r2_3),
      Singular_fit = c(isSingular(mod_intercept), isSingular(mod_slope), isSingular(mod_slope_intercept)))
    
    # Create a list to store all models
    models_list <- list(
      mod_intercept = mod_intercept,
      mod_slope = mod_slope,
      mod_slope_intercept = mod_slope_intercept,
      AICmodel_results = AICmodel_results
    )
    
    # Save the list as an RDS file in the appropriate directory
    LDMC_str <- as.character(trans$LDMC[1])
    LMA_str <- as.character(trans$LMA[1])
    CHL_str <- as.character(trans$CHL[1])
    
    # Define the directory path
    dir_path <- paste0("./candidate_mod_objects/", grt_transform, "_grt/")
    
    # Check if the directory exists, and create it if not
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)  # 'recursive = TRUE' ensures nested directories are created
    }
    
    # Save the models list with correct filename formatting
    saveRDS(models_list, paste0(
      dir_path,
      LDMC_str, "LDMC_", LMA_str, "LMA_", CHL_str, "CHL_", "_AICmodels.rds"
    ))
    
    all_model_results <- bind_rows(all_model_results, AICmodel_results)
    
    # Clear memory prior to next iteration
    rm(models_list, mod_intercept, mod_slope, mod_slope_intercept, aic1, aic2, aic3, r2_1, r2_2, r2_3)
  }
}

# View results
print(all_model_results)

# A total of models were fit
all_model_results <- all_model_results %>%
  distinct(LDMC_trans, LMA_trans, CHL_trans, grt_trans,
           Model, Singular_fit, .keep_all = TRUE)

write_csv(all_model_results, "all_AICresults.csv")
