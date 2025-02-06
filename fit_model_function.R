##' Linear model fitting with mixed effects.
##' 
##' Reference Guide for implementation:
##'   https://www.zoology.ubc.ca/~bio501/R/workshops/lme.html
##'   
##' @author [Nathan D. Malamud]
##' @date [2025-01-29]

library(tidyverse)
library(lmerTest)
library(lme4)
library(MuMIn)
library(visreg)

fit_mixed_models <- function(data, trait, nitrogen, species_list,
                             plot_visreg = FALSE, plot_dir = ".") {
  
  # Filter data for trait value , nitrogen column, and species
  df <- data %>% select(species, all_of(nitrogen), all_of(trait)) %>%
    filter(species %in% species_list)
  
  # Fit models, using REML = FALSE for AIC comparison
  z1 <- lmer(as.formula(paste(trait, "~", nitrogen, "+ (1 | species)")), data = df, REML = FALSE)
  z2 <- lmer(as.formula(paste(trait, "~", nitrogen, "+ (0 +", nitrogen, "| species)")), data = df, REML = FALSE)
  z3 <- lmer(as.formula(paste(trait, "~", nitrogen, "+ (", nitrogen, "| species)")), data = df, REML = FALSE)
  
  models <- list("Random Intercept" = z1,
                 "Random Slope" = z2,
                 "Random Intercept + Slope" = z3)
  
  model_stats <- function(model, model_name) {
    summ <- summary(model)
    
    # Check for singular fit
    singular_fit <- lme4::isSingular(model)
    
    # Extract R² using MuMIn
    r2_vals <- MuMIn::r.squaredGLMM(model)
    
    data.frame(
      Trait = trait,
      Model = model_name,
      AIC = AIC(model),
      Singular = singular_fit,  # Reports TRUE if the model is singular
      R2_marginal = r2_vals[1],  # Fixed effects only
      R2_conditional = r2_vals[2],  # Fixed + random effects
      #Slope = if (nitrogen %in% rownames(summ$coefficients)) summ$coefficients[nitrogen, "Estimate"] else NA,
      P_Value = if (nitrogen %in% rownames(summ$coefficients)) anova(model)$`Pr(>F)`[1] else NA # TODO: check this
    )
    
  }
  
  # If plot_visreg is TRUE, create visreg plots
  # Check if z1, z2, and z3 are properly fitted
  if (plot_visreg) {
    if (!is.null(z1) && inherits(z1, "merMod")) {
      pdf(paste0(plot_dir, "/", trait, "_Random_Intercept.pdf"))
      visreg(z1, by = "species", ylab=trait)
      dev.off()
    }
    if (!is.null(z2) && inherits(z2, "merMod")) {
      pdf(paste0(plot_dir, "/", trait, "_Random_Slope.pdf"))
      visreg(z2, by = "species", ylab=trait)
      dev.off()
    }
    if (!is.null(z3) && inherits(z3, "merMod")) {
      pdf(paste0(plot_dir, "/", trait, "_Random_Intercept_Slope.pdf"))
      visreg(z3, by = "species", ylab=trait)
      dev.off()
    }
  }
  
  stats_table <- bind_rows(lapply(names(models), function(name) model_stats(models[[name]], name)))
  
  return(stats_table)
}
