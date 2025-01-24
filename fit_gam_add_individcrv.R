# Purpose: Fit gam models with all combinations of relevant variables
# Plots:  Plot partial effects for four gam models
# Dependencies: configure.datfile.R
# Outputs: "GAM_Model_Fits.pdf"



# Function to fit GAM models
fit_gam <- function(data, response, predictors) {
  gam_formula <- as.formula(paste(response, "~", paste(sapply(predictors, function(x) paste("s(", x, ", bs='ts', k=5)")), collapse = "+")))
  fit <- gam(gam_formula, data = data, select = TRUE, method = "REML")
  return(fit)
}

# Fit four GAM models for each group or condition
curve_ids <- unique(dat$curveid)
gam_models <- lapply(curve_ids, function(curve) {
  data_subset <- dat[curveid == curve]
  model_tleaf <- fit_gam(data_subset, "photo", "tleaf")
  model_tleaf_cond <- fit_gam(data_subset, "photo", c("tleaf", "cond"))
  model_tleaf_vpdl <- fit_gam(data_subset, "photo", c("tleaf", "vpdl"))
  model_tleaf_cond_vpdl <- fit_gam(data_subset, "photo", c("tleaf", "cond", "vpdl"))
  list(tleaf = model_tleaf, tleaf_cond = model_tleaf_cond, tleaf_vpdl = model_tleaf_vpdl, tleaf_cond_vpdl = model_tleaf_cond_vpdl)
})

# Extract smooth estimates and SE for tleaf for each group and model
smooth_estimates <- lapply(gam_models, function(models) {
  tleaf_estimates <- gratia::smooth_estimates(models$tleaf, unconditional = TRUE, smooth = 's(tleaf)')
  tleaf_cond_estimates <- gratia::smooth_estimates(models$tleaf_cond, unconditional = TRUE, smooth = 's(tleaf)')
  tleaf_vpdl_estimates <- gratia::smooth_estimates(models$tleaf_vpdl, unconditional = TRUE, smooth = 's(tleaf)')
  tleaf_cond_vpdl_estimates <- gratia::smooth_estimates(models$tleaf_cond_vpdl, unconditional = TRUE, smooth = 's(tleaf)')
  list(tleaf = tleaf_estimates, tleaf_cond = tleaf_cond_estimates, tleaf_vpdl = tleaf_vpdl_estimates, tleaf_cond_vpdl = tleaf_cond_vpdl_estimates)
})

# Ensure consistent column names and combine smooth estimates into a single data table with consistent column names
smooth_data <- rbindlist(lapply(seq_along(curve_ids), function(i) {
  curveid <- curve_ids[i]
  tleaf_data <- data.table(curveid = curveid, model = "tleaf", smooth_estimates[[i]]$tleaf)
  tleaf_cond_data <- data.table(curveid = curveid, model = "tleaf_cond", smooth_estimates[[i]]$tleaf_cond)
  tleaf_vpdl_data <- data.table(curveid = curveid, model = "tleaf_vpdl", smooth_estimates[[i]]$tleaf_vpdl)
  tleaf_cond_vpdl_data <- data.table(curveid = curveid, model = "tleaf_cond_vpdl", smooth_estimates[[i]]$tleaf_cond_vpdl)
  
  # Ensure consistent column names
  tleaf_data <- tleaf_data[, .(curveid, model, tleaf, .estimate, .se)]
  tleaf_cond_data <- tleaf_cond_data[, .(curveid, model, tleaf, .estimate, .se)]
  tleaf_vpdl_data <- tleaf_vpdl_data[, .(curveid, model, tleaf, .estimate, .se)]
  tleaf_cond_vpdl_data <- tleaf_cond_vpdl_data[, .(curveid, model, tleaf, .estimate, .se)]
  
  rbind(tleaf_data, tleaf_cond_data, tleaf_vpdl_data, tleaf_cond_vpdl_data, fill = TRUE)
}), fill = TRUE)

# Convert curveid and model to factor
smooth_data[, curveid := factor(curveid, levels = curve_ids)]
smooth_data[, model := factor(model, levels = c("tleaf", "tleaf_cond", "tleaf_vpdl", "tleaf_cond_vpdl"))]

# Function to find the optimal temperature (tleaf at max estimate)
find_optimal_tleaf <- function(smooth_data) {
  optimal_row <- smooth_data[which.max(.estimate)]
  return(optimal_row$tleaf)
}

# Function to calculate breadth
get_breadth <- function(photo, tleaf, threshold = 0.95) {
  # Find the indices where photo drops below the threshold
  idx_left <- max(which(photo >= threshold * max(photo)), default = 1)
  idx_right <- min(which(photo >= threshold * max(photo)), default = length(photo))
  
  # Calculate breadth as absolute difference
  breadth <- abs(tleaf[idx_right] - tleaf[idx_left])
  
  return(list(idx_left = idx_left, idx_right = idx_right, breadth = breadth))
}

# Data frame to store optimal temperatures, breadth, AIC, and concurvity
optimal_temps <- data.table(curveid = character(),
                            tleaf_opt = numeric(),
                            tleaf_cond_opt = numeric(),
                            tleaf_vpdl_opt = numeric(),
                            tleaf_cond_vpdl_opt = numeric(),
                            weib_Topt = numeric(),
                            weib_breadth = numeric(),
                            tleaf_breadth = numeric(),
                            tleaf_cond_breadth = numeric(),
                            tleaf_vpdl_breadth = numeric(),
                            tleaf_cond_vpdl_breadth = numeric(),
                            tleaf_AIC = numeric(),
                            tleaf_cond_AIC = numeric(),
                            tleaf_vpdl_AIC = numeric(),
                            tleaf_cond_vpdl_AIC = numeric(),
                            tleaf_concurvity = numeric(),
                            tleaf_cond_concurvity = numeric(),
                            tleaf_vpdl_concurvity = numeric(),
                            tleaf_cond_vpdl_concurvity = numeric())

# Loop through each curveid
pdf("GAM_Model_Fits.pdf")
for (i in seq_along(curve_ids)) {
  curve_data <- dat[curveid == curve_ids[i]]
  smooth_data_curve <- smooth_data[curveid == curve_ids[i]]
  
  # Calculate breadths and get indices using get_breadth function
  tleaf_breadth_data <- get_breadth(smooth_data_curve[model == "tleaf"]$.estimate, smooth_data_curve[model == "tleaf"]$tleaf)
  tleaf_cond_breadth_data <- get_breadth(smooth_data_curve[model == "tleaf_cond"]$.estimate, smooth_data_curve[model == "tleaf_cond"]$tleaf)
  tleaf_vpdl_breadth_data <- get_breadth(smooth_data_curve[model == "tleaf_vpdl"]$.estimate, smooth_data_curve[model == "tleaf_vpdl"]$tleaf)
  tleaf_cond_vpdl_breadth_data <- get_breadth(smooth_data_curve[model == "tleaf_cond_vpdl"]$.estimate, smooth_data_curve[model == "tleaf_cond_vpdl"]$tleaf)
  
  # Extract breadths
  tleaf_breadth <- tleaf_breadth_data$breadth
  tleaf_cond_breadth <- tleaf_cond_breadth_data$breadth
  tleaf_vpdl_breadth <- tleaf_vpdl_breadth_data$breadth
  tleaf_cond_vpdl_breadth <- tleaf_cond_vpdl_breadth_data$breadth
  
  # Extract unique weib.topt and weib.breadth values for the current curveid
  weib_Topt_value <- unique(curve_data$weib.topt)
  weib_breadth_value <- unique(curve_data$weib.breadth)
  
  # Extract AIC and concurvity values
  model_tleaf <- gam_models[[i]]$tleaf
  model_tleaf_cond <- gam_models[[i]]$tleaf_cond
  model_tleaf_vpdl <- gam_models[[i]]$tleaf_vpdl
  model_tleaf_cond_vpdl <- gam_models[[i]]$tleaf_cond_vpdl
  
  tleaf_AIC <- AIC(model_tleaf)
  tleaf_cond_AIC <- AIC(model_tleaf_cond)
  tleaf_vpdl_AIC <- AIC(model_tleaf_vpdl)
  tleaf_cond_vpdl_AIC <- AIC(model_tleaf_cond_vpdl)
  
  tleaf_concurvity <- concurvity(model_tleaf)
  tleaf_cond_concurvity <- concurvity(model_tleaf_cond)
  tleaf_vpdl_concurvity <- concurvity(model_tleaf_vpdl)
  tleaf_cond_vpdl_concurvity <- concurvity(model_tleaf_cond_vpdl)
  
  # Append optimal temperatures, breadths, AIC, and concurvity to the data frame
  optimal_temps <- rbind(optimal_temps, data.table(
    curveid = curve_ids[i],
    tleaf_opt = find_optimal_tleaf(smooth_data_curve[model == "tleaf"]),
    tleaf_cond_opt = find_optimal_tleaf(smooth_data_curve[model == "tleaf_cond"]),
    tleaf_vpdl_opt = find_optimal_tleaf(smooth_data_curve[model == "tleaf_vpdl"]),
    tleaf_cond_vpdl_opt = find_optimal_tleaf(smooth_data_curve[model == "tleaf_cond_vpdl"]),
    weib_Topt = ifelse(length(weib_Topt_value) > 0, weib_Topt_value, NA),
    weib_breadth = ifelse(length(weib_breadth_value) > 0, weib_breadth_value, NA),
    tleaf_breadth = tleaf_breadth,
    tleaf_cond_breadth = tleaf_cond_breadth,
    tleaf_vpdl_breadth = tleaf_vpdl_breadth,
    tleaf_cond_vpdl_breadth = tleaf_cond_vpdl_breadth,
    tleaf_AIC = tleaf_AIC,
    tleaf_cond_AIC = tleaf_cond_AIC,
    tleaf_vpdl_AIC = tleaf_vpdl_AIC,
    tleaf_cond_vpdl_AIC = tleaf_cond_vpdl_AIC,
    tleaf_concurvity = max(tleaf_concurvity),
    tleaf_cond_concurvity = max(tleaf_cond_concurvity),
    tleaf_vpdl_concurvity = max(tleaf_vpdl_concurvity),
    tleaf_cond_vpdl_concurvity = max(tleaf_cond_vpdl_concurvity)
  ))
  
  # Extract optimal temperatures for the current curve
  current_optimal <- optimal_temps[curveid == curve_ids[i]]
  
  p <- ggplot() +
    geom_line(data = smooth_data_curve[model == "tleaf"], 
              aes(x = tleaf, y = .estimate, color = "tleaf"), size = 1.2) +
    geom_line(data = smooth_data_curve[model == "tleaf_cond"], 
              aes(x = tleaf, y = .estimate, color = "tleaf + cond"), size = 1.2) +
    geom_line(data = smooth_data_curve[model == "tleaf_vpdl"], 
              aes(x = tleaf, y = .estimate, color = "tleaf + vpdl"), size = 1.2) +
    geom_line(data = smooth_data_curve[model == "tleaf_cond_vpdl"], 
              aes(x = tleaf, y = .estimate, color = "tleaf + cond + vpdl"), size = 1.2) +
    geom_ribbon(data = smooth_data_curve[model == "tleaf"], 
                aes(x = tleaf, ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = "tleaf"), alpha = 0.2) +
    geom_ribbon(data = smooth_data_curve[model == "tleaf_cond"], 
                aes(x = tleaf, ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = "tleaf + cond"), alpha = 0.2) +
    geom_ribbon(data = smooth_data_curve[model == "tleaf_vpdl"], 
                aes(x = tleaf, ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = "tleaf + vpdl"), alpha = 0.2) +
    geom_ribbon(data = smooth_data_curve[model == "tleaf_cond_vpdl"], 
                aes(x = tleaf, ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = "tleaf + cond + vpdl"), alpha = 0.2) +
    geom_vline(xintercept = current_optimal$tleaf_opt, color = "blue", linetype = "dashed", size = 1) +
    geom_vline(xintercept = current_optimal$tleaf_cond_opt, color = "deeppink", linetype = "dashed", size = 1) +
    geom_vline(xintercept = current_optimal$tleaf_vpdl_opt, color = "forestgreen", linetype = "dashed", size = 1) +
    geom_vline(xintercept = current_optimal$tleaf_cond_vpdl_opt, color = "purple", linetype = "dashed", size = 1) +
    geom_vline(xintercept = weib_Topt_value, color = "grey", linetype = "dashed", size = 1) +
    labs(title = paste("Curve:", curve_ids[i]), x = "Leaf Temperature (°C)", y = "Partial Effect on Photosynthesis (µmol m^-2 s^-1)") +
    scale_color_manual(name = "Model", 
                       values = c("tleaf" = "blue", "tleaf + cond" = "deeppink", "tleaf + vpdl" = "forestgreen", "tleaf + cond + vpdl" = "purple"),
                       labels = c("tleaf" = "tleaf", "tleaf + cond" = "tleaf + cond", "tleaf + vpdl" = "tleaf + vpdl", "tleaf + cond + vpdl" = "tleaf + cond + vpdl")) +
    scale_fill_manual(name = "Model", 
                      values = c("tleaf" = "blue", "tleaf + cond" = "deeppink", "tleaf + vpdl" = "forestgreen", "tleaf + cond + vpdl" = "purple"),
                      labels = c("tleaf" = "tleaf", "tleaf + cond" = "tleaf + cond", "tleaf + vpdl" = "tleaf + vpdl", "tleaf + cond + vpdl" = "tleaf + cond + vpdl")) +
    theme_minimal() +
    annotate("text", x = Inf, y = Inf, label = paste("tleaf_AIC:", round(tleaf_AIC, 2), "\nconcurvity:", round(max(tleaf_concurvity), 2)), hjust = 1.1, vjust = 2, size = 3.5, color = "blue") +
    annotate("text", x = Inf, y = Inf, label = paste("tleaf_cond_AIC:", round(tleaf_cond_AIC, 2), "\nconcurvity:", round(max(tleaf_cond_concurvity), 2)), hjust = 1.1, vjust = 4, size = 3.5, color = "red") +
    annotate("text", x = Inf, y = Inf, label = paste("tleaf_vpdl_AIC:", round(tleaf_vpdl_AIC, 2), "\nconcurvity:", round(max(tleaf_vpdl_concurvity), 2)), hjust = 1.1, vjust = 6, size = 3.5, color = "green") +
    annotate("text", x = Inf, y = Inf, label = paste("tleaf_cond_vpdl_AIC:", round(tleaf_cond_vpdl_AIC, 2), "\nconcurvity:", round(max(tleaf_cond_vpdl_concurvity), 2)), hjust = 1.1, vjust = 8, size = 3.5, color = "purple")
  
  print(p)
}
dev.off()
write.csv(optimal_temps, "gam_fits_add.csv", row.names = FALSE)
