# ------------------------------------------------------ #
# Objective: Merge all CSV files together from raw data
#            to produce a single CSV file of fluorometry
#            measurements.
#
# Reference:
#   https://github.com/MichaletzLab/Excision-Photosynthesis-Methods/blob/main/Read.and.Modify.Data.R
#
# @author Nathan Malamud
# @date   2024.05.16
# ------------------------------------------------------ #

# ------------------------------------------------------ #
# LIBRARIES
library(dplyr)
library(purrr)
library(readxl)
library(tidyverse)

# ------------------------------------------------------ #
# Get the files in the directory
# REMEMBER: Set working directory to source file location
files <- list.files(path = "./to_merge", full.names = TRUE)

# ------------------------------------------------------ #
# Process files
processed_data <- lapply(files, function(file) {
  # Read in the file and skip the first header
  data <- read_csv(file, skip = 1)
  
  # Remove the first row containing units
  data <- data[-1, ]
  
  # Return the processed data frame
  return(data)
})

# ------------------------------------------------------ #
# Find common column names
# Check the column names of all data frames in processed_data
all_col_names <- lapply(processed_data, colnames)

# Find the common column names across all data frames
common_col_names <- Reduce(intersect, all_col_names)

# Filter the data frames to keep only the common columns
processed_data_filtered <- lapply(processed_data, function(df) {
  df[, common_col_names, drop = FALSE]
})

# ------------------------------------------------------ #
# Combine into a single data frame
full.exin.data <- do.call(rbind, processed_data_filtered)

# ------------------------------------------------------ #
# Write to CSV
write_csv(full.exin.data, "./fluor_data.csv")
