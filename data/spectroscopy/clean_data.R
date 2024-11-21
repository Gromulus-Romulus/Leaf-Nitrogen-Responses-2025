## ------------------------------------------------------ ##
## Objective: Process and clean raw spectroscopy data files
##            to produce a single dataset of corrected 
##            spectroscopy measurements.
##
## Notes:
## - Corrects for sensor overlap (990-1900 nm range).
## - Merges with metadata for additional attributes.
## - Outputs data in both wide and long formats.
##
## Adapted from code written for the 7th Annual Plant 
## Functional Traits Course (PFTC7).
##
## @author [Nicole Bison, Nathan Malamud]
## @date   2024.10.02
##
## Source:
##   https://github.com/MichaletzLab/pftc7_spectroscopy
## ------------------------------------------------------ ##

# LIBRARIES
library(readxl)
library(tidyverse) 
library(dplyr)
library(readr)
library(lubridate)

# This script depends on the "spectrolab" library.
# Source: https://github.com/meireles/spectrolab
# Install using: devtools::install_github("meireles/spectrolab")
library(spectrolab)

# Define useful dplyr "macro" here
`%notin%` <- Negate(`%in%`)

# ------------------------------------------------------ #
# SET WORKING DIRECTORIES //
# Ensure working directory contains folder with the raw data files
INPUT_DIR <- "./raw_spec" 
OUTPUT_DIR <- "../"

if (!file.exists(INPUT_DIR)) {
  print("Please set working directory to source file location in RStudio.")
  stop("Script terminated.")
}

# ------------------------------------------------------ #
# LOAD SPECTROSCOPY FILES //
all.files <- list.files(
  path = INPUT_DIR,
  full.names = TRUE,
  recursive = TRUE,
  pattern = 'sig'
)

# Remove metadata file from paths variable
all.files <- all.files[
  all.files != paste(INPUT_DIR, "/", "/metadata.xlsx")
]

# ------------------------------------------------------ #
# PROCESS SPECTROSCOPY FILES //
spectra_list <- lapply(all.files, function(file) {
  # Read data and import as matrix
  x <- read_spectra(file, format = "sig")
  
  # Correct for sensor overlap in 990-1900 range
  x2 <- match_sensors(x, c(990, 1900)) %>% as.data.frame()
  
  # Convert x to data frame
  x <- as.matrix(x) %>% as.data.frame()
  
  # Extract barcodeID and other info using RegEx
  pattern <- "(?<date>[0-9]{8})_(?<barcodeID>.*?)_(?<scan>[0-9]{4})\\.sig"
  match <- str_match(rownames(x)[1], pattern) %>% as.data.frame()
  
  x2$date <- match$date[1] %>% as.character()
  x2$barcodeID <- match$barcodeID[1] %>% as.character()
  x2$scan <- match$scan[1] %>% as.character()
  
  return(x2)
})

# Filter columns to keep only common column names
all_col_names <- lapply(spectra_list, colnames)
common_col_names <- Reduce(intersect, all_col_names)
processed_data_filtered <- lapply(
  spectra_list, 
  function(df) df[, common_col_names, drop = FALSE]
)

# Combine into single data frame
all.spectra <- do.call(rbind, processed_data_filtered) %>%
  select(barcodeID, scan, date, everything())

# ------------------------------------------------------ #
# DATA CLEANING + AGGREGATION //
# Remove white references and practice test scans
all.spectra$barcodeID <- toupper(all.spectra$barcodeID)
all.spectra <- all.spectra %>%
  filter(barcodeID %notin% c("WR", "WR_", "TEST", "SIG", "BARCODE"))

# Convert dates to standard format
all.spectra$date <- ymd(all.spectra$date)

# Remove specific barcode ID based on metadata file
all.spectra <- all.spectra %>% filter(barcodeID != "27826")

# Average spectra for all scans per barcode ID
all.spectra <- all.spectra %>%
  select(-c(scan, sample_name)) %>% 
  group_by(barcodeID, date) %>%
  summarize(across(everything(), mean)) %>%
  ungroup()

# Merge spectral measurements with sample metadata
metadata <- read.csv("./metadata.csv")
all.spectra <- merge(metadata, all.spectra, by = "barcodeID")

# ------------------------------------------------------ #
# WRITE OUTPUTS (csv format)
all.spectra.wide <- all.spectra
all.spectra.long <- all.spectra %>%
  pivot_longer(
    cols = -c(colnames(metadata), date),
    names_to = "wavelength",
    values_to = "reflectance"
  )

# Output as wide or long table format depending on preference
write_csv(all.spectra.wide, paste(OUTPUT_DIR, "spec_data_wide.csv", sep = ""))
write_csv(all.spectra.long, paste(OUTPUT_DIR, "spec_data_long.csv", sep = ""))
