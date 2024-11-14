# How are these crops different
# Let's describe them by NUE

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
    treatment_mmol <= 5 ~ "0 - 5 mmol",
    treatment_mmol <= 15 ~ "10 - 15 mmol",
    treatment_mmol <= 25 ~ "20 - 25 mmol",
    treatment_mmol <= 35 ~ "30 - 35 mmol",
  ))

# // -------------------------------