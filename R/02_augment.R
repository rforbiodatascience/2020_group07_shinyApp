# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
clinical_data_clean <- read_csv(file = "data/01_clincal_data_clean.csv")
PAM50_clean <- read_csv(file = "data/01_PAM50_clean.csv")
proteome_data_wide_clean <- read_csv(file = "data/01_proteome_data_wide_clean.csv")


# Wrangle data
# ------------------------------------------------------------------------------
### Join clinical and proteome data
joined_data <- left_join(clinical_data_clean, proteome_data_wide_clean, by="patient_ID")



### Make long version of proteome data



# Write data
# ------------------------------------------------------------------------------
write_csv(x = joined_data,
          path = "data/02_joined_data_clean.csv")