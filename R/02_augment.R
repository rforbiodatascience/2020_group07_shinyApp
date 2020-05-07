# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Load data
# ------------------------------------------------------------------------------
clinical_data_clean <- read_csv(file = "data/01_clinical_data_clean.csv")
PAM50_clean <- read_csv(file = "data/01_PAM50_clean.csv")
proteome_data_clean <- read_csv(file = "data/01_proteome_data_clean.csv")

# Wrangle data
# ------------------------------------------------------------------------------
## Handle NA values in proteome data
proteome_data_aug <- proteome_data_clean %>% 
  # Remove genes with too many NAs (over 50%)
  discard(~sum(is.na(.x))/length(.x) >= 0.5) %>% 
  # Take median value to replace NA values
  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) 

## Join clinical and proteome data (wide version)
joined_data_aug <- proteome_data_aug %>%
  right_join(clinical_data_clean, ., by = "patient_ID") %>% 
  # Add PAM50_mRNA class for control samples
  mutate(PAM50_mRNA = replace_na(PAM50_mRNA, "Control"))
  


# Write data
# ------------------------------------------------------------------------------
write_csv(x = joined_data_aug,
          path = "data/02_joined_data_aug.csv")

