# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Load data
# ------------------------------------------------------------------------------
clinical_data_clean <- read_csv(file = "data/01_clincal_data_clean.csv")
PAM50_clean <- read_csv(file = "data/01_PAM50_clean.csv")
proteome_data_wide_clean <- read_csv(file = "data/01_proteome_data_wide_clean.csv")


# Wrangle data
# ------------------------------------------------------------------------------
## Handle NA values in proteome data
proteome_data_wide_aug = proteome_data_wide_clean %>% 
  discard(~sum(is.na(.x))/length(.x) >= 0.5) %>% # Remove genes with too many na's (over 50%)
  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) # Take median value to replace NA values

### Join clinical and proteome data
joined_data_aug <- proteome_data_wide_clean %>%
  left_join(clinical_data_clean, ., by="patient_ID")


### Make long version of proteome data
proteome_data_long_aug <- proteome_data_wide_aug


# Write data
# ------------------------------------------------------------------------------
write_csv(x = joined_data_aug,
          path = "data/02_joined_data_aug.csv")

write_csv(x = proteome_data_wide_aug,
          path = "data/02_proteome_data_wide_aug.csv")

write_csv(x = proteome_data_long_aug,
          path = "data/02_proteome_data_long_aug.csv")

