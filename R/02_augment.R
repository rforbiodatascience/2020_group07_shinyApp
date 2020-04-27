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
proteome_data_clean <- read_csv(file = "data/01_proteome_data_clean.csv")


# Wrangle data
# ------------------------------------------------------------------------------
## Handle NA values in proteome data
proteome_data_aug = proteome_data_clean %>% 
  discard(~sum(is.na(.x))/length(.x) >= 0.5) %>% # Remove genes with too many na's (over 50%)
  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) # Take median value to replace NA values

### Join clinical and proteome data (wide version)
joined_data_aug <- proteome_data_aug %>%
  right_join(clinical_data_clean,. , by="patient_ID") 

### Add PAM50_mRNA class for control samples
joined_data_aug$PAM50_mRNA <- joined_data_aug$PAM50_mRNA %>% replace_na("Control") 


# Write data
# ------------------------------------------------------------------------------
write_csv(x = joined_data_wide_aug,
          path = "data/02_joined_data_aug.csv")

