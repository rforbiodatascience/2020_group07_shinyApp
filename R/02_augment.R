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

### Join clinical and proteome data (wide version)
joined_data_wide_aug <- proteome_data_wide_aug %>%
  right_join(clinical_data_clean,. , by="patient_ID") 

### Add PAM50_mRNA class for control samples
joined_data_wide_aug$PAM50_mRNA <- joined_data_wide_aug$PAM50_mRNA %>% replace_na("Control") 

### Make long version of proteome data
proteome_data_long_aug <- proteome_data_wide_aug %>%
  pivot_longer(-patient_ID, names_to = "RefSeq_accession_number", values_to = "value")

### Join clinical and proteome data (long version)
joined_data_long_aug <- proteome_data_long_aug %>%
  full_join(clinical_data_clean,. , by="patient_ID") %>%
  select(patient_ID, RefSeq_accession_number, value, everything())

### Add PAM50_mRNA class for control samples
joined_data_long_aug$PAM50_mRNA <- joined_data_long_aug$PAM50_mRNA %>% replace_na("Control") 


# Write data
# ------------------------------------------------------------------------------
write_csv(x = joined_data_wide_aug,
          path = "data/02_joined_data_wide_aug.csv")

write_csv(x = joined_data_long_aug,
          path = "data/02_joined_data_long_aug.csv")

write_csv(x = proteome_data_wide_aug,
          path = "data/02_proteome_data_wide_aug.csv")

write_csv(x = proteome_data_long_aug,
          path = "data/02_proteome_data_long_aug.csv")

