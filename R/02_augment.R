# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Load data
# ------------------------------------------------------------------------------
clinical_data_clean <- read_csv(file = "data/01_clinical_data_clean.csv")
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
  # Rename Class names
  rename(Class = PAM50_mRNA ) %>% 
  mutate(Class = replace_na(Class, "Control"),
         Class = str_replace(Class,
                                  pattern = "HER2-enriched",
                                  replacement = "HER2"),
         Class = str_replace(Class,
                                  pattern = "Luminal A",
                                  replacement = "LumA"),
         Class = str_replace(Class,
                                  pattern = "Luminal B",
                                  replacement = "LumB"),
         Class = str_replace(Class,
                                  pattern = "Basal-like",
                                  replacement = "Basal")) %>%
  # Transform class to factor
  mutate(Class = factor(Class, levels = c("Basal", "HER2", "LumA", "LumB", "Control")))


# Write data
# ------------------------------------------------------------------------------
write_csv(x = joined_data_aug,
          path = "data/02_joined_data_aug.csv")

