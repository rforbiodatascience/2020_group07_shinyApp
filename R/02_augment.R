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
PAM50_clean <- read_csv(file = "data/01_PAM50_clean.csv")

# Wrangle data
# ------------------------------------------------------------------------------
# Rename Class groups in clinical data
clinical_data_aug <- clinical_data_clean %>% 
  rename(Class = PAM50_mRNA ) %>% 
  mutate(Class = str_replace(Class,
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
  mutate(Class = factor(Class, 
                        levels = c("Basal", "HER2", "LumA", "LumB", "Control")))


## Handle NA values in proteome data
proteome_data_aug <- proteome_data_clean %>% 
  # Remove genes with too many NAs (over 50%)
  discard(~sum(is.na(.x))/length(.x) >= 0.5) %>% 
  # Take median value to replace NA values
  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) 


## Create a PAM50genes-filtered version of the proteome data
PAM50genes <- PAM50_clean %>% 
  select(RefSeq_accession_number) %>% 
  pull()
  
proteome_data_PAM50_aug <- proteome_data_aug %>% 
  select(any_of(PAM50genes)) %>% 
  mutate(patient_ID = proteome_data_aug$patient_ID)



## Join clinical and proteome data (full version)
joined_data_full_aug <- proteome_data_aug %>%
  right_join(clinical_data_aug, 
             ., 
             by = "patient_ID") %>% 
  # Add control labels
  mutate(Class = replace_na(Class, 
                            replace = "Control"))


## Join clinical and proteome data (PAM50genes filtered version)
joined_data_PAM50_aug <- proteome_data_PAM50_aug %>%
  right_join(clinical_data_aug, 
             ., 
             by = "patient_ID")  %>% 
  # Add control labels
  mutate(Class = replace_na(Class, "Control")) %>% 
  mutate(Class = factor(Class, levels = c("Basal", "HER2", "LumA", "LumB", "Control")))


# Write data
# ------------------------------------------------------------------------------
## Full version
write_csv(x = joined_data_full_aug,
          path = "data/02_joined_data_full_aug.csv")

## PAM50genes filtered version
write_csv(x = joined_data_PAM50_aug,
          path = "data/02_joined_data_PAM50_aug.csv")
