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
proteome_data_clean <- read_csv(file = "data/01_proteome_data_clean.csv")


# Wrangle data
# ------------------------------------------------------------------------------
### Clean file
proteome_data_aug <- proteome_data_clean %>%
  select(-gene_symbol, -gene_name) %>% # Remove redundant columns
  select(-replicates)  %>% # Remove replicate columns
  select(-c(ends_with("CPTAC")))  %>%  # Remove healthy patients (we have no clincal information on them)
  #semi_join(., PAM50_clean, by = "RefSeq_accession_number") %>% # Remove non-PAM50 proteins
  rename_all(funs(stringr::str_replace_all(., '\\..*', ''))) #%>% # Simplify ID name
#pivot_longer(cols = -c("RefSeq_accession_number"),
#             names_to = "patient_ID",
#             values_to = "value") %>% # Make proteins (variables) the columns
#pivot_wider(names_from = "RefSeq_accession_number",
#            values_from = "value") #%>% # Make patient_ID (observations) the rows



### Join clinical and proteome data
joined_data <- left_join(clinical_data_clean, proteome_data_clean, by="patient_ID")


# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data_clean_aug,
          path = "data/03_my_data_clean_aug.tsv")