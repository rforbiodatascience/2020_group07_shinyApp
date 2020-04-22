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
clinical_data <- read_csv(file = "data/_raw/clinical_data_breast_cancer.csv")
PAM50 <- read_csv(file = "data/_raw/PAM50_proteins.csv")
proteome_data <- read_csv(file = "data/_raw/77_cancer_proteomes_CPTAC_itraq.csv")

# Wrangle data
# ------------------------------------------------------------------------------
## Clean PAM50 data
PAM50_clean <- PAM50 %>%
  select(-Species) %>% # Remove redundant column
  select(RefSeqProteinID, everything()) %>% # ID column first
  rename(RefSeq_accession_number = "RefSeqProteinID",
         gene_symbol = "GeneSymbol",
         gene_name = "Gene Name") # Rename columns to similiar column names as in proteome_data


## Clean proteome data
### Identify patient_ID replicates 
replicates <- colnames(proteome_data) %>% 
  str_replace_all(., '\\..*', '') %>% # Simplify ID name
  duplicated() %>% # Find replicates (true/false)
  colnames(proteome_data)[.] # Extract replicate column names (excluding first apperance)

### Clean file
proteome_data_clean <- proteome_data %>%
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


## Clean clinical data
clinical_data_clean <- clinical_data %>%
  mutate(patient_ID = str_remove_all(.$`Complete TCGA ID`, "TCGA-")) %>% # Simplify ID name
  select(-`Complete TCGA ID`) %>% # Remove old ID column
  semi_join(., proteome_data_clean, by = "patient_ID") %>% # Remove clinical data with no protein information
  select(patient_ID, everything()) # ID column first

names(clinical_data_clean) <- gsub(" ", "_", names(clinical_data_clean)) # Remove whitespaces in column names
names(clinical_data_clean) <- gsub("-", "_", names(clinical_data_clean)) # Remove dash in column names


# Write data
# ------------------------------------------------------------------------------
write_csv(x = clinical_data_clean,
          path = "data/02_clincal_data_clean.csv")

write_csv(x = PAM50_clean,
          path = "data/02_PAM50_clean.csv")

write_csv(x = proteome_data_clean,
          path = "data/02_proteome_data_clean.csv")