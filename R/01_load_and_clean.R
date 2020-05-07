# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

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
  select(-replicates) %>% # Remove replicate columns
  # Remove non-PAM50 proteins
  semi_join(PAM50_clean, 
            by = "RefSeq_accession_number") %>%
  # Simplify ID name
  rename_all(funs(stringr::str_replace_all(., 
                                           pattern = '\\..*', 
                                           replacement = ''))) %>% 
  # Make proteins (variables) the columns
  pivot_longer(cols = -c("RefSeq_accession_number"),
               names_to = "patient_ID",
               values_to = "value" ) %>% 
  # Make patient_ID (observations) the rows
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") 

## Clean clinical data
clinical_data_clean <- clinical_data %>% 
  # Change non-syntactic column names
  rename_all(funs(str_replace_all(.,
                                  pattern = c(" "), 
                                  replacement = "_"))) %>%
  rename_all(funs(str_replace_all(.,
                                  pattern = c("-"), 
                                  replacement = "_"))) %>%
  # Simplify ID name
  mutate(patient_ID = str_sub(Complete_TCGA_ID,
                              start = 6, 
                              end = -1)) %>%
  # Remove old ID column
  select(-Complete_TCGA_ID) %>% 
  semi_join(proteome_data_clean, 
            by = "patient_ID") %>% 
  # Remove clinical data with no protein information
  # ID column first
  select(patient_ID, 
         everything()) 



# Write data
# ------------------------------------------------------------------------------
write_csv(x = clinical_data_clean,
          path = "data/01_clinical_data_clean.csv")

write_csv(x = PAM50_clean,
          path = "data/01_PAM50_clean.csv")

write_csv(x = proteome_data_clean,
          path = "data/01_proteome_data_clean.csv")

