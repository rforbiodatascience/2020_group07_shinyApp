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
## Clean proteome data

### The IDs needs cleaning up:
#names(proteome_data)[-1] <- sub("\\..*", "", names(proteome_data)[-1]) # Make ID simple
### Probably not how it should be done?? Since it doesnt use >%>.
### Also due to replicates in patient_IDs, this makes the next step fail. 
### So it is commented out for now


proteome_data_clean <- proteome_data %>%
  select(-gene_symbol, -gene_name) %>% # Remove redundant columns
  pivot_longer(cols = -c("RefSeq_accession_number"),
               names_to = "patient_ID",
               values_to = "value") %>% # Make tibble long
  mutate(RefSeq_accession_number = factor(RefSeq_accession_number)) # Factor ID

# Should NAs be replaced with zeroes? Or maybe estimate them with median values?
colSums(is.na(proteome_data_clean))

### TODO: Drop unused columns (We only know this for sure at the end of the analysis)

## Clean clinical data
clinical_data_clean <- clinical_data %>%
  mutate(patient_ID = str_remove_all(.$`Complete TCGA ID`, "TCGA-")) %>% # Make ID simple
  select(-`Complete TCGA ID`) %>% # Remove old ID column
  select(patient_ID, everything())  # ID column first

# There is only NA values in days to date of death for living people
colSums(is.na(clinical_data_clean))

### TODO: Remove clinical data with no protein information
### TODO: Drop unused columns (We only know this for sure at the end of the analysis)
### TODO: Make columns that are female/male and positive/negative binary? (if it make sense)

## Clean PAM50 data
PAM50_clean <- PAM50 %>%
  select(-Species) %>% # Remove redundant column
  mutate(RefSeqProteinID = factor(RefSeqProteinID)) %>% # Factor ID
  select(RefSeqProteinID, everything()) %>% # ID column first
  rename(RefSeq_accession_number = "RefSeqProteinID",
         gene_symbol = "GeneSymbol",
         gene_name = "Gene Name"
         )

# No NA numbers
colSums(is.na(PAM50_clean))

# Include observations missing, which is present in the proteome data.
PAM50_proteome_clean <- proteome_data %>% 
  select(RefSeq_accession_number, gene_symbol, gene_name) 

### I am a bit unsure of the use of PAM50_clean. 
### Since it is incomplete, and the full information is already found in PAM50_proteome.

#________________________________________________________________________

### should we join tables? 
### when proteome_data_clean is long as it should, with an ID column.
#joined_data <- left_join(clinical_data_clean, proteome_data_clean, by="ID")

#________________________________________________________________________
## Useful code for later:
### NA values
#data %>%
#  filter(complete.cases(.))
# negate with ! to check
# delete obs with na in some columns
##data %>%
#  filter(is.na(ID))

### Make something binary
#mutate(gender_bin = case_when(gender == "female" ~ 1,
#                              gender == "male" ~ 0))
#________________________________________________________________________


# Write data
# ------------------------------------------------------------------------------
write_csv(x = clincal_data_clean,
          path = "data/02_clincal_data_clean.csv")

write_csv(x = PAM50_clean,
          path = "data/02_PAM50_clean.csv")

write_csv(x = PAM50_proteome_clean,
          path = "data/02_PAM50_proteome_clean.csv")

write_csv(x = proteome_data_clean,
          path = "data/proteome_data_clean.csv")
