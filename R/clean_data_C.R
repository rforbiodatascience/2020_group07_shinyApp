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

## Clean clinical data
clinical_data_clean <- clinical_data %>%
  mutate(ID = str_remove_all(.$`Complete TCGA ID`, "TCGA-")) %>% # Make ID similar
  select(-`Complete TCGA ID`) %>% # Remove old column
  select(ID, everything())        # ID column first

### TODO: Drop unused columns
### TODO: Make columns that er female/male positive/negative binary?
### TODO: Remove clinical data with no protein information
### TODO: NAN values (delete obs with na values, or try to estimate it (median value?))

## Clean proteome data
proteome_data_clean <- proteome_data
names(proteome_data_clean)[-1] <- sub("\\..*", "", names(proteome_data_clean)[-1])
### Probably not how it should be done?? Since it doesnt use >%>.
#proteome_data_clean <- proteome_data %>%
  
### TODO: IDs is NOT in long format!
### TODO: NAN values (delete obs with na values, or try to estimate it (median value?))

## Clean PAM50 data
#PAM50_clean <- PAM50 %>%


## join tables? when proteome_data_clean is long as it should, with an ID column.
#joined_data <- left_join(clinical_data_clean, proteome_data_clean, by="ID")


## Useful code for later:
### NA values
#data %>%
#  filter(complete.cases(.))
# negate with ! to check
# delete obs with na in some columns
##data %>%
#  filter(is.na(ID))


# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data_clean,
          path = "data/02_my_data_clean.tsv")
