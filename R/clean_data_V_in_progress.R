# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Define functions (DELETE IF NOT USED)
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
df_clin <- read_csv(file = "data/_raw/clinical_data_breast_cancer.csv")
df_PAM <- read_csv(file = "data/_raw/PAM50_proteins.csv")
df_prot <- read_csv(file = "data/_raw/77_cancer_proteomes_CPTAC_itraq.csv")

# Data prep
# ------------------------------------------------------------------------------
# PROTEOME DATA
# ------------------------------------------------------------------------------
# Find NORMAL samples
normal_patients <- df_prot %>%
  select(ends_with("CPTAC"))

# Remove those patients from the rest - there are duplicates!!!
# df_prot %>%
#  select(-colnames(normal_patients)) %>% # Extract the patient ID and modify
#  rename_all(funs(str_replace_all(., "(\\...TCGA)", ""))) %>%

# Find the duplicates
duplicates <- df_prot %>% 
              colnames(.) %>%
              str_replace_all(., "(\\...TCGA)", "")

# control the duplicates 
df_prot[str_detect(colnames(df_prot), "AO-A12D\\.|C8-A131\\.|AO-A12B\\.")]

# or the tidy way
df_prot_dup <- colnames(df_prot) %>%
                  str_detect(string = ., "AO-A12D\\.|C8-A131\\.|AO-A12B\\.") %>%
                  df_prot[.]

# side by side
df_prot_dup %>% select(sort(colnames(.)))
  
# Remove normal patients and duplicates, truncate the IDs
df_prot_clean <- df_prot %>%
                  select(-colnames(.)[duplicates], -gene_name) %>%
                  select(-colnames(normal_patients)) %>% 
                  rename(RefSeqProteinID  = "RefSeq_accession_number") %>% 
                  rename_all(funs(str_replace(., "(\\...TCGA)", ""))) 
                  
# Try and transpose the data
df_prot_clean %>%
  gather(pat_ID, val, 3:ncol(.)) %>%
  spread(RefSeqProteinID,val)

# Clinical Data
# ------------------------------------------------------------------------------
# Truncate the "Complete TCGA ID" column
df_clin1 <- df_clin %>%
              rename(ID = `Complete TCGA ID`) %>%
              mutate(ID = str_replace(ID,"(TCGA-)", ""))

# Maybe transform all columns into factors
# DF[sapply(DF, is.character)] <- lapply(DF[sapply(DF, is.character)], as.factor)


# PAM50 Data
# ------------------------------------------------------------------------------

df_PAM_clean <- df_PAM %>% 
                select(-c(Species, `Gene Name`)) %>%
                select(RefSeqProteinID,GeneSymbol)

# Any missing values  
indx <- apply(df_PAM_clean, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames(df_PAM_clean)[indx]
