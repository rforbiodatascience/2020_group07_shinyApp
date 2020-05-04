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
  select(ends_with("CPTAC")) %>%
  colnames(.)

# Remove those patients from the rest - there are duplicates!!!
# df_prot %>%
#  select(-colnames(normal_patients)) %>% # Extract the patient ID and modify
#  rename_all(funs(str_replace_all(., "(\\...TCGA)", ""))) %>%

# Find the duplicates
duplicates <- df_prot %>% 
  colnames(.) %>%
  str_replace_all(., "(\\...TCGA)", "") %>%
  duplicated()



# Remove normal patients and duplicates, truncate the IDs
df_prot_clean <- df_prot %>%
  select(-c(colnames(.)[duplicates],
            gene_name,
            normal_patients)) %>%
  rename(RefSeqProteinID  = "RefSeq_accession_number") %>% 
  rename_all(list(~str_replace(., "(\\...TCGA)", ""))) # something about funs() being depricated soon
# solution was based on this example: Before: funs(name = f(.)) After: list(name = ~f(.))

# save space and remove the old df
rm(df_prot)



# ------------------------------------------------------------------------------
# Clinical Data
# ------------------------------------------------------------------------------
# Truncate the "Complete TCGA ID" column
df_clin1 <- df_clin %>%
  rename(ID = `Complete TCGA ID`) %>%
  mutate(ID = str_replace(ID,"(TCGA-)", ""))

df_clin1 <- df_clin1 %>% 
  rename_all(list(~(str_replace_all(., c(" " = '_',"--" = '_', "-" = '_'))))) %>%
  rename_all(list(~(str_c("clin_",.)))) # maybe to identify this as a separate data source 

write_csv(df_clin1,path = "data/02_clinical_V.csv")
  
  
# Maybe transform all columns into factors
# DF[sapply(DF, is.character)] <- lapply(DF[sapply(DF, is.character)], as.factor)

# ------------------------------------------------------------------------------
# PAM50 Data
# ------------------------------------------------------------------------------

df_PAM_clean <- df_PAM %>% 
  select(-c(Species, `Gene Name`)) %>%
  select(RefSeqProteinID,GeneSymbol) %>% 
  rename(gene_symbol_PAM = GeneSymbol)

# Any missing values  
indx <- apply(df_PAM_clean, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames(df_PAM_clean)[indx]

# save space
rm(df_PAM)
# ------------------------------------------------------------------------------
# JOIN datasets
# ------------------------------------------------------------------------------

# Prot data and PAM
df_prot_PAM <- left_join(df_prot_clean, df_PAM_clean, "RefSeqProteinID") %>%
  select("RefSeqProteinID",
         "gene_symbol", 
         "gene_symbol_PAM", 
         everything())
write_csv(x = df_prot_PAM,
          path = "data/02_prot_PAM_all_V.csv")

df_prot_PAM <- df_prot_PAM%>%
  pivot_longer(data = ., 
               cols = -c("RefSeqProteinID",
                         "gene_symbol",
                         "gene_symbol_PAM"),
               names_to = "patient_ID", 
               values_to =  "values"  ) 

# how to get it into longer format???
#df_prot_transposed %>%
#  pivot_wider(names_from =  "RefSeqProteinID", values_from = "values",values_fill = 0)

