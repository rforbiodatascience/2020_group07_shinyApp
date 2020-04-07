library("tidyverse")

# Load raw data
setwd("../data/_raw")
metadata_raw = read_csv(file = "clinical_data_breast_cancer.csv")
proteome_raw = read_csv(file = "77_cancer_proteomes_CPTAC_itraq.csv")


# Clean metadata
metadata <- metadata_raw %>% 
  mutate(ID = gsub("TCGA-", "", `Complete TCGA ID`)) %>% 
  select(-`Complete TCGA ID`) %>% # change name of patients' ID
  select(ID, everything())        # put ID column in first position


# Clean proteome
proteome <- proteome_raw %>% 
  select(-gene_symbol, -gene_name) %>%    # remove gene info
  select(-c(ends_with("CPTAC")))          # remove healthy patients 

# identify replicates
replicates <- colnames(proteome) %>% 
  str_replace_all(., '\\...TCGA', '') %>% 
  sort() %>%  
  duplicated()
cols_proteome <- colnames(proteome) %>% 
  sort()

proteome <- proteome %>% 
  select (-cols_proteome[replicates]) %>%    # remove replicates
  rename_all(funs(stringr::str_replace_all(., '\\...TCGA', '')))   # change ID patients

# Transpose proteome
proteome_transposed <- proteome %>%
  gather(ID, value, -RefSeq_accession_number) %>% 
  spread(RefSeq_accession_number, value)


# Join two data sets
data <- inner_join(metadata, proteome_transposed, by="ID")


# Keep gene data - just in case
protein_data <- proteome_raw %>% 
  select(RefSeq_accession_number, gene_symbol, gene_name)


# Save clean data
setwd("../")
save(data, file = "clean_data_P.RData")
save(protein_data, file = "protein_data_P.RData")
