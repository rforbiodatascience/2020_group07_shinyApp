# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())


# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")


# Set working directory
# ------------------------------------------------------------------------------
# setwd("~/DTU_Studies/4-Semester/R_for_Bioscience/Exam_Project/2020_group07/R/")


# Define functions
# ------------------------------------------------------------------------------
# source(file = "R/99_project_functions.R")


# Load data
# ------------------------------------------------------------------------------
# check the datasets in tables
proteome_raw_data <- as_tibble(read_csv(file = "../data/_raw/77_cancer_proteomes_CPTAC_itraq.csv")) # -> rows x cols - 12,553 x 86
clinical_raw_data <- as_tibble(read_csv(file = "../data/_raw/clinical_data_breast_cancer.csv"))  # -> rows x cols - 105 x 30
pam50_raw_data <- as_tibble(read_csv(file = "../data/_raw/PAM50_proteins.csv"))


# Wrangle data
# ------------------------------------------------------------------------------

##### PROTEOMIC DATA #####
# 1. - pivot: sample names spread over in the header -> single column ("sampleID") containing sample names 
#           expression values in cells -> single column ("protein_expression") containing expr.values
#    - create a simplified ID column -> also help to find replicates (duplicated IDs)
# 3. identify replicates -> one patient has been tested for 12553 genes, more occurance implies replicates
# 4. take the average of replicates
# 5. create a data where replicates are replaced with avarge expressions
# 6. drop columns -> gene_symbol and gene_name refer the same, however gene_symbol column is uncomplete
# 7. create a df for healty patient
# 8. diselect healty patients from proteome data

proteome_long <- proteome_raw_data %>%
  pivot_longer(cols = "AO-A12D.01TCGA":"c4155b-C.CPTAC", names_to = "TCGA_sampleID", values_to = "protein_expression") %>%
  mutate(sampleID = str_sub(TCGA_sampleID, 1, 7)) %>% 
  select(-TCGA_sampleID)
  
replicates <- proteome_long %>%
  group_by(sampleID) %>% 
  summarise(n = n()) %>%
  filter(n > 12553) 

replicates_averaged <- proteome_long %>% 
  filter(sampleID %in% replicates$sampleID) %>% 
  group_by(sampleID, RefSeq_accession_number, gene_symbol, gene_name) %>% 
  summarise(average = mean(protein_expression, na.rm = TRUE)) 

replicates_free <- proteome_long %>% 
  filter(!sampleID %in% replicates$sampleID) %>% 
  full_join(replicates_averaged) %>%  # proteome_long without relicates + average of replicates
  mutate(protein_expression = as.character(protein_expression),
         average = as.character(average)) %>% # na.rm removes only chr NAs in the next line
  unite(expression, protein_expression, average, na.rm = TRUE, remove = TRUE) %>% 
  mutate(expression = as.double(expression)) # convert expression values back to doubles, to have the NAs recognizable

proteome_reduced <- replicates_free %>% 
  select(-gene_symbol)

healthy_data <- proteome_reduced  %>% 
  filter(grepl("^.+(-)$", sampleID)) # after renaming the original ID col, healthy patients ended up like this: c4155b-

proteome_cleaned <- anti_join(proteome_reduced, healthy_data)


##### CLINICAL DATA #####
# 1. change spaces to _ in header
# 2. create a simplified ID column -> will used to join datasets

colnames(clinical_raw_data) <- str_replace_all(colnames(clinical_raw_data)," ", "_")

clinical_cleaned <- clinical_raw_data %>%
  mutate(sampleID = str_sub(Complete_TCGA_ID, 6, -1)) %>% 
  select(-Complete_TCGA_ID)
  
  
##### JOIN DATA #####

data_joined <- right_join(clinical_cleaned, proteome_cleaned, by = "sampleID") # this will drop 28 samples from the clinical data which did not pass the quality control according to the paper

# Write data
# ------------------------------------------------------------------------------
write_csv(x = proteome_cleaned, path = "../data/02_proteome_data_clean_Reb.csv")

write_csv(x = healthy_data, path = "../data/02_healthy_data_cleaned_Reb.csv")

write_csv(x = clinical_cleaned, path = "../data/02_clinical_data_cleaned_Reb.csv")

write_csv(x = data_joined, path = "../data/02_joined_data_Reb.csv")
