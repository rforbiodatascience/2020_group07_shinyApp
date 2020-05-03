# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# if (!requireNamespace("clusterProfiler", quietly = TRUE))
#   install.packages("clusterProfiler")

# The package is not available for R version 3.6
# library(clusterProfiler)

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
joined_data_aug <- read_csv(file = "data/02_joined_data_aug.csv")
df_prot_raw <- read_csv(file = "data/_raw/77_cancer_proteomes_CPTAC_itraq.csv")
df_clinical <- read_csv(file = "data/01_clinical_data_clean.csv")
PAM50 <- read_csv(file = "data/01_PAM50_clean.csv")
# df_prot <- read_csv(file = "data/01_proteome_data_clean.csv")


# Data wrangling
# ------------------------------------------------------------------------------
# COmbine the full gene set with clinical data
## Handle NA values in proteome data
df_prot_raw <- df_prot_raw %>% 
                      discard(~sum(is.na(.x))/length(.x) >= 0.5) %>% # Remove genes with too many na's (over 50%)
                      mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) # Take median value to replace NA values

# Identify patient_ID replicates 
replicates <- colnames(df_prot_raw) %>% 
  str_replace_all(., '\\..*', '') %>% # Simplify ID name
  duplicated() %>% # Find replicates (true/false)
  colnames(df_prot_raw)[.] # Extract replicate column names (excluding first apperance)

# Prepare proteome df
proteome_data_clean <- df_prot_raw %>%
  select(-gene_symbol, -gene_name) %>% # Remove redundant columns
  select(-replicates)  %>% # Remove replicate columns
  rename_all(funs(stringr::str_replace_all(., '\\..*', ''))) %>% # Simplify ID name
  pivot_longer(cols = -c("RefSeq_accession_number"),
               names_to = "patient_ID",
               values_to = "value") %>% # Make proteins (variables) the columns
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") # Make patient_ID (observations) the rows

# Join clinical and proteome data (wide version with ALL genes)
df_big <- proteome_data_clean %>%
  right_join(df_clinical,. , by="patient_ID") 




# Data split for each cancer group: without CONTROL samples
df_small <- joined_data_aug %>%
                  select(starts_with("NP_"),
                         PAM50_mRNA, patient_ID) %>%
                  filter(PAM50_mRNA != "Control") 


# For the SMALL df
# Subset for each group extracting the mean across all patients and then sorted absolute values
top_genes = 30

df_HER2 <- df_small %>% return_top_genes(filter_by = "HER2-enriched", new_label = "HER2" , n = top_genes)
                  
df_Basal <- df_small %>% return_top_genes(filter_by = "Basal-like", new_label = "Basal", n = top_genes)

df_LumA <- df_small %>% return_top_genes(filter_by = "Luminal A", new_label = "LumA", n = top_genes)

df_LumB <- df_small %>% return_top_genes(filter_by = "Luminal B", new_label = "LumB", n = top_genes)

# Find a common set between all top genes in each PAM5_mRNA group
top_gene_intersection <- intersect(intersect(df_HER2,df_LumB ),intersect(df_LumA,df_Basal))

# Want to map them to gene symbols (PANTHER GO did not recongnize some of the transcript names)
top_gene_names <- PAM50 %>% rename(new_label = RefSeq_accession_number) %>%
                  left_join(top_gene_intersection,.) %>%
                  select(gene_symbol)

# Output results to usego terms interface
write_csv(x = top_gene_names,
          path = paste0("results/03_GOTerms_top",top_genes,"_intersection.csv"))


# For the BIG df
# Subset for each group extracting the mean across all patients and then sorted absolute values
top_genes2 = 50

df_HER2_full <- df_big %>% return_top_genes(filter_by = "HER2-enriched", new_label = "HER2" , n = top_genes2)

df_Basal_full <- df_big %>% return_top_genes(filter_by = "Basal-like", new_label = "Basal", n = top_genes2)

df_LumA_full <- df_big %>% return_top_genes(filter_by = "Luminal A", new_label = "LumA", n = top_genes2)

df_LumB_full <- df_big %>% return_top_genes(filter_by = "Luminal B", new_label = "LumB", n = top_genes2)
  

# Find a common set between all top genes in each PAM5_mRNA group
top_gene_intersection_full <- intersect(intersect(df_HER2_full, df_Basal_full),intersect(df_LumA_full,df_LumB_full))


# Any intersection between the two sets ?
intersect(top_gene_intersection, top_gene_intersection_full)

# Want to map them to gene symbols (PANTHER GO did not recongnize some of the transcript names)
# top_gene_names_full <- PAM50 %>% rename(new_label = RefSeq_accession_number) %>%
#                         left_join(top_gene_intersection_full,.) %>%
#                         select(gene_symbol)

# Output results to usego terms interface
write_csv(x = top_gene_intersection_full,
          path = paste0("results//03_GOTerms_top",top_genes2,"_intersection_full.csv"))


#df_HER2_full %>% 