# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  install.packages("clusterProfiler")

# The package is not available for R version 3.6
# library(clusterProfiler)

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
joined_data_aug <- read_csv(file = "data/02_joined_data_aug.csv")

# Data wrangling
# ------------------------------------------------------------------------------
# Data split for each cancer group: without CONTROL samples
df <- joined_data_aug %>%
                  select(starts_with("NP_"),
                         PAM50_mRNA, patient_ID
                         ) %>%
                  filter(PAM50_mRNA != "Control") 

###################################
# Maybe wrap this in a function !!!
# Subset for each group extracting the mean across all patients and then sorted absolute values
top_genes = 10

df_HER2 <- df %>% return_top_genes(filter_by = "HER2-enriched", new_label = "LumB", n = top_genes)
                  
df_Basal<- df %>% return_top_genes(filter_by = "Basal-like", new_label = "LumB", n = top_genes)

df_LumA <- df %>% return_top_genes(filter_by = "Luminal B", new_label = "LumB", n = top_genes)

df_LumB <- df %>% return_top_genes(filter_by = "Luminal B", new_label = "LumB", n = top_genes)

############ TO CONTINUE : look at VEnn diagram scrip and check how to create sets, join by gene symbol if necesssary
# sets are done through intersect()


x <- df %>%
  return_top_genes(filter_by = "Luminal B", new_label = "LumB", n = 15)
