######## TO BE MOVED TO 03_EDA ######## 

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
clinical <- read_csv(file = "data/02_clinical_V.csv")


# Preliminary fast hand plots
# ------------------------------------------------------------------------------
clinical %>% 
  ggplot(mapping = aes(clin_Gender, fill = clin_Tumor)) +
  geom_bar()
ggsave("results/gender_vs_tumortype.png", device = "png")    

clinical %>% 
  ggplot(mapping = aes(y = clin_Tumor, x =clin_Metastasis_Coded, colour = clin_Tumor)) +
  geom_jitter(width = 0.1)
ggsave("results/metastatis_vs_tumortype.png", device = "png")

clinical %>% 
  ggplot(mapping = aes(y = clin_Tumor, x =clin_methylation_Clusters, colour = clin_Tumor)) +
  geom_jitter(width = 0.1)
ggsave("results/methylCluster_vs_tumorType.png", device = "png")

clinical %>% 
  ggplot(mapping = aes(y = clin_Tumor, x =clin_Age_at_Initial_Pathologic_Diagnosis, colour = clin_PAM50_mRNA)) +
  geom_jitter(width = 0.1)
ggsave("results/age_vs_cellTyp.png", device = "png")
  
