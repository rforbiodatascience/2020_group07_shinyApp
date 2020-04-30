# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("broom")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
clincal_data_aug <- read_csv(file = "data/01_clinical_data_clean.csv")
PAM50_clean <- read_csv(file = "data/01_PAM50_clean.csv")
joined_data_aug <- read_csv(file = "data/02_joined_data_aug.csv")

# Catrines plots
# ------------------------------------------------------------------------------
ggplot(data = clincal_data_aug) +
  geom_histogram(mapping = aes(x = Age_at_Initial_Pathologic_Diagnosis), binwidth = 5)
ggsave(filename = "results/03_age_distribution.png",device = "png")

ggplot(data = joined_data_aug, mapping = aes(x = NP_057427, y = NP_002408)) + 
  geom_point()
ggsave(filename = "results/03_random_gene_correlation.png",device = "png")

# View class distribution
joined_data_aug %>% count(PAM50_mRNA) %>% print
joined_data_aug %>% count(Tumor) %>% print

# Preliminary fast hand plots
# ------------------------------------------------------------------------------
clincal_data_aug %>% 
  ggplot(mapping = aes(Gender, fill = Tumor)) +
  geom_bar()
ggsave("results/gender_vs_tumortype.png", device = "png")    

clincal_data_aug %>% 
  ggplot(mapping = aes(y = Tumor, x =Metastasis_Coded, colour = Tumor)) +
  geom_jitter(width = 0.1)
ggsave("results/metastatis_vs_tumortype.png", device = "png")

clincal_data_aug %>% 
  ggplot(mapping = aes(y = Tumor, x =methylation_Clusters, colour = Tumor)) +
  geom_jitter(width = 0.1)
ggsave("results/methylCluster_vs_tumorType.png", device = "png")

clincal_data_aug %>% 
  ggplot(mapping = aes(y = Tumor, x =Age_at_Initial_Pathologic_Diagnosis, colour = PAM50_mRNA)) +
  geom_jitter(width = 0.1)
ggsave("results/age_vs_cellTyp.png", device = "png")