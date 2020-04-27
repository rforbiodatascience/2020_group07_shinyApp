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
clincal_data_aug <- read_csv(file = "data/01_clincal_data_clean.csv")
PAM50_aug <- read_csv(file = "data/01_PAM50_clean.csv")
proteome_data_aug <- read_csv(file = "data/02_proteome_data_wide_aug.csv")
joined_data_aug <- read_csv(file = "data/02_joined_data_aug.csv")

# Wrangle data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...

# Model data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
ggplot(data = clincal_data_aug) +
  geom_histogram(mapping = aes(x = Age_at_Initial_Pathologic_Diagnosis), binwidth = 5)
ggsave(filename = "results/03_age_distribution.png",device = "png")

ggplot(data = proteome_data_aug, mapping = aes(x = NP_057427, y = NP_002408)) + 
  geom_point()
ggsave(filename = "results/03_random_gene_correlation.png",device = "png")

# View class distribution
joined_data_aug %>% count(PAM50_mRNA) %>% print
joined_data_aug %>% count(Tumor) %>% print

# PCA ---------------------------------------------------------------------------

proteome_pca <- proteome_data_aug %>%
  select(starts_with("NP"))  %>%
  prcomp(center = TRUE, scale = TRUE) 
  
proteome_pca %>%
  tidy("pcs") %>% 
  ggplot(aes(x = PC, y = percent)) +
  geom_col() +
  theme_bw()
ggsave(filename = "results/04_scree.png",device = "png")

proteome_pca_aug <- proteome_pca %>%
  augment(proteome_data_aug) %>%
  mutate(PAM50_mRNA = joined_data_aug$PAM50_mRNA)

proteome_pca_aug %>% 
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, label = patient_ID, color = PAM50_mRNA)) +
  geom_text() +
  theme(legend.position = "bottom")
ggsave(filename = "results/04_PCA.png",device = "png")
