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
clincal_data_aug <- read_csv(file = "data/02_clincal_data_clean.csv")
PAM50_aug <- read_csv(file = "data/02_PAM50_clean.csv")
proteome_data_aug <- read_csv(file = "data/02_proteome_data_clean.csv")
joined_data_aug <- read_csv(file = "data/02_joined_data_clean.csv")

# Wrangle data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% ...

# Model data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
ggplot(data = clincal_data_aug) +
  geom_histogram(mapping = aes(x = Age_at_Initial_Pathologic_Diagnosis), binwidth = 5)

ggplot(data = clincal_data_aug, mapping = aes(x = Age_at_Initial_Pathologic_Diagnosis, y = Days_to_Date_of_Last_Contact)) + 
  geom_point()


ggplot(data = proteome_data_aug, mapping = aes(x = NP_057427, y = NP_002408)) + 
  geom_point()

# View class distribution
joined_data_aug %>% count(PAM50_mRNA) %>% print
joined_data_aug %>% count(Tumor) %>% print

# PCA ---------------------------------------------------------------------------
proteome_pca <- proteome_data_aug %>%
  select(-patient_ID) %>%
  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% # Take median value of NA values
  prcomp(center = TRUE, scale = TRUE)

proteome_pca %>%
  tidy("pcs") %>% 
  ggplot(aes(x = PC, y = percent)) +
  geom_col() +
  theme_bw()

proteome_pca_aug <- proteome_pca %>%
  augment(proteome_data_aug)

proteome_pca_aug %>% 
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, label = patient_ID)) +
  geom_text() +
  theme(legend.position = "bottom")


#%%%%
proteome_pca <- proteome_data_aug %>%
  select(starts_with("NP"))  %>%
  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
  prcomp(center = TRUE, scale = TRUE) 
  

proteome_pca_aug <- proteome_pca %>%
  augment(proteome_data_aug) %>%
  mutate(PAM50_mRNA = joined_data_aug$Tumor)

proteome_pca_aug %>% 
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, label = patient_ID, color = PAM50_mRNA)) +
  geom_text() +
  theme(legend.position = "bottom")
  
# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data_clean_aug,
          path = "data/04_my_data_clean_aug_anl.tsv")
ggsave(...)