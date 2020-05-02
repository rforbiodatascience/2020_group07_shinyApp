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
joined_data_aug <- read.csv(file = "data/02_joined_data_aug.csv")

# Wrangle data
# ------------------------------------------------------------------------------
proteome_data <- joined_data_aug %>% 
  select(starts_with("NP"))

# PCA 
# ---------------------------------------------------------------------------
### Compute PCA
proteome_pca <- proteome_data %>% 
  prcomp(center = TRUE, scale = TRUE) 

### Scree plot
proteome_pca %>%
  tidy("pcs") %>% 
  ggplot(aes(x = PC, y = percent)) +
  geom_col() +
  labs(title = "Scree plot - PCA proteome data", 
       x = "Principal components",
       y = "variance explained (%)") +
  scale_y_continuous(labels = scales::percent)
ggsave(filename = "results/04_scree_plot.png",device = "png")

### Augment and add y class (PAM50_mRNA clusters)
proteome_pca_aug <- proteome_pca %>%
  augment(proteome_data) %>%
  mutate(PAM50_mRNA = joined_data_aug$PAM50_mRNA)


### Scatter proteome data - PC1/PC2
proteome_pca_aug %>% 
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = PAM50_mRNA)) +
  geom_point() +
  labs(title = "PCA - proteome data", 
       x = "PC1",
       y = "PC2") +
ggsave(filename = "results/04_PCA.png",device = "png")


# Clustering
# ------------------------------------------------------------------------------


