# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("broom")
library(patchwork)

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
pca <- proteome_data %>% 
  prcomp(center = TRUE, scale = TRUE) 

### Scree plot
pca %>%
  tidy("pcs") %>% 
  ggplot(aes(x = PC, y = percent)) +
  geom_col() +
  labs(title = "Scree plot - PCA proteome data", 
       x = "Principal components",
       y = "variance explained (%)") +
  theme(plot.title = element_text(hjust = 0.5))
  scale_y_continuous(labels = scales::percent)
ggsave(filename = "results/04_scree_plot.png",device = "png")

### Augment and add y class (PAM50_mRNA clusters)
proteome_pca_aug <- pca %>%
  augment(proteome_data) %>%
  mutate(class = joined_data_aug$PAM50_mRNA)


### Scatter proteome data - PC1/PC2
proteome_pca_aug %>% 
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = class)) +
  geom_point() +
  labs(title = "PCA - proteome data", 
       x = "PC1",
       y = "PC2") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "results/04_PCA.png",device = "png")



# K-means clustering
# ------------------------------------------------------------------------------
k = length(levels(proteome_pca_aug$class)) # 5 levels

### Clustering on original data
cluster_original <- proteome_data %>%
  kmeans(centers = k)

### Augment to PCA data
proteome_pca_cluster_aug <- 
  cluster_original %>%
  broom::augment (proteome_pca_aug) %>% 
  rename(cluster_original = .cluster)


### Clustering on dimensionality-reduced data (2 first PCs)
cluster_pca <- proteome_pca_aug %>%
  select(.fittedPC1, .fittedPC2) %>%
  kmeans(centers = k)

### Augment to PCA/kmeans data
proteome_pca_cluster_aug <- 
  cluster_pca %>% 
  broom::augment (proteome_pca_cluster_aug) %>% 
  rename(cluster_pca = .cluster)



# Visualization of clusters
# ------------------------------------------------------------------------------
### Original classes
plot1 <- proteome_pca_cluster_aug %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = class)) +
  geom_point() +
  labs(title = "Original data",
       x = 'PC1',
       y = 'PC2',
       colour = "class") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title.align = 0.5,
        legend.title = element_text (size = 10),
        legend.text = element_text (size = 8),
        legend.key.size = unit (0.1, "cm")) +
  guides(colour = guide_legend( title.position = "top",
                                nrow = 2,
                                byrow = TRUE))
  



### Clusters on original data
plot2 <- proteome_pca_cluster_aug %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2,
             colour = cluster_original)) +
  geom_point() +
  labs(title = "Clusters on\noriginal data",
       x = 'PC1',
       y = 'PC2',
       colour = "clusters") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        #legend.box = "horizontal",
        legend.title.align=0.5,
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.1,"cm")) +
  guides(colour = guide_legend(title.position="top")
  )



### Clusters on dimensionality-reduced data (first 2 PCs)
plot3 <- proteome_pca_cluster_aug %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = cluster_pca)) +
  geom_point() +
  labs(title = "Clusters on\nPCA data",
       x = 'PC1',
       y = 'PC2',
       colour = "clusters") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        #legend.box = "horizontal",
        legend.title.align=0.5,
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.1,"cm")) +
  guides(colour = guide_legend(title.position="top"))

(plot1 + plot2 + plot3)
ggsave(filename = "results/04_PCA_kmeans.png",
       device = "png")
