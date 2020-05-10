# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(broom)
library(patchwork)
library(RColorBrewer)


# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
joined_data_aug <- read_csv(file = "data/02_joined_data_PAM50_aug.csv")

# Wrangle data
# ------------------------------------------------------------------------------
# Remove healthy control samples
joined_data_aug <- joined_data_aug %>% 
  filter(Class != "Control")

proteome_data <- joined_data_aug %>%
  # select only proteome-count columns
  select(starts_with("NP"))            


# Custom colors
# ---------------------------------------------------------------------------
n_levels <- joined_data_aug %>% 
  select(Class) %>% 
  n_distinct

custom_colors <- brewer.pal(n_levels, "Set1")

names(custom_colors) <- joined_data_aug %>% 
  count(Class) %>% 
  select(Class) %>% 
  pull()

  

# PCA 
# ---------------------------------------------------------------------------
## Compute PCA
pca <- proteome_data %>% 
  prcomp(center = TRUE, scale = TRUE) 

## Scree plot
pca %>%
  tidy("pcs") %>% 
  ggplot(aes(x = PC, y = percent)) +
  geom_col() +
  labs(title = "Scree plot - PCA proteome data", 
       x = "Principal components",
       y = "variance explained (%)") +
  theme(base_size = 18,
        plot.title = element_text(hjust = 1.5, size = 25)) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent)

ggsave(filename = "results/04_scree_plot.png", device = "png",
       height = 6)


## Augment and add y class
proteome_pca_aug <- pca %>%
  augment(proteome_data) %>%
  mutate(Class = factor(joined_data_aug$Class, levels = c("Basal", "HER2", "LumA", "LumB")))

## Get PC percent
PC1_perc <- pca %>% 
  tidy("pcs") %>% 
  filter(PC==1) %>% 
  pull(percent) 

PC2_perc <- pca %>% 
  tidy("pcs") %>% 
  filter(PC==2) %>% 
  pull(percent) 

## Scatter proteome data - PC1/PC2
proteome_pca_aug %>% 
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = Class)) +
  geom_point() +
  labs(title = "PCA plot of proteome data", 
       x = str_c("PC1 (", round(PC1_perc * 100, 2), "%)" ),
       y = str_c("PC2 (", round(PC2_perc * 100, 2), "%)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = custom_colors)
ggsave(filename = "results/04_PCA.png", device = "png",
       height = 6)



# K-means clustering
# ------------------------------------------------------------------------------
k = n_levels # 4 levels

### Clustering on original data
set.seed(12)
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


# Which clustering technique performs better
# ------------------------------------------------------------------------------

accuracy <- proteome_pca_cluster_aug %>%
  
  select(Class, cluster_original, cluster_pca) %>%
  
  mutate(cluster_original = case_when(cluster_original == 1 ~ "HER2",
                                      cluster_original == 2 ~ "LumA",
                                      cluster_original == 3 ~ "LumB",
                                      cluster_original == 4 ~ "Basal"),
         cluster_pca = case_when(cluster_pca == 1 ~ "LumA",
                                 cluster_pca == 2 ~ "Basal",
                                 cluster_pca == 3 ~ "HER2",
                                 cluster_pca == 4 ~ "LumB"),
         
         cluster_original_correct = case_when(Class == cluster_original ~ 1,
                                              Class != cluster_original ~ 0),
         cluster_pca_correct = case_when(Class == cluster_pca ~ 1,
                                         Class != cluster_pca ~ 0)) %>% 
  
  summarise(score_original = mean(cluster_original_correct) * 100,
            score_pca = mean(cluster_pca_correct) * 100)



# Visualization of clusters on PCs
# ------------------------------------------------------------------------------
### Original classes
plot1 <- proteome_pca_cluster_aug %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = Class)) +
  geom_point() +
  labs(title = "Original data",
       x = 'PC1',
       y = 'PC2',
       colour = "real class") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title.align = 0.5,
        legend.title = element_text (size = 10),
        legend.text = element_text (size = 8),
        legend.key.size = unit (0.1, "cm")) +
  scale_fill_manual(values = custom_colors) +
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
       subtitle = paste0("accuracy = ", round(accuracy[[1]], 1), "%"),
       x = 'PC1',
       y = 'PC2',
       colour = "clusters") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        legend.position = "bottom",
        legend.title.align = 0.5,
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.1,"cm")) +
  scale_fill_manual(values = custom_colors) +
  guides(colour = guide_legend(title.position="top")
  )



### Clusters on dimensionality-reduced data (first 2 PCs)
plot3 <- proteome_pca_cluster_aug %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = cluster_pca)) +
  geom_point() +
  labs(title = "Clusters on\nPCA data",
       subtitle = paste0("accuracy = ", round(accuracy[[2]], 1), "%"),
       x = 'PC1',
       y = 'PC2',
       colour = "clusters") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        legend.position = "bottom",
        legend.title.align = 0.5,
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.1,"cm")) +
  scale_fill_manual(values = custom_colors) +
  guides(colour = guide_legend(title.position="top"))

(plot1 + plot2 + plot3) 
ggsave(filename = "results/04_PCA_kmeans.png", device = "png",
       height = 5)


