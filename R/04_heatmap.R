# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(plotly)


# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data
# ------------------------------------------------------------------------------
data <- read_csv(file = "data/02_joined_data_full_aug.csv")

# Wrangle data
# ------------------------------------------------------------------------------
### Create long data for heatmap
gene_set <- data %>% 
  # keep only relevant columns for heatmap
  select(patient_ID, Class, NP_000537, NP_006209, NP_001002295, NP_000116, NP_000917, NP_004439) %>% 
  # give gene names to RefSeq numbers
  rename("TP53" = "NP_000537",
         "PIK3CA" = "NP_006209",
         "GATA3" = "NP_001002295",
         "ESR1" = "NP_000116",
         "PGR" = "NP_000917",
         "ERBB2" = "NP_004439") %>% 
  # create long data for heatmap
  pivot_longer(cols = c("TP53", "PIK3CA", "GATA3", "ESR1", "PGR", "ERBB2"), 
               names_to = "RefSeq",
               values_to = "ITRAQ_log2_ratio") %>% 
  # set factors and levels fot the plot
  mutate(Class = factor(Class, levels = c("Basal", "HER2", "LumA", "LumB", "Control"))) %>% 
  mutate(RefSeq = factor(RefSeq, levels = c("TP53", "PIK3CA", "GATA3", "ESR1", "PGR", "ERBB2")))

# Create heatmap
# ------------------------------------------------------------------------------
heatmap <- gene_set %>% 
  ggplot(mapping = aes(RefSeq, patient_ID, fill = ITRAQ_log2_ratio)) +
  geom_tile() +
  facet_grid(Class ~ ., scales = "free") +
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       name = "ITRAQ log2 ratio") +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 15),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        panel.spacing.y = unit(0.1, "cm")) +
  labs(title = "Expression of common biomarkers in breast cancer across tumor classes",
       x = NULL,
       y = NULL)
  
ggsave(filename = "results/04_heatmap_specific_genes.png", device = "png")

