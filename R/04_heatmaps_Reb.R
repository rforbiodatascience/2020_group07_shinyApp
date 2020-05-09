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
proteome_raw_data <- read_csv(file = "data/_raw/77_cancer_proteomes_CPTAC_itraq.csv")
clinical_raw_data <- read_csv(file = "data/_raw/77_cancer_proteomes_CPTAC_itraq.csv")


# Wrangle data
# ------------------------------------------------------------------------------
### Extract class information from clinical data
class_data <- clinical_raw_data %>% 
  # unify patient IDs
  rename(patient_ID = "Complete TCGA ID") %>% 
  mutate(patient_ID = str_sub(patient_ID,
                              start = 6,
                              end = -1)) %>% 
  select(patient_ID, "PAM50 mRNA")

  
### Identify patient_ID replicates 
replicates <- colnames(proteome_raw_data) %>% 
  # Simplify ID name
  str_replace_all(., '\\..*', '') %>% 
  # Find replicates (true/false)
  duplicated() %>% 
  # Extract replicate column names (excluding first apperance)
  colnames(proteome_raw_data)[.]


### Create long data for heatmap
gene_set <- proteome_raw_data %>% 
  # filter for genes of interest
  filter(gene_symbol == "TP53" |
           gene_symbol == "PIK3CA" | 
           gene_symbol == "GATA3" |
           gene_symbol == "ESR1" |
           gene_symbol == "PGR" |
           gene_symbol == "ERBB2") %>% 
 
  # keep only relevant columns for heatmap
  select(-replicates, gene_symbol, ends_with("TCGA"), ends_with("CPTAC")) %>% 
  # create long data for heatmap
  pivot_longer(cols = "AO-A12D.01TCGA":"c4155b-C.CPTAC", 
               names_to = "patient_ID",
               values_to = "ITRAQ_log2_ratio") %>% 
  # unify patient IDs
  mutate(patient_ID = str_sub(patient_ID,
                              start = 1,
                              end = 7)) %>% 
  # add Class information
  left_join(class_data, by = "patient_ID") %>% 
  # Rename Class names
  rename(Class = "PAM50 mRNA" ) %>% 
  mutate(Class = replace_na(Class, "Control"),
         Class = str_replace(Class,
                             pattern = "HER2-enriched",
                             replacement = "HER2"),
         Class = str_replace(Class,
                             pattern = "Luminal A",
                             replacement = "LumA"),
         Class = str_replace(Class,
                             pattern = "Luminal B",
                             replacement = "LumB"),
         Class = str_replace(Class,
                             pattern = "Basal-like",
                             replacement = "Basal")) %>%
  # Transform variables to factor for the plot
  mutate(Class = factor(Class, levels = c("Basal", "HER2", "LumA", "LumB", "Control"))) %>% 
  mutate(gene_symbol = factor(gene_symbol, levels = c("ESR1", "PGR", "ERBB2", "GATA3", "PIK3CA", "TP53")))



# Create heatmap
# ------------------------------------------------------------------------------
heatmap <- gene_set %>% 
  ggplot(mapping = aes(gene_symbol, patient_ID, fill = ITRAQ_log2_ratio, colour = "")) +
  geom_tile() +
  facet_grid(Class ~ ., scales = "free") +
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       name = "ITRAQ log2 ratio",
                       na.value = "black",
                       space = "Lab",
                       guide = "colourbar",
                       aesthetics = "fill") +
  scale_colour_manual(values = NA) +
  guides(colour = guide_legend("No data", override.aes = list(colour = "black"))) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        panel.spacing.y = unit(0.1, "cm")) +
  labs(title = "Frequently mutated genes (TP53, GATA3, PIK3CA) in breast cancer + biomarker genes across classes",
       x = NULL,
       y = NULL) 
  
ggplotly(heatmap)
ggsave(filename = "results/04_heatmap_specific_genes.png", device = "png")

