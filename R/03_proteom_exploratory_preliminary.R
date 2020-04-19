# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("pheatmap")

# Define functions (DELETE IF NOT USED)
# ------------------------------------------------------------------------------


# Load data
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
df <- read_csv(file = "data/02_joined_data_clean.csv")

# ------------------------------------------------------------------------------
proteom_mtx <- df %>% 
  column_to_rownames(., "patient_ID") %>%
  select(starts_with("NP_"),
         PAM50_mRNA)  

proteom_mtx[is.na(proteom_mtx)] <- 0  
mtx <- proteom_mtx %>%
  select(-PAM50_mRNA) %>%
  dist(proteom_mtx,method ="euclidean")

annotation <- df %>% column_to_rownames(., "patient_ID") %>%
  select(PAM50_mRNA) 


png(filename = "results/heatmap_tissueAnnotation.png")
mtx %>%  pheatmap(.,
                  clustering_distance_rows = "correlation",
                  clustering_distance_cols = "correlation",
                  annotation = annotation,angle_col = 45)
dev.off()

# ------------------------------------------------------------------------------
df %>%   select(patient_ID,
                starts_with("NP_"),
                PAM50_mRNA) %>%
  pivot_longer(cols = starts_with('NP_')) %>%
  ggplot(aes(x=value, fill=patient_ID)) + geom_density(alpha=0.5, bins=100) + 
  geom_vline(xintercept = c(-0.58, 0.58), linetype="dashed") + 
  ggtitle("") +
  theme_bw(base_family = "Times", base_size = 10) +
  theme(legend.position = "none")
ggsave(filename = "results/density_sample_expression.png",device = "png")


# ------------------------------------------------------------------------------
df %>%   select(patient_ID,
                starts_with("NP_"),
                PAM50_mRNA) %>%
  pivot_longer(cols = starts_with('NP_')) %>%
  ggplot(aes(x=value, fill=PAM50_mRNA)) + geom_density(alpha=0.5, bins=100) + 
  geom_vline(xintercept = c(-0.58, 0.58), linetype="dashed") + 
  ggtitle("") +
  theme_bw(base_family = "Times", base_size = 10) +
  theme(legend.position = "none")
ggsave(filename = "results/density_tissue_expression.png",device = "png")
