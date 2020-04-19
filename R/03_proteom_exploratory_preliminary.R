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
df <- read_csv(file = "data/02_joined_data_clean.csv")

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



mtx %>%
  pheatmap(.,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           annotation = annotation,angle_col = 45)


