# HEATMAP MATRIX
# ------------------------------------------------------------------------------
# Subset the expression data matrix & assign the "row names"
# Overwrite any remaining NAs
joined_data_aug_mat <- joined_data_aug
joined_data_aug_mat[is.na(joined_data_aug_mat)] <- 0

mtx <- joined_data_aug_mat %>%
  filter(Class != "Control") %>%
  column_to_rownames( ., "patient_ID") %>%
  select(starts_with("NP_"))

# Create a distance matrix using "euclidean" method
mtx <- mtx %>%
  dist(., method ="euclidean")

# Create an annotation file: assign row names and their respective tissue types
annotation <- joined_data_aug_mat %>% 
  filter(Class != "Control") %>%
  column_to_rownames(., "patient_ID") %>%
  rename('PAM50 Profile' = Class) %>%
  select('PAM50 Profile')

mtx %>%  pheatmap( .,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   annotation = annotation,
                   angle_col = 45,
                   scale = "none",
                   cellwidth = 5,
                   cellheight = 5,
                   #show_rownames = TRUE, 
                   main = "Heatmap using Euclidean distance:\nPAM50 genes' expressions across samples",
                   filename = "results/03_EDA_heatmap_annotated.png"
) 
dev.off()
