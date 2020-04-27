######## TO BE MOVED TO 03_EDA ######## 

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
# ------------------------------------------------------------------------------

# Overwrite missing values
# ------------------------------------------------------------------------------
# Overwrite any remaining NAs
df[is.na(df)] <- 0
# ------------------------------------------------------------------------------



# HEATMAP: consists of: DISTANCE matrix and ANNOTATION
# ------------------------------------------------------------------------------
# Subset the expression data matrix & assign the "row names"
proteom_mtx <- df %>% 
  column_to_rownames(., "patient_ID") %>%
  select(starts_with("NP_"),
         PAM50_mRNA)  




# Create a distance matrix using "euclidean" method
mtx <- proteom_mtx %>%
       select(-PAM50_mRNA) %>%
       dist(.,method ="euclidean")

# Create an annotation file: assign row names and their respective tissue types
annotation <- df %>% 
              column_to_rownames(., "patient_ID") %>%
              rename('PAM50 Profile' = PAM50_mRNA) %>%
              select('PAM50 Profile')

mtx %>%  pheatmap(.,
                  clustering_distance_rows = "correlation",
                  clustering_distance_cols = "correlation",
                  annotation = annotation,
                  angle_col = 45,
                  scale = "none",
                  cellwidth = 5,
                  cellheight = 5,
                  #show_rownames = TRUE, 
                  main = "Heatmap using Euclidean distance:\nPAM50 genes' expressions across samples",
                  filename = "results/03_exploratory_heatmap_tissueAnnotation1.png",
                  ) 
dev.off()


# Distribution of expression for all samples
# ------------------------------------------------------------------------------
df %>%   select(patient_ID,
                starts_with("NP_"),
                PAM50_mRNA) %>%
  pivot_longer(cols = starts_with('NP_')) %>%
  ggplot(aes(x=value, fill=patient_ID)) + geom_density(alpha=0.5) + 
  geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
  ggtitle("Density") +
  theme_bw(base_family = "Times", base_size = 10) +
  theme(legend.position = "none")
ggsave(filename = "results/density_sample_expression.png",device = "png")



# Try subsetting by iteration using the "tidy" way
df_PAM50_split <- df %>%
        group_split(PAM50_mRNA)

# Function for the plot that was used previously
plot_PAM50_density <- function (data) {
  title <- data %>% select(PAM50_mRNA) %>% unique(.)
  data %>%   select(patient_ID,
                starts_with("NP_"),
                PAM50_mRNA) %>%
  pivot_longer(cols = starts_with('NP_')) %>%
  ggplot(aes(x=value, fill=patient_ID)) + geom_density(alpha=0.5) + 
  geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
  ggtitle(paste0(title,": Density")) +
  theme_bw(base_family = "Times", base_size = 10) +
  theme(legend.position = "none")
  ggsave(filename = paste0("results/03_exploratory_tissue_",title,"_density.png"),device = "png")
}

# Mappings between the subsets and the plot

plots <- 
  map(df_PAM50_split, ~plot_PAM50_density(.x))
# # All subtypes
# PAM50_subtypes <- df %>% select(PAM50_mRNA) %>% unique(.) 
# df %>%   select(patient_ID,
#                   starts_with("NP_"),
#                   PAM50_mRNA) %>%
#            filter(.,PAM50_mRNA == as.character(PAM50_subtypes[1,1])) %>%
#     pivot_longer(cols = starts_with('NP_')) %>%
#     ggplot(aes(x=value, fill=patient_ID)) + geom_density(alpha=0.5) + 
#     geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
#     ggtitle(paste("PAM50: ",as.character(PAM50_subtypes[1,1])," expression density plot.")) +
#     theme_bw(base_family = "Times", base_size = 10) +
#     theme(legend.position = "center")
# ggsave(filename = paste("03_exploratory_density_",as.character(PAM50_subtypes[1,1]),".png"),device = "png")
#   
# df %>%   select(patient_ID,
#                   starts_with("NP_"),
#                   PAM50_mRNA) %>%
#     filter(.,PAM50_mRNA ==as.character(PAM50_subtypes[2,1])) %>%
#     pivot_longer(cols = starts_with('NP_')) %>%
#     ggplot(aes(x=value, fill=patient_ID)) + geom_density(alpha=0.5) + 
#     geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
#     ggtitle(paste("PAM50: ",as.character(PAM50_subtypes[2,1])," expression density plot.")) +
#     theme_bw(base_family = "Times", base_size = 10) +
#     theme(legend.position = "center")
# ggsave(filename = paste("03_exploratory_density_",as.character(PAM50_subtypes[1,1]),".png"),device = "png")
#   
# df %>%   select(patient_ID,
#                   starts_with("NP_"),
#                   PAM50_mRNA) %>%
#     filter(.,PAM50_mRNA == as.character(PAM50_subtypes[3,1])) %>%
#     pivot_longer(cols = starts_with('NP_')) %>%
#     ggplot(aes(x=value, fill=patient_ID)) + geom_density(alpha=0.5) + 
#     geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
#     ggtitle(paste("PAM50: ",as.character(PAM50_subtypes[3,1])," expression density plot.")) +
#     theme_bw(base_family = "Times", base_size = 10) +
#     theme(legend.position = "center")
#   ggsave(filename = paste("03_exploratory_density_",as.character(PAM50_subtypes[3,1]),".png"),device = "png")
#   
#   df %>%   select(patient_ID,
#                   starts_with("NP_"),
#                   PAM50_mRNA) %>%
#     filter(.,PAM50_mRNA == as.character(PAM50_subtypes[4,1])) %>%
#     pivot_longer(cols = starts_with('NP_')) %>%
#     ggplot(aes(x=value, fill=patient_ID)) + geom_density(alpha=0.5) + 
#     geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
#     ggtitle(paste("PAM50: ",as.character(PAM50_subtypes[4,1])," expression density plot.")) +
#     theme_bw(base_family = "Times", base_size = 10) +
#     theme(legend.position = "center")
#   ggsave(filename = paste("03_exploratory_density_",as.character(PAM50_subtypes[4,1]),".png"),device = "png")
#   

# ------------------------------------------------------------------------------
df %>%   select(patient_ID,
                starts_with("NP_"),
                PAM50_mRNA) %>%
  pivot_longer(cols = starts_with('NP_')) %>%
  ggplot(aes(x=value, fill=PAM50_mRNA)) + geom_density(alpha=0.5,) + 
  geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
  ggtitle("") +
  theme_bw(base_family = "Times", base_size = 10) +
  theme(legend.position = "none")
ggsave(filename = "results/03_exloratory_density_tissue_expression.png",device = "png")
