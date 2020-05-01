# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("broom")
if (!requireNamespace("gridExtra", quietly = TRUE))
  install.packages("gridExtra")
library(gridExtra)
# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
clincal_data_aug <- read_csv(file = "data/01_clinical_data_clean.csv")
PAM50_clean <- read_csv(file = "data/01_PAM50_clean.csv")
joined_data_aug <- read_csv(file = "data/02_joined_data_aug.csv")

# Catrines plots
# ------------------------------------------------------------------------------
ggplot(data = clincal_data_aug) +
  geom_histogram(mapping = aes(x = Age_at_Initial_Pathologic_Diagnosis), binwidth = 5)
ggsave(filename = "results/03_age_distribution.png",device = "png")

ggplot(data = joined_data_aug, mapping = aes(x = NP_057427, y = NP_002408)) + 
  geom_point()
ggsave(filename = "results/03_random_gene_correlation.png",device = "png")

# View class distribution
joined_data_aug %>% count(PAM50_mRNA) %>% print
joined_data_aug %>% count(Tumor) %>% print


######################## TAG: VLAD #################
# Preliminary fast hand plots
# ------------------------------------------------------------------------------
clincal_data_aug %>% 
  ggplot(mapping = aes(Gender, fill = Tumor)) +
  geom_bar() +
  labs(y = "Count",
       title = "Tumor subtype: occurance based on gender") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5,
                                  size = 20,)
        ) +
  theme_bw(base_family = "Times", 
           base_size = 15)
ggsave("results/gender_vs_tumortype.png", device = "png")    

clincal_data_aug %>% 
  ggplot(mapping = aes(y = Tumor, x =Metastasis_Coded, colour = Tumor)) +
  geom_jitter(width = 0.15, size=4) +
  theme_bw(base_family = "Times", 
           base_size = 15) +
  labs(x = "Metastasis",
       y = "Tumor Type",
       title = "Tumor metastasis based on the subtype of tumor")+
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5,
                                  size = 20,)
  )
ggsave("results/metastatis_vs_tumortype.png", device = "png")

clincal_data_aug %>% 
  ggplot(mapping = aes(y = Tumor, x =methylation_Clusters, colour = Tumor)) +
  geom_jitter(width = 0.1, size = 4) +
  theme_bw(base_family = "Times", 
           base_size = 15) +
  labs(y = "Count",
       x = "Methylation cluster",
       title = "Tumor subtype: occurance based on gender") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5,
                                  size = 20,)
  ) 

ggsave("results/methylCluster_vs_tumorType.png", device = "png")

clincal_data_aug %>% 
  ggplot(mapping = aes(y = Tumor, x =Age_at_Initial_Pathologic_Diagnosis, colour = PAM50_mRNA)) +
  geom_jitter(width = 0.1) +
  labs(y = "Count",
       title = "Tumor subtype: occurance based on gender") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5,
                                  size = 20,)
  ) +
  theme_bw(base_family = "Times", 
           base_size = 15)
ggsave("results/age_vs_cellTyp.png", device = "png")


# HEATMAP MATRIX prep
# Subset the expression data matrix & assign the "row names"
# Overwrite any remaining NAs
joined_data_aug_mat <- joined_data_aug
joined_data_aug_mat[is.na(joined_data_aug_mat)] <- 0

mtx <- joined_data_aug_mat %>%
                filter(PAM50_mRNA != "Control") %>%
                column_to_rownames( ., "patient_ID") %>%
                select(starts_with("NP_"))
                
  
# Create a distance matrix using "euclidean" method
mtx <- mtx %>%
        dist(., method ="euclidean")

# Create an annotation file: assign row names and their respective tissue types
annotation <- joined_data_aug_mat %>% 
                filter(PAM50_mRNA != "Control") %>%
                column_to_rownames(., "patient_ID") %>%
                rename('PAM50 Profile' = PAM50_mRNA) %>%
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


# Distribution of expression for all samples
# ------------------------------------------------------------------------------
joined_data_aug %>%   select(patient_ID,
                             starts_with("NP_"), 
                             PAM50_mRNA
                             ) %>%
                      pivot_longer(cols = starts_with('NP_')) %>%
                        ggplot(aes(x=value, fill=patient_ID)) + geom_density(alpha=0.5) + 
                        geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
                        ggtitle("Density") +
                        theme_bw(base_family = "Times", base_size = 10) +
                        theme(legend.position = "none")
ggsave(filename = "results/density_sample_expression.png",device = "png")



# Try subsetting by iteration using the "tidy" way
df_PAM50_split <- joined_data_aug %>%
                    group_split(PAM50_mRNA)

# Function: Use the data splits based on PAM50 subtype to generate plots iteratively
plotting_PAM50_density <- function (data) {
                        subset <- data %>% select(PAM50_mRNA) %>% unique(.)
                        title <- paste0(subset," tissue: Density")
                        file_prefex <- "results/03_EDA_tissue_"
                        file_suffix <- "_density.png"
                        data %>%   select(patient_ID,
                                          starts_with("NP_"),
                                          PAM50_mRNA
                                          ) %>%
                                   pivot_longer(cols = starts_with('NP_')) %>%
                                   ggplot(aes( x=value, 
                                               fill=patient_ID
                                            )) + 
                                   geom_density(alpha=0.5) + 
                                   geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
                                   ggtitle(title) +
                                   theme_bw(base_family = "Times", 
                                            base_size = 10) +
                                   theme(legend.position = "none") +
                                   labs(x = "Log2 Expression",
                                        y = "Density")
                        ggsave( filename = paste0(file_prefex, subset,file_suffix),
                                device = "png")
}

# Mappings between the subsets and the plot

plots_density <- map(df_PAM50_split, ~plotting_PAM50_density(.x))


# Combined densities per Tissue group
joined_data_aug %>%   select( patient_ID,
                              starts_with("NP_"),
                              PAM50_mRNA
                             ) %>%
                              pivot_longer(cols = starts_with('NP_')) %>%
                                ggplot(aes(x=value, fill=PAM50_mRNA)) + 
                                geom_density(alpha=0.5,) + 
                                geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
                                ggtitle("") +
                                theme_bw(base_family = "Times", base_size = 10)

ggsave(filename = "results/03_EDA_density_tissue_perGroup.png",
       device = "png")


#### Boxplot

joined_data_aug %>%   filter(PAM50_mRNA != "Control") %>%
                      select(patient_ID,
                             starts_with("NP_"),
                             PAM50_mRNA
                             ) %>%
                        pivot_longer(cols = starts_with('NP_')) %>%
                        ggplot(aes(y = patient_ID, 
                                   x = value,
                                   fill = PAM50_mRNA
                                   )
                               ) + 
                        geom_boxplot(alpha=0.5,
                                     varwidth = TRUE,
                                     outlier.shape = NA
                                    ) + 
                        facet_grid(rows = "PAM50_mRNA",
                                   scales = "free",
                                   space = "free_y") +
                        labs(x = "Log2 Expression",
                             y = "Density", 
                             title = "Boxplot: Expression profiles between PAM50 identified groups" ) + 
                        geom_vline(xintercept = c(-1.5, 0.5), linetype="dashed") + # based on Control sample profiles
                        theme_bw(base_family = "Times", 
                                 base_size = 10) +
                        theme(legend.position = "none",
                              axis.text.y = element_text(angle = 30,
                                                           size = 7,
                                                           face = 'italic'
                                                         ),
                              strip.text.y = element_text(size = 12)
                              )

ggsave(filename = "results/03_EDA_boxplot.png",
       device = "png",
       height = 10,
       )


# Boxplot: if I try combining them individually into a grid

# Data split: without CONTROL samples
df_PAM50_split <- joined_data_aug %>%
                  filter(PAM50_mRNA != "Control") %>%
                  group_split(PAM50_mRNA)

plotting_boxplot <- function(data, subset_term, color) {
                    data %>%   
                    select(patient_ID,
                           starts_with("NP_"),
                           PAM50_mRNA
                           ) %>%
                    subset(PAM50_mRNA == subset_term) %>%
                    pivot_longer(cols = starts_with('NP_')) %>%
                    ggplot(aes( y = patient_ID, 
                                x = value,
                              )
                           ) + 
                    geom_boxplot(alpha=0.5,
                                 varwidth = TRUE,
                                 outlier.shape = NA,
                                 fill = color,
                                 ) + 
                    labs(x = "Log2 Expression levels",
                         y = "Patients",
                         title = subset_term) +
                    xlim(-10,10) +
                    geom_vline(xintercept = c(-1.5, 0.5), linetype="dashed") + # based on Control sample profiles
                    theme_bw(base_family = "Times", 
                             base_size = 10) +
                    theme(legend.position = "none", 
                          axis.text.y = element_blank(),
                          plot.title = element_text(hjust = 0.5,
                                                    size = 20,)
                          )
                  }


# Mappings between the subsets and the plot: create separate plots , each colored individually to be combined later
p1_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "Basal-like", color = "red")
p2_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "HER2-enriched", color = "green")
p3_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "Luminal A", color = "turquoise3")
p4_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "Luminal B", color = "purple")

combo_plot <- grid.arrange(p1_boxplot, p3_boxplot, p2_boxplot, p4_boxplot, ncol =2)

ggsave(plot = combo_plot, filename = "results/03_EDA_boxplot_combined.png",
       device = "png",
       height = 5,
)
######################## TAG: VLAD #################
