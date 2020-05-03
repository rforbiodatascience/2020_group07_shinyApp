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

# Combining x4 Exploratory analysis plots into a single canvas
# ------------------------------------------------------------------------------
# Plot 1/4
p11 <- clincal_data_aug %>% 
  ggplot(mapping = aes(Gender, 
                       fill = Tumor
                       )
         ) +
  geom_bar() +
  theme_bw(base_family = "Times", 
           base_size = 12
           ) +
  labs(y = "Count",
       title = "Tumor subtype: occurance based on gender"
       ) +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14
                                  )
        ) 
# Save the plot then remove the legend for the canvas
ggsave(plot = p11,"results/gender_vs_tumortype.png", device = "png")    
p11 <- p11 + theme(legend.position = "none")

# Plot 2/4
p22 <- clincal_data_aug %>% 
  ggplot(mapping = aes(y = Tumor, 
                       x = Metastasis_Coded,
                       colour = Tumor
                       )
         ) +
  geom_jitter(width = 0.15, 
              size=2
              ) +
  theme_bw(base_family = "Times", 
           base_size = 12
           ) +
  labs(x = "Metastasis",
       y = "Tumor Type",
       title = "Metastasis of different tumor subtypes"
       ) +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14
                                  )
        )

# Save then remove the legend
ggsave(plot = p22, "results/metastatis_vs_tumortype.png", device = "png")
p22 <- p22 + theme(legend.position = "none")

# Plot 3/4
p33 <- clincal_data_aug %>% 
  ggplot(mapping = aes(y = Tumor, 
                       x = methylation_Clusters,
                       colour = Tumor
                       )
         ) +
  geom_jitter(width = 0.1, 
              size = 2
              ) +
  theme_bw(base_family = "Times", 
           base_size = 12
           ) +
  labs(y = "Tumor Type",
       x = "Methylation cluster",
       title = "Methylation clusters of different tumor subtypes"
       ) +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14
                                  )
        ) 

# Save then remove the legend
ggsave(plot = p33, "results/methylCluster_vs_tumorType.png", device = "png")
p33 <- p33 + theme(legend.position = "none")

# Plot 4/4
p44 <- clincal_data_aug %>% 
        ggplot(mapping = aes(y = Tumor, 
                             x =Age_at_Initial_Pathologic_Diagnosis,
                             colour = Tumor
        )
        ) +
        geom_violin(scale = "area", 
                    trim = TRUE,
                    mapping = aes(fill=Tumor)
        ) +
        #geom_dotplot(binaxis='x', dotsize=1) +
        stat_summary(fun=median, geom="point", 
                     size=2, color="black") +
        theme_bw(base_family = "Times", 
                 base_size = 12
        ) +
        labs(y = "Tumor Type",
             x = " Age at initial diagnosis",
             title = "Tumor subtype based on person's age"
        ) +
        theme(plot.title = element_text(hjust = 0.5, 
                                        size = 14
        )
        )
# Save then adjust the legend,  save it as a variable and remove from the plot itself
ggsave(plot = p44, "results/age_vs_tumorType.png", device = "png")
legend_temp <- p44 + theme(legend.text = element_text(size = 20),
            legend.title = element_text(size = 20),
            legend.key.size = unit(20,"point"),
            legend.position = "bottom",
            #legend.box.background = element_rect(colour = "black")
            ) 

# Call the helper function for legend extraction
shared_legend <- get_legend(legend_temp)

# Remove the legend from the remaing plot
p44 <- p44 + theme(legend.position = "none")

# Combine the 4 plots and the shared legend
plot_EDA1 <- grid.arrange(p11, p22, p33, p44, shared_legend,
             ncol= 2,
             widths = c(2.7, 2.7),
             nrow = 3,
             heights = c(2.7, 2.7, 0.5),
             layout_matrix = rbind(c(1,2), c(3,4), c(5,5)) # to bind the legend as a shared 3rd row c(el.5, el.5)
             )
              
# Save the combined plot
ggsave(plot=plot_EDA1, "results/03_EDA_grid4_combo.png", device = "png" )




# HEATMAP MATRIX
# ------------------------------------------------------------------------------
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
joined_data_aug %>%   select(patient_ID,
                            starts_with("NP_"),
                            PAM50_mRNA
                            ) %>%
                            pivot_longer(cols = starts_with('NP_')) %>%
                            ggplot(aes(x=value, fill=PAM50_mRNA)) + 
                            geom_density(alpha=0.5,) + 
                            geom_vline(xintercept = c(-1, 1), linetype="dashed") + 
                            ggtitle("") +
                            theme_bw(base_family = "Times", base_size = 10)

ggsave(filename = "results/03_EDA_density_tissue_perGroup.png", device = "png")


# BOXPLOTS: combining 4x plots into one canvas
# ------------------------------------------------------------------------------
# Tidy approacj to orthogonal data processing:

# Data split for each cancer group: without CONTROL samples
df_PAM50_split <- joined_data_aug %>%
                  filter(PAM50_mRNA != "Control") %>%
                  group_split(PAM50_mRNA)

# ggplot function for handling a single instance
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
                    geom_vline(xintercept = c(-1.5, 0.5), linetype="dashed", co) + # based on Control sample profiles
                    theme_bw(base_family = "Times", 
                             base_size = 10) +
                    theme(legend.position = "none", 
                          axis.text.y = element_blank(),
                          plot.title = element_text(hjust = 0.5,
                                                    size = 20,)
                          )
}


# Mappings between the subsets and the plotting function:
# Creating a separate plot while gettign colored individually, to be combined later
p1_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "Basal-like", color = "red")
p2_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "HER2-enriched", color = "green")
p3_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "Luminal A", color = "turquoise3")
p4_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "Luminal B", color = "purple")

plot_EDA2_boxplot_combo <- grid.arrange(p1_boxplot, p3_boxplot, p2_boxplot, p4_boxplot, ncol =2)

ggsave(plot = plot_EDA2_boxplot_combo, filename = "results/03_EDA_boxplot_combined.png",
       device = "png",
       height = 5,
)

# ------------------------------------------------------------------------------
plot_TvsPAM50_boxplot <- joined_data_aug %>%   
                          select(patient_ID,
                                 Tumor,
                                 starts_with("NP_"),
                                 PAM50_mRNA
                          ) %>%
                          filter(PAM50_mRNA != "Control") %>%
                          pivot_longer(cols = starts_with('NP_')) %>%
                          ggplot(aes( y = patient_ID, 
                                      x = value,
                          )
                          ) + 
                          geom_boxplot(alpha=0.5,
                                       varwidth = TRUE,
                                       outlier.shape = NA,
                                       mapping = aes(fill=PAM50_mRNA)
                        
                          ) + 
                          labs(x = "Log2 Expression levels",
                               y = "Patients"
                               ) +
                          xlim(-10,10
                               ) +
                          geom_vline(xintercept = c(-1.5, 0.5), 
                                     linetype="dashed"
                                     ) + # based on Control sample profiles
                          theme_bw(base_family = "Times", 
                                   base_size = 10) +
                          theme(legend.position = "none", 
                                axis.text.y = element_blank(),
                                plot.title = element_text(hjust = 0.5,
                                                          size = 20,)
                                ) +
                          facet_grid(~Tumor
                                     ) + 
                          theme(legend.position = "right")
  
ggsave(plot = plot_TvsPAM50_boxplot, filename = "results/03_EDA_boxplot_T_vs_PAM50_mapping.png",
       device = "png",
)
