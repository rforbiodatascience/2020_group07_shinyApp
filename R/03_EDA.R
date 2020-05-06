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

# Catrine: if we do the install packages (which is pretty smart),
# We should do it in every script.

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
clincal_data_aug <- read_csv(file = "data/01_clinical_data_clean.csv")
PAM50_clean <- read_csv(file = "data/01_PAM50_clean.csv")
joined_data_aug <- read_csv(file = "data/02_joined_data_aug.csv")

# Rewrite a label of HER2-enriched category
### Catrine: Why ?
joined_data_aug <- joined_data_aug %>% 
                    mutate(PAM50_mRNA = str_replace(PAM50_mRNA,
                                                    pattern = "HER2-enriched",
                                                    replacement = "HER2"))

# Catrines plots
# ------------------------------------------------------------------------------
ggplot(data = clincal_data_aug) +
  geom_histogram(mapping = aes(x = Age_at_Initial_Pathologic_Diagnosis), 
                 binwidth = 5
                 ) +
  labs(title = "Age at diagnosis",
       x = 'Age',
       y = 'Count')
ggsave(filename = "results/03_EDA_age_distribution.png",device = "png")

clincal_data_aug %>%
  ggplot(aes(x = AJCC_Stage, 
             y = Age_at_Initial_Pathologic_Diagnosis, 
             colour = AJCC_Stage)) +
  geom_violin() + 
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5) +
  stat_summary(fun=median, geom="point", 
               size=2, color="black") +
  theme_bw(base_family = "Times", 
           base_size = 12) +
  labs(title = "Age versus cancer severity",
       x = 'AJJC stage',
       y = 'Age at initial diagnosis')
ggsave(filename = "results/03_EDA_age_cancer_severity.png",device = "png")

# ------------------------------------------------------------------------------
# Combining x4 Exploratory analysis plots into a single canvas
# ------------------------------------------------------------------------------
### Catrine: I think it should be seperate plot.
# Plot 1/4
p11 <- clincal_data_aug %>% 
  ggplot(mapping = aes(Gender, 
                       fill = PAM50_mRNA
                       )
         ) +
  geom_bar() +
  theme_bw(base_family = "Times", 
           base_size = 12
           ) +
  labs(y = "Count",
       title = "PAM50 tumor subtype: occurance based on gender"
       ) +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14
                                  )
        ) 
# Save the plot then remove the legend for the canvas
ggsave(plot = p11,"results/03_EDA_gender_vs_tumortype.png", device = "png")    
p11 <- p11 + theme(legend.position = "none")

# Plot 2/4
p22 <- clincal_data_aug %>% 
  ggplot(mapping = aes(y = PAM50_mRNA, 
                       x = Metastasis_Coded,
                       colour = PAM50_mRNA
                       )
         ) +
  geom_jitter(width = 0.15, 
              size=2
              ) +
  theme_bw(base_family = "Times", 
           base_size = 12
           ) +
  labs(x = "Metastasis",
       y = "PAM50 tumor subtype",
       title = "Metastasis of different tumor subtypes"
       ) +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14
                                  )
        )

# Save then remove the legend
ggsave(plot = p22, "results/03_EDA_metastatis_vs_tumortype.png", device = "png")
p22 <- p22 + theme(legend.position = "none")

# Plot 3/4
p33 <- clincal_data_aug %>% 
  ggplot(mapping = aes(y = PAM50_mRNA, 
                       x = methylation_Clusters,
                       colour = PAM50_mRNA
                       )
         ) +
  geom_jitter(width = 0.1, 
              size = 2
              ) +
  theme_bw(base_family = "Times", 
           base_size = 12
           ) +
  labs(y = "PAM50 tumor subtype",
       x = "Methylation cluster",
       title = "Methylation clusters of different tumor subtypes"
       ) +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14
                                  )
        ) 

# Save then remove the legend
ggsave(plot = p33, "results/03_EDA_methylCluster_vs_tumorType.png", device = "png")
p33 <- p33 + theme(legend.position = "none")

# Plot 4/4
p44 <- clincal_data_aug %>% 
        ggplot(mapping = aes(y = PAM50_mRNA, 
                             x = Age_at_Initial_Pathologic_Diagnosis,
                             colour = PAM50_mRNA)) +
        geom_violin(scale = "area", 
                    trim = TRUE,
                    mapping = aes(fill=Tumor)) +
        #geom_dotplot(binaxis='x', dotsize=1) +
        stat_summary(fun=median, geom="point", 
                     size=2, color="black") +
        theme_bw(base_family = "Times", 
                 base_size = 12) +
        labs(y = "PAM50 tumor subtype",
             x = " Age at initial diagnosis",
             title = "Tumor subtype based on person's age"
        ) +
        theme(plot.title = element_text(hjust = 0.5, 
                                        size = 14))

# Save
ggsave(plot = p44, "results/03_EDA_age_vs_tumorType.png", device = "png")

# Make legend adjustments
legend_temp <- p44 + theme(legend.text = element_text(size = 20),
            legend.title = element_text(size = 20),
            legend.key.size = unit(20,"point"),
            legend.position = "bottom") 

# Remove the "legend object"
shared_legend <- get_legend(legend_temp)

# Remove the legend from the original plot
p44 <- p44 + theme(legend.position = "none")

# Combine the 4 plots and the shared legend
plot_EDA1 <- grid.arrange(p11, p22, p33, p44, shared_legend,
             ncol= 2,
             widths = c(2.7, 2.7),
             nrow = 3,
             heights = c(2.7, 2.7, 0.5),
             layout_matrix = rbind(c(1,2), 
                                   c(3,4), 
                                   c(5,5)) # to bind the legend as a shared 3rd row c(el.5, el.5)
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





# BOXPLOTS: combining 4x plots into one canvas
# ------------------------------------------------------------------------------
# Mappings between the subsets and the plotting function:
# Creating a separate plot while gettign colored individually, to be combined later
source(file = "R/99_project_functions.R")
p1_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "Basal-like", color = "red")
p2_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "HER2", color = "green")
p3_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "Luminal A", color = "turquoise3")
p4_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "Luminal B", color = "purple")


# Call the helper function for legend extraction
shared_legend <- get_legend(p4_boxplot)

# Remove the legend from the remaining plot
p1_boxplot <- p1_boxplot + theme(legend.position = "none")
p2_boxplot <- p2_boxplot + theme(legend.position = "none")
p3_boxplot <- p3_boxplot + theme(legend.position = "none")
p4_boxplot <- p4_boxplot + theme(legend.position = "none")

# Combine the 4 plots and the shared legend
plot_EDA2_boxplot_combo <- grid.arrange(p1_boxplot,
                                        p2_boxplot, 
                                        p3_boxplot, 
                                        p4_boxplot,
                                        shared_legend,
                                        ncol= 2,
                                        widths = c(2.7, 2.7),
                                        nrow = 3,
                                        heights = c(2.7, 2.7, 0.5),
                                        layout_matrix = rbind(c(1,2), 
                                                              c(3,4), 
                                                              c(5,5)) # to bind the legend as a shared 3rd row c(el.5, el.5)
                                        )


ggsave(plot = plot_EDA2_boxplot_combo, filename = "results/03_EDA_boxplot_combined.png",
       device = "png",
       height = 5,
)

### Catrine: the black line is the control samples... median +/- something?
### Why is some samples outside the range?

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
                          geom_boxplot(alpha= 0.5,
                                       varwidth = TRUE,
                                       outlier.shape = NA,
                                       mapping = aes(fill=PAM50_mRNA)
                        
                          ) + 
                          labs(x = "Log2 Expression levels",
                               y = "Patients",
                               title = "Patient gene expression profiles faceted by diagnosed tumor type"
                               ) +
                          xlim(-10, 10 ) +
                          geom_vline(xintercept = c(-1.5, 0.5), 
                                     linetype="dashed"
                                     ) + # based on Control sample profiles
                          theme_bw(base_family = "Times", 
                                   base_size = 10) +
                          theme(legend.position = "right", 
                                axis.text.y = element_blank(),
                                plot.title = element_text(hjust = 0.5,
                                                          size = 15
                                                          )
                                ) +
                          facet_grid(~Tumor
                                     ) 
  
ggsave(plot = plot_TvsPAM50_boxplot, filename = "results/03_EDA_boxplot_T_vs_PAM50_mapping.png",
       device = "png",
)

# TEST for violins
source(file = "R/99_project_functions.R")
p1_boxplot <- plotting_violinplot(data = joined_data_aug, subset_term = "Basal-like", color = "red")
p2_boxplot <- plotting_violinplot(data = joined_data_aug, subset_term = "HER2", color = "green")
p3_boxplot <- plotting_violinplot(data = joined_data_aug, subset_term = "Luminal A", color = "turquoise3")
p4_boxplot <- plotting_violinplot(data = joined_data_aug, subset_term = "Luminal B", color = "purple")


# Call the helper function for legend extraction
shared_legend <- get_legend(p4_boxplot)

# Remove the legend from the remaining plot
p1_boxplot <- p1_boxplot + theme(legend.position = "none")
p2_boxplot <- p2_boxplot + theme(legend.position = "none")
p3_boxplot <- p3_boxplot + theme(legend.position = "none")
p4_boxplot <- p4_boxplot + theme(legend.position = "none")

# Combine the 4 plots and the shared legend
plot_EDA2_boxplot_combo <- grid.arrange(p1_boxplot,
                                        p2_boxplot, 
                                        p3_boxplot, 
                                        p4_boxplot,
                                        shared_legend,
                                        ncol= 2,
                                        widths = c(2.7, 2.7),
                                        nrow = 3,
                                        heights = c(2.7, 2.7, 0.5),
                                        layout_matrix = rbind(c(1,2), 
                                                              c(3,4), 
                                                              c(5,5)) # to bind the legend as a shared 3rd row c(el.5, el.5)
)
