# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("broom")
library("gridExtra")


# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data
# ------------------------------------------------------------------------------
#clincal_data_aug <- read_csv(file = "data/01_clinical_data_clean.csv")
#PAM50_clean <- read_csv(file = "data/01_PAM50_clean.csv")
joined_data_aug <- read_csv(file = "data/02_joined_data_PAM50_aug.csv")


# Wrangle data
# ------------------------------------------------------------------------------
joined_data_aug <- joined_data_aug %>%
  filter(Class != "Control") 


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



# Plot: Age distribution in different cancer types
# ------------------------------------------------------------------------------
joined_data_aug %>% 
  ggplot(mapping = aes(x = Age_at_Initial_Pathologic_Diagnosis, 
                       fill = Class)) +
  geom_histogram(binwidth = 10) +
  scale_x_continuous(breaks = seq(20, 100, 10)) + 
  labs(title = "Age distribution in different cancer types",
       x = 'Age',
       y = 'Count') +
  #scale_fill_manual(values = custom_colors, ) +
  theme_bw(base_family = "Times", 
           base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, 
                                size = 25)) 

ggsave(filename = "results/03_EDA_age_distribution.png", device = "png")


# Plot: COME BACK TO IT
# ------------------------------------------------------------------------------
joined_data_aug %>%
  ggplot(aes(x = AJCC_Stage, 
             y = Age_at_Initial_Pathologic_Diagnosis, 
             colour = AJCC_Stage)) +
  geom_violin(scale = "area") + 
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5) +
  stat_summary(fun=median, geom="point", 
               size=2, color="black") +
  theme_bw(base_family = "Times", 
           base_size = 12) +
  labs(title = "Age versus cancer severity",
       x = 'AJJC stage',
       y = 'Age at initial diagnosis')

ggsave(filename = "results/03_EDA_age_cancer_severity.png",device = "png")


# Plot: Gender distribution 
# ------------------------------------------------------------------------------
joined_data_aug %>% 
  ggplot(mapping = aes(Gender)) +
  geom_bar() +
  theme_bw(base_family = "Times", 
           base_size = 12) +
  labs(y = "Count",
       title = "Gender distribution") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14))

ggsave(filename = "results/03_EDA_gender_vs_tumortype.png", device = "png")   

# Plot: Class distribution across patients: REORDER THE BARS
# ------------------------------------------------------------------------------
joined_data_aug %>% 
  ggplot(mapping = aes(Class, fill = Class)) +
  geom_bar() +
  theme_bw(base_family = "Times", 
           base_size = 12) +
  labs(y = "Count",
       title = "Class distribution across patients") +
  #scale_fill_manual(values = custom_colors, ) +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14))

ggsave(filename = "results/03_EDA_class_distribution.png", device = "png")    


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


# BOXPLOTS: combining 4x plots into one canvas
# ------------------------------------------------------------------------------
# Mappings between the subsets and the plotting function:
# Creating a separate plot while gettign colored individually, to be combined later
p1_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "Basal", color = "red")
p2_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "HER2", color = "green")
p3_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "LumA", color = "turquoise3")
p4_boxplot <- plotting_boxplot(data = joined_data_aug, subset_term = "LumB", color = "purple")


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
                                                              c(5,5))) # to bind the legend as a shared 3rd row c(el.5, el.5)

ggsave(plot = plot_EDA2_boxplot_combo, filename = "results/03_EDA_boxplot_combined.png",
       device = "png",
       height = 5)

### Catrine: the black line is the control samples... median +/- something?
### Why is some samples outside the range?

# ------------------------------------------------------------------------------
plot_TvsPAM50_boxplot <- joined_data_aug %>%   
  select(patient_ID,
         Tumor,
         starts_with("NP_"),
         Class
  ) %>%
  filter(Class != "Control") %>%
  pivot_longer(cols = starts_with('NP_')) %>%
  ggplot(aes( y = patient_ID, 
              x = value,
  )
  ) + 
  geom_boxplot(alpha= 0.5,
               varwidth = TRUE,
               outlier.shape = NA,
               mapping = aes(fill=Class)
               
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
p1_boxplot <- plotting_violinplot(data = joined_data_aug, subset_term = "Basal", color = "red")
p2_boxplot <- plotting_violinplot(data = joined_data_aug, subset_term = "HER2", color = "green")
p3_boxplot <- plotting_violinplot(data = joined_data_aug, subset_term = "LumA", color = "turquoise3")
p4_boxplot <- plotting_violinplot(data = joined_data_aug, subset_term = "LumB", color = "purple")


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
                                                              c(5,5))) # to bind the legend as a shared 3rd row c(el.5, el.5)