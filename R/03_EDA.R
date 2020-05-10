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
joined_data_aug <- read_csv(file = "data/02_joined_data_PAM50_aug.csv")
joined_data_full_aug <-  read_csv(file = "data/02_joined_data_full_aug.csv")

# Wrangle data
# ------------------------------------------------------------------------------
# Get the Control sample range
control_range <- joined_data_full_aug %>%
  filter(Class != "Control") %>%
  select(patient_ID, starts_with("NP_")) %>%
  pivot_longer(cols = -patient_ID) %>% 
  group_by(patient_ID) %>% 
  summarise(lower = quantile(value, 0.25), 
            upper = quantile(value, 0.75)) %>% 
  summarise(min = min(lower), 
            max = max(upper)) %>%
  unlist()

# Remove controls, since there are too few samples in this group
joined_data_aug <- joined_data_aug %>%
  filter(Class != "Control") 


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
           base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, 
                                size = 18)) 

ggsave(filename = "results/03_EDA_age_distribution.png", device = "png")


# Plot: Gender distribution 
# ------------------------------------------------------------------------------
joined_data_aug %>% 
  ggplot(mapping = aes(Gender)) +
  geom_bar() +
  theme_bw(base_family = "Times", 
           base_size = 15) +
  labs(y = "Count",
       title = "Gender distribution") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 18))

ggsave(filename = "results/03_EDA_gender_vs_tumortype.png", device = "png")   

# Plot: Class distribution across patients
# ------------------------------------------------------------------------------
joined_data_aug %>% 
  ggplot(mapping = aes(Class, fill = Class)) +
  geom_bar() +
  theme_bw(base_family = "Times", 
           base_size = 12) +
  labs(y = "Count",
       title = "Class distribution across patients") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 18))

ggsave(filename = "results/03_EDA_class_distribution.png", device = "png")    


# BOXPLOTS: combining 4x plots into one canvas
# ------------------------------------------------------------------------------
# Mappings between the subsets and the plotting function:
# Creating a separate plot while gettign colored individually, to be combined later
p1_boxplot <- plotting_boxplot(data = joined_data_full_aug, 
                               subset_term = "Basal", 
                               color = "red",
                               control_range = control_range)
p2_boxplot <- plotting_boxplot(data = joined_data_full_aug, 
                               subset_term = "HER2", 
                               color = "green",
                               control_range = control_range)
p3_boxplot <- plotting_boxplot(data = joined_data_full_aug, 
                               subset_term = "LumA", 
                               color = "turquoise3",
                               control_range = control_range)
p4_boxplot <- plotting_boxplot(data = joined_data_full_aug, 
                               subset_term = "LumB", 
                               color = "purple",
                               control_range = control_range)


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
                                        # To bind the legend as a shared 3rd row c(el.5, el.5)
                                        layout_matrix = rbind(c(1,2), 
                                                              c(3,4), 
                                                              c(5,5))) 

ggsave(plot = plot_EDA2_boxplot_combo, filename = "results/03_EDA_boxplot_combined.png",
       device = "png",
       height = 5)
