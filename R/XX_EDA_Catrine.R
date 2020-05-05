# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library('tidyverse')

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
clinical_data_clean <- read_csv(file = "data/01_clinical_data_clean.csv")




clinical_data_clean %>%
  ggplot(aes(x = AJCC_Stage, 
             y = Age_at_Initial_Pathologic_Diagnosis, 
             colour = AJCC_Stage)) +
  geom_violin() + 
  geom_jitter(alpha = 0.5) +
  stat_summary(fun=median, geom="point", 
               size=2, color="black") +
  theme_bw(base_family = "Times", 
           base_size = 12) +
  labs(title = "Age versus cancer severity",
       x = 'AJJC stage',
       y = 'Age at initial diagnosis')

ggsave(filename = "results/03_EDA_age_cancer_severity.png",device = "png")

