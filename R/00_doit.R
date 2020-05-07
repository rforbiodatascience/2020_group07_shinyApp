# Install packages
# ------------------------------------------------------------------------------
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse") 
if (!requireNamespace("broom", quietly = TRUE))
  install.packages("broom") 
if (!requireNamespace("gridExtra", quietly = TRUE))
  install.packages("gridExtra") 


# Run all scripts
# ------------------------------------------------------------------------------
source(file = "R/01_load_and_clean.R.R")
source(file = "R/02_augment.R")
source(file = "R/03_EDA.R")
source(file = "R/04_PCA_clustering.R")
source(file = "R/05_ANN_model.R")


# Catrine:
# Do install packages for all scripts?

# should we use same theme for all plots?:
#theme_bw(base_family = "Times", 
#         base_size = 12)

### Tumor groups
# Doesnt make sense to plot them
#https://cancerstaging.org/references-tools/quickreferences/Documents/BreastMedium.pdf?fbclid=IwAR0nEcUZfzS7DrIhi1GiUXyBO-uyxtFS5PocRxpP7INBnW8HnrLY9oBnYq4
#https://www.cancer.org/cancer/breast-cancer/understanding-a-breast-cancer-diagnosis/breast-cancer-hormone-receptor-status.html?fbclid=IwAR3kMPbjhFbOxB87h4_soWBaVUHf4P1kZsEDZDyM_9_TqrgbBYuHfUaB-wg
