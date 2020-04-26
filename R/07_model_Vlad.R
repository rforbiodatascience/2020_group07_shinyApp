# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library(tidymodels)
library(tidyverse)
library(workflows)
library(keras)
library(tune)
if (!require('tidymodels')) install.packages('tidymodels'); library('tidymodels')
if (!require('keras')) install.packages('keras'); library('keras')
if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('tensorflow')) install.packages('tensorflow'); library('tensorflow')

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
proteome_data_aug <- read_csv(file = "data/02_proteome_data_wide_aug.csv")
df <- read_csv(file = "data/02_joined_data_aug.csv")

# Splitting data
# ------------------------------------------------------------------------------
# create a data partition
set.seed(234589)
# split the data into trainng (75%) and testing (25%)
df_split <- initial_split(df, 
                         prop = 1/5)


# extract training and testing sets
train_idx <- training(df_split)
test_idx <- testing(df_split)


# At some point we’re going to want to do some parameter tuning, and to do that 
# we’re going to want to use cross-validation. So we can create a cross-validated 
# version of the training set in preparation for that moment using vfold_cv().

# create CV object from training data
df_cv <- vfold_cv(train_idx)

# Define the recipy()

# define the recipe
logReg_recipy <- df %>%
 select(Tumor, starts_with("NP_")) %>%
 recipe(Tumor ~ ., data =  df) %>%
 step_normalize(all_numeric()) %>%
 step_knnimpute(all_predictors())



# Specify the model
# So far we’ve split our data into training/testing, and we’ve specified our pre-processing
# steps using a recipe. The next thing we want to specify is our model (using the parsnip package).
# Parsnip offers a unified interface for the massive variety of models that exist in R.
# This means that you only have to learn one way of specifying a model, and you can use 
# this specification and have it generate a linear model, a random forest model, a support vector machine model, and more with a single line of code.


# I want to try a logistic regression based on "Tumor" class & expression data
# with  "L2" penalty regularization using "keras" library ?!
logReg_model <- 
 # specify that the model is a random forest
 logistic_reg() %>%
 # specify that the `mtry` parameter needs to be tuned
 set_args(penalty="L2") %>% # additionally a tune() could be tried for for tuning regularization weights of 'parameters" parameter
 # select the engine/package that underlies the model
 set_engine("glmnet") %>%
 # choose either the continuous regression or binary classification mode
 set_mode("classification")


# Put it all together in a workflow
# We’re now ready to put the model and recipes together into a workflow. 
# You initiate a workflow using workflow() (from the workflows package) and then you
# can add a recipe and add a model to it.

# set the workflow
logReg_workflow <- workflow() %>%
 # add the recipe
 add_recipe(logReg_recipy) %>%
 # add the model
 add_model(logReg_model)

# An additional step can be performed if we are using the tuning workflow for our parameters
# which chould be included here
# Example code just in acase
# specify which values eant to try
# 
# rf_grid <- expand.grid(mtry = c(3, 4, 5))
# # extract results
# rf_tune_results <- rf_workflow %>%
#   tune_grid(resamples = diabetes_cv, #CV object
#             grid = rf_grid, # grid of values to try
#             metrics = metric_set(accuracy, roc_auc) # metrics we care about
#   )
# 
# Then you would need to add that tuned parameter
# param_final <- rf_tune_results %>%
#
# select_best(metric = "accuracy", maximize = TRUE)
# param_final
# 
# rf_workflow <- rf_workflow %>%
#   finalize_workflow(param_final)
#

# Fit the final model
# Now we’ve defined our recipe, our model, and tuned the model’s parameters, 
# we’re ready to actually fit the final model. Since all of this information is 
# contained within the workflow object, we will apply the last_fit() function 
# to our workflow and our train/test split object. This will automatically train the model 
# specified by the workflow using the training data, 
# and produce evaluations based on the test set.

logReg_fit <- logReg_workflow %>%
 # fit on entire training set and evaluate on test set
 last_fit(df_split)

# Wrangle data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% ...

# Model data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% ...

# Write data
# ------------------------------------------------------------------------------
write_tsv(file = "data/05_my_data_clean_aug_anl_modelX.tsv")
ggsave(...)