# Scaled down version of this blog post on the TensorFlow for R Blog at RStudio
# https://blogs.rstudio.com/tensorflow/posts/2018-01-29-dl-for-cancer-immunotherapy/

# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library('tidyverse')
library('keras')

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
clincal_data_aug <- read_csv(file = "data/01_clincal_data_clean.csv")
PAM50_aug <- read_csv(file = "data/01_PAM50_clean.csv")
proteome_data_aug <- read_csv(file = "data/02_proteome_data_wide_aug.csv")
joined_data_aug <- read_csv(file = "data/02_joined_data_clean.csv")

# View the data
proteome_data_aug %>% print

# View class distribution
joined_data_aug %>% count(PAM50_mRNA) %>% print

# Prepare data
# ------------------------------------------------------------------------------
# Partition into data_type where test (~10%) and training set (~90%)

joined_data_aug <- joined_data_aug  %>% 
  mutate(data_type = sample(10, size = nrow(.), replace = TRUE))  %>% 
  mutate(PAM50_mRNA_bin = case_when(PAM50_mRNA == "Basal-like" ~ 0,
                                    PAM50_mRNA == "HER2-enriched" ~ 1,
                                    PAM50_mRNA == "Luminal A" ~ 2,
                                    PAM50_mRNA == "Luminal B" ~ 3)) %>%
  select(patient_ID, starts_with("NP"),data_type, PAM50_mRNA_bin) %>%
  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) # Take median value of NA values


joined_data_aug %>% count(PAM50_mRNA_bin) %>% print


# Define training and test feature matrices
X_train = joined_data_aug %>%
  filter(data_type > 3) %>%
  select(patient_ID, starts_with("NP")) %>%
  as_matrix()
         
X_test = joined_data_aug %>%
  filter(data_type <= 3) %>%
  select(patient_ID, starts_with("NP")) %>%
  as_matrix()
  



# Define known target classes for training and test data
y_train = joined_data_aug %>%
  filter(data_type > 3 ) %>%
  pull(PAM50_mRNA_bin) %>% 
  to_categorical
y_test = joined_data_aug %>%
  filter(data_type <= 3 ) %>%
  pull(PAM50_mRNA_bin) %>% 
  to_categorical


# Define ANN model
# ------------------------------------------------------------------------------

# Set hyperparameters
n_hidden_1 = 43
h1_activate = 'relu'
drop_out_1 = 0.4
n_hidden_2 = 40
h2_activate = 'relu'
drop_out_2 = 0.3
n_hidden_3 = 20
h3_activate = 'relu'
drop_out_3 = 0.2
n_hidden_4 = 10
h4_activate = 'relu'
drop_out_4 = 0.1
n_output   = 4
o_ativate  = 'softmax'
n_epochs = 150
batch_size = 50
loss_func = 'categorical_crossentropy'
learn_rate = 0.001

# Set architecture
model = keras_model_sequential() %>% 
  layer_dense(units = n_hidden_1, activation = h1_activate, input_shape = 43) %>% 
  layer_dropout(rate = drop_out_1) %>% 
  layer_dense(units = n_hidden_2, activation = h2_activate) %>%
  layer_dropout(rate = drop_out_2) %>%
  layer_dense(units = n_hidden_3, activation = h3_activate) %>%
  layer_dropout(rate = drop_out_3) %>%
  layer_dense(units = n_hidden_4, activation = h4_activate) %>%
  layer_dropout(rate = drop_out_4) %>%
  layer_dense(units = n_output, activation = o_ativate)

# Compile model
model %>%
  compile(loss = loss_func,
          optimizer = optimizer_rmsprop(lr = learn_rate),
          metrics = c('accuracy')
  )

# View model
model %>% summary %>% print

# Train model
# ------------------------------------------------------------------------------

history = model %>%
  fit(x = X_train,
      y = y_train,
      epochs = n_epochs,
      batch_size = batch_size,
      validation_split = 0
  )

# Evaluate model
# ------------------------------------------------------------------------------
# OBS THIS ONLY WORKS IF ALL CLASSES IS PRESENT IN y_test 
# OTHERWISE, THE FACTORING GOES WRONG.
perf_test = model %>% evaluate(X_test, y_test)
acc_test = perf_test %>% pluck('acc') %>% round(3) * 100
perf_train = model %>% evaluate(X_train, y_train)
acc_train = perf_train %>% pluck('acc') %>% round(3) * 100
results = bind_rows(
  tibble(y_true = y_test %>%
           apply(1, function(x){ return( which(x==1) - 1) }) %>%
           factor,
         y_pred = model %>%
           predict_classes(X_test) %>%
           factor,
         Correct = ifelse(y_true == y_pred ,"yes", "no") %>%
           factor,
         data_type = 'test')
  ,
  tibble(y_true = y_train %>%
           apply(1, function(x){ return( which(x==1) - 1) }) %>%
           factor,
         y_pred = model %>%
           predict_classes(X_train) %>%
           factor,
         Correct = ifelse(y_true == y_pred ,"yes", "no") %>%
           factor,
         data_type = 'train'))
my_counts = results %>% count(y_pred, y_true, data_type)

# Visualise model performance
# ------------------------------------------------------------------------------
title = paste0('Performance of Deep Feed Forward Neural Network (',
               'Total number of model parameters = ', count_params(model), ').')
sub_title = paste0("Test Accuracy = ", acc_test, "%, n = ", nrow(X_test), ". ",
                   "Training Accuracy = ", acc_train, "%, n = ", nrow(X_train), ".")
xlab  = 'Predicted (Class assigned by Keras/TensorFlow deep FFN)'
ylab  = 'Measured (Real class)'
results %>%
  ggplot(aes(x = y_pred, y = y_true, fill = Correct)) +
  geom_jitter(pch = 21, size = 4, alpha = 0.4, colour = 'black') +
  geom_text(data = my_counts, aes(x = y_pred, y = y_true, label = n),
            size = 20, inherit.aes = FALSE) +
  xlab(xlab) +
  ylab(ylab) +
  ggtitle(label = title, subtitle = sub_title) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(labels = c('No', 'Yes'),
                     values = c('tomato','cornflowerblue')) +
  facet_wrap(~data_type, nrow = 1)

# Save results
#ggsave(filename = "results/ANN_performance.png",device = "png")

# Save model
# ------------------------------------------------------------------------------
#save_model_hdf5(object = model,
#                filepath = "Models/05_peptide_model.h5")
