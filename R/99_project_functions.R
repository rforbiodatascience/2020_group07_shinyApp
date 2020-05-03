# Load libraries
# ------------------------------------------------------------------------------

if (!requireNamespace("docstring", quietly = TRUE))
  install.packages("docstring")
library(docstring)

# Define project functions
# ------------------------------------------------------------------------------
as_matrix <- function(x){
  if(!tibble::is_tibble(x) ) stop("x must be a tibble")
  y <- as.matrix.data.frame(x[,-1])
  rownames(y) <- x[[1]]
  y
}
# Credit: https://rdrr.io/github/HuntsmanCancerInstitute/hciR/man/as_matrix.html
# ------------------------------------------------------------------------------
# Extract the plot's legend
library(gridExtra)
get_legend <- function(myggplot){
  #' Returns the "legend" object from a plot object:
  #' 
  #' @param myggplot a ggplot object
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Wrap for `?function_name``
docstring(get_legend)

# Extract top genes from PAM50_mRNA grouped data:
# takes the filterby as astring, a new label for the column and a number of genes to return
# ------------------------------------------------------------------------------
# Data split for each cancer group: without CONTROL samples
return_top_genes <- function(data, filter_by, new_label, n = 10){
  #' Data wrangling helper function:
  #' 
  #' Returns the top n differentially expressed genes for a subset of data
  #' 
  #' @param filter_by specifying the subset string for filter()
  #' @param new_label specify the new column name in the final tible
  #' @param n the number of top expressed genes
  data <- data %>% 
    select(patient_ID, PAM50_mRNA, starts_with("NP_")) %>%
    filter(PAM50_mRNA == filter_by) %>%
    pivot_longer(data = ., cols = -c(PAM50_mRNA, patient_ID)) %>% 
    pivot_wider(data = ., values_from = value, names_from = patient_ID) %>%
    mutate(avg = rowMeans(select(.,-PAM50_mRNA,-name), na.rm = TRUE)) %>%
    mutate(avg = abs(avg)) %>%
    arrange(desc(avg)) %>%
    select(name) %>%
    rename(., new_label = name) %>%
    slice(.,1:n)
  return(data)
}

# Wrap for `?function_name``
docstring(return_top_genes)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------