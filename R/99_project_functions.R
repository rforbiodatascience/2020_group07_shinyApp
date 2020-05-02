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
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract top genes from PAM50_mRNA grouped data:
# takes the filterby as astring, a new label for the column and a number of genes to return
# ------------------------------------------------------------------------------
# Data split for each cancer group: without CONTROL samples
return_top_genes <- function(data, filter_by, new_label, n = 10){
  data <- data %>% 
    filter(PAM50_mRNA == filter_by) %>%
    pivot_longer(data = ., cols = -c(PAM50_mRNA, patient_ID)) %>% 
    pivot_wider(data = ., values_from = value, names_from = patient_ID) %>%
    mutate(mean = rowMeans(select(.,-PAM50_mRNA,-name), na.rm = TRUE)) %>%
    mutate(mean = abs(mean)) %>%
    arrange(desc(mean)) %>%
    rename(new_label = name) %>%
    select(new_label) %>%
    slice(.,1:n)
  return(data)
}
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------