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
#docstring(get_legend)

### Catrine: this makes the helpy thing show up everytime you 
### load this script.

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
#docstring(return_top_genes)

### Catrine: this makes the helpy thing show up everytime you 
### load this script.


# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# ggplot function for handling a single instance
plotting_violinplot <- function(data, subset_term, color) {
  plot <- data %>%   
    select(patient_ID,
           starts_with("NP_"),
           PAM50_mRNA) %>%
    subset(PAM50_mRNA == subset_term) %>%
    pivot_longer(cols = starts_with('NP_')) %>%
    ggplot(aes(y = reorder(patient_ID, value,FUN = median), 
               x = value)
          ) + 
    geom_violin(alpha=0.5,
                scale = "count",
                fill = color
                 ) + 
    labs(x = "Log2 Expression levels",
         y = "Patients",
         title = subset_term
         ) +
    xlim(-10,10) +
    geom_vline(xintercept = c(-1.5, 0.5), 
               linetype="dashed") +  # based on Control sample profiles
    theme_bw(base_family = "Times", 
             base_size = 14) +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 20
                                    )
          )
    return(plot)
}

# ------------------------------------------------------------------------------
# ggplot function for handling a single instance
# This one would add the reference line legend
plotting_boxplot <- function(data, subset_term, color) {
  plot <- data %>%   
    select(patient_ID,
           starts_with("NP_"),
           PAM50_mRNA) %>%
    subset(PAM50_mRNA == subset_term) %>%
    pivot_longer(cols = starts_with('NP_')) %>%
    ggplot(aes(y = reorder(patient_ID, value,FUN = median), 
               x = value)
    ) + 
    geom_boxplot(alpha=0.5,
                 varwidth = TRUE,
                 outlier.shape = NA,
                 fill = color
    ) + 
    labs(x = "Log2 Expression levels",
         y = "Patients",
         title = subset_term
    ) +
    xlim(-10,10) +
    geom_vline(aes(xintercept=c(0.5),
                   color="Range"), 
               linetype="solid",
               size=1) +
    geom_vline(aes(xintercept=c(-1.5),
                   #color="Lower Bound"
                   ), 
               linetype="solid",
               size=1) +
    scale_color_manual(name = "Control samples:", 
                       values = c('Range' = "black", 
                                  'Lower Bound' = "grey12")
                                  ) +
    stat_summary(fun = "median", 
                 geom = "point", 
                 colour = "white", 
                 shape = 23, 
                 size = 1,) +
    theme_bw(base_family = "Times", 
           base_size = 20,) +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 20
          )
    )
  return(plot)
}

# ------------------------------------------------------------------------------



# # Diverging Barcharts
# joined_data_aug %>%
#   select(patient_ID,
#          starts_with("NP_"),
#          PAM50_mRNA) %>%
#   subset(PAM50_mRNA == subset_term) %>%
#   pivot_longer(cols = starts_with('NP_')) %>%


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------