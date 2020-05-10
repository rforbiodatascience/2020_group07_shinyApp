# DEFINE PROJECT FUNCTIONS
# ------------------------------------------------------------------------------


# Tibble to matrix
# ------------------------------------------------------------------------------
# Credit: https://rdrr.io/github/HuntsmanCancerInstitute/hciR/man/as_matrix.html
as_matrix <- function(x){
  if(!tibble::is_tibble(x) ) stop("x must be a tibble")
  y <- as.matrix.data.frame(x[,-1])
  rownames(y) <- x[[1]]
  y
}


# Extract the plot's legend
# ------------------------------------------------------------------------------
get_legend <- function(myggplot){
  #' Returns the "legend" object from a plot object:
  #' 
  #' @param myggplot a ggplot object
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# Customized boxplot
# ------------------------------------------------------------------------------
plotting_boxplot <- function(data, subset_term, color, 
                             control_range = control_range, 
                             control_range_colour = "yellow1") {
  plot <- data %>%   
    select(patient_ID,
           starts_with("NP_"),
           Class) %>%
    subset(Class == subset_term) %>%
    pivot_longer(cols = starts_with('NP_')) %>%
    ggplot(aes(y = reorder(patient_ID, value,FUN = median), 
               x = value)) + 
    geom_boxplot(alpha=0.5,
                 varwidth = TRUE,
                 outlier.shape = NA,
                 fill = color) + 
    labs(x = "Log2 Expression levels",
         y = "Patients",
         title = subset_term) +
    xlim(-10,10) +
    geom_vline(aes(xintercept=control_range[2],
                   color="Control samples"), 
               linetype="solid",
               size=1) +
    geom_vline(aes(xintercept=control_range[1]), 
               linetype="solid",
               size=1) +
    scale_color_manual(name = "Inter Quartile Range:", 
                       values = c('Control samples' = "black", 
                                  'Lower Bound' = "grey12")) +
    stat_summary(fun = "median", 
                 geom = "point", 
                 colour = "white", 
                 shape = 23, 
                 size = 1) +
    theme_bw(base_family = "Times", 
             base_size = 20) +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 20)) +
    annotate("rect", 
             xmin= control_range[1], xmax=control_range[2], # supply from data
             ymin=-Inf, ymax=Inf,
             alpha=0.3, fill= control_range_colour) 
  return(plot)
}



# Set up the plot's theme to be consistent across all plots
# ------------------------------------------------------------------------------
myplot_aes <- theme_bw(base_family = "Times", 
                         base_size = 18) +
              theme(plot.title = element_text(hjust = 0.5, 
                                              size = 25),
                    legend.position = "bottom")

