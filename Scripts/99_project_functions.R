################################
# Functions for Wilcox
################################

# Make a function to do a MWW U test
row_mww <- function(x) {
  subtype <- x[lung_pheno$Expression_Subtype == i]
  rest <- x[!lung_pheno$Expression_Subtype == i]
  res <- wilcox.test(subtype, rest)
  return(res$p.value)
}
# Make a function to calculate log2 fc
row_fc <- function(x) {
  subtype <- x[lung_pheno$Expression_Subtype == i]
  rest <- x[!lung_pheno$Expression_Subtype == i]
  res <- log2(median(subtype)/median(rest))
  return(res)
}

# Check for and summarize NA values in a dataframe
na_count <- function(df) {
  df %>% 
    summarise_all(~ sum(is.na(.)))
}

# Rotate x-axis labels by 45 deg
rotate_x <- function() {
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
}

# Rotate y-axis labels by 45 deg
rotate_y <- function() {
  theme(axis.text.y = element_text(angle = 45,
                                   hjust = 1))
}

# Plotting distribution of variables
datadistribution_plot <- function(x, y, df){
  ggplot(df, aes(x = .data[[x]], 
                 y = .data[[y]],
                 fill = Expression_Subtype)) + 
    geom_violin() +
    geom_boxplot(aes(x = .data[[x]], 
                     y = .data[[y]]), 
                 width = 0.3)
}

# Save plots in images to plots folder in results 

folder <- "results/plots"

if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}

image_path = "results/plots"

save_plot_list <- function(prefix, plot_list, .x) {
  ggsave(filename = paste0(prefix, .x, ".png"),
         path = image_path,
         plot = plot_list[[.x]],
         height = 7,
         width = 7,
         unit = "in")
}