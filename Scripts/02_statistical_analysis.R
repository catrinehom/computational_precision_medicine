library(tidyverse)
library(dplyr)
library(purrr)
library(ggplot2)

# Load functions
source(file = "Scripts/99_project_functions.R")

# Load data ---------------------------------------------------------------
lung_pheno <- read_tsv(file = "data/01_lung_pheno.tsv",
                          show_col_types = FALSE)
lung_data <- read_tsv(file = "data/01_lung_data.tsv",
                       show_col_types = FALSE)
# Check for NAs --------------------------------------------------------
NAs_pheno <- na_count(lung_pheno) 
#Finds a few NAs

NAs_data <- na_count(as_tibble(lung_data)) 
#No NAs found

# Check distributions --------------------------------------------------------
# Generate a list with the names of the columns that contain numeric values
numeric_pheno <- lung_pheno %>%
  select(where(is.numeric)) %>%
  colnames() %>%
  set_names() #this function belongs to purr package and uses the values of vector as names

# categ_pheno <- lung_pheno %>%
#   select(where(is.character)) %>%
#   colnames() %>%
#   set_names()

# Generate plot iterated for all variables
plots = map(numeric_pheno, 
            ~datadistribution_plot("Expression_Subtype",
                                   ., 
                                   lung_pheno))

# Generate images of the plots
map(names(plots),
    ~save_plot_list("02_",
                    plots, .x)) 

remove(plots)

# PCA plots
pca <- prcomp(t(lung_data), scale=TRUE)
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
df$condition <- lung_pheno$Expression_Subtype

png(file="/Results/plots/PCA_lung_data",)
ggplot(df, aes(x = PC1, y = PC2, color=condition)) +
  geom_point()
dev.off()
