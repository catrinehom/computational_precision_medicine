library("tidyverse")
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
numeric_ones <- lung_pheno %>%
  select(where(is.numeric)) %>%
  colnames() %>%
  set_names() #this function belongs to purr package and uses the values of vector as names

# Generate plot iterated for all variables
plots = map(numeric_ones, 
            ~datadistribution_plot("Expression_Subtype",
                                   ., 
                                   lung_pheno))

# Generate images of the plots
map(names(plots),
    ~save_plot_list("02_",
                    plots, .x)) 

remove(plots)