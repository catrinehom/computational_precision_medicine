library("tidyverse")
library(dplyr)
library(purrr)
library(ggplot2)

# Load functions
source(file = "Scripts/99_project_functions.R")

## Load count data
lung_data <- read.delim("Data/_raw/HiSeqV2-2", row.names=1)
# Rename columns to fit pheno data
colnames(lung_data) <- gsub(pattern = "\\.", replacement = "-", x = colnames(lung_data))
colnames(raw) <- gsub(pattern = "\\.", replacement = "-", x = colnames(raw))

#load pheno data
pheno <- read.delim("Data/_raw/TCGA.LUAD.sampleMap-LUAD_clinicalMatrix",row.names=1)
lung_pheno <- pheno
lung_pheno <- lung_pheno[rownames(lung_pheno) %in% colnames(lung_data),]
lung_pheno <- lung_pheno[match(colnames(lung_data), rownames(lung_pheno)),]

#Clean data 
table(rownames(lung_pheno)==colnames(lung_data))
lung_data <- lung_data[, lung_pheno$Expression_Subtype %in% c("Bronchioid", "Squamoid","Magnoid")]
lung_data <- lung_data[ rowSums(lung_data != 0) > 0,]
lung_pheno <- lung_pheno[lung_pheno$Expression_Subtype %in% c("Bronchioid", "Squamoid","Magnoid"),]
lung_data <- apply(lung_data, 2, function(x) x*1) 

#Replace above with a direct loading of the datasets, after adding a 
#function to save them in the previous part

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