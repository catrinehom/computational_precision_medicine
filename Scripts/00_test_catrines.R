## Load necessary packages
library("tidyverse")
library("ggplot2")
library("DESeq2")

## Load and clean expression data
pheno <- read_tsv(file = "/Users/catrinehom/Desktop/Computational_Precision_Medicine/computational_precision_medicine/Data/_raw/TCGA-SKCM.GDC_phenotype.tsv", col_names = TRUE)


pheno <- read_tsv(file = "/Users/catrinehom/Downloads/TCGA.LUNG.sampleMap_LUNG_clinicalMatrix", col_names = TRUE)
