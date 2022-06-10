library("DESeq2")
library("ggplot2")
library("tidyverse")
library("dplyr")
library("readr")

raw <- read.table(file = gzfile('Data/_raw/TCGA-SKCM.htseq_counts.tsv.gz','rt'), sep = '\t', header = TRUE)
phenotype <- read_tsv(file = gzfile('Data/_raw/TCGA-SKCM.GDC_phenotype.tsv.gz'), col_names=TRUE)
