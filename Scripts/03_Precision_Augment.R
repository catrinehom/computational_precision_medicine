# Instructions --------------------------------------------------------------------------------------------------------------------------------------------
# Add new variables to your data


# Load libraries 
library(tidyverse)
#library(data.table)
library(stringr)
# load the data. 
load("Data/lung_data_clean.Rdata")

lung_count_data_clean

# set the rownames of the count data to be the gene names. 
lung_count_data_aug <- lung_count_data_clean %>% 
  tibble::column_to_rownames('gene_names')

as_tibble(lung_count_data_aug)




# count data
ldc <- lung_count_data_clean %>% 
  select(-gene_names) 

# What are colnames for the counts data? 
col_counts <- lung_count_data_clean %>% 
  select(-gene_names) %>% 
  colnames()

col_counts

  # 2. Make the data long.
  #pivot_longer(cols = starts_with("TCGA"), names_to = c("samples")) %>% 
  #distinct(samples, .keep_all = TRUE) 


as_tibble(ldc)
#The count data has 585 columns

lung_pheno_data

# Change the sample type column so it does not include spaces.  
lung_pheno_data <- lung_pheno_data  %>% 
  mutate(sample_type.samples, sample_type.samples =  gsub("y ", "y_", sample_type.samples),
         sample_type.samples, sample_type.samples =  gsub("d ", "d_", sample_type.samples),
         sample_type.samples, sample_type.samples =  gsub("e ", "e_", sample_type.samples),
         sample_type.samples, sample_type.samples =  gsub("t ", "t_", sample_type.samples))

lung_pheno_data

# Phenotype data
lung_pheno_data_1 <- lung_pheno_data %>% 
  select(samples, sample_type.samples) %>% 
  filter(samples %in% colnames(ldc)) %>% 
  mutate(samples= as.factor(samples))

lung_pheno_data_1 


# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# Use the column to rownamees function to set the rownames up!!!. 
lung_pheno_data_aug <- lung_pheno_data_1 %>% 
  tibble::column_to_rownames('samples')

# Note: https://stackoverflow.com/questions/35518218/how-to-set-the-row-names-of-a-data-frame-passed-on-with-the-pipe-operator


# The row_names of the phenotype data 
# rownames(pheno)
rownames(lung_pheno_data_aug)

# The column names of the count data 
# colnames(ensemble_gene_counts)
colnames(lung_count_data_aug)

# save these indices for how to reorder the column names of the count data such that it matches the rownames of the metadata
# genomic_idx <- base::match(rownames(pheno), colnames(ensemble_gene_counts))
# genomic_idx
genomic_idx <- base::match(rownames(lung_pheno_data_aug), colnames(lung_count_data_aug))
genomic_idx

# view the ordered count data. 
# example:
# lung_count_data_ordered  <- ensemble_gene_counts[ , genomic_idx]
# lung_count_data_ordered 
lung_count_data_ordered  <- lung_count_data_aug[ , genomic_idx]
lung_count_data_ordered 


as_tibble(lung_count_data_ordered)
# REORDERING THE DATA SET!!
# Note: https://hbctraining.github.io/Intro-to-R-flipped/lessons/09_reordering-to-match-datasets.html


# View the data. 
head(lung_pheno_data_aug)
head(as_tibble(lung_count_data_ordered))

# Save the augmented data as an RData file. 
save(lung_pheno_data_aug, lung_count_data_ordered, lung_fpkm_data, file = "Data/lung_data_aug.Rdata")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# phenotype data. 
# sample, subtype condition 

# count data 
# samples as columns genes as rows

