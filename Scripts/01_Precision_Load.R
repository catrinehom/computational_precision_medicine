

### Original DataSet Chosen : TCGA lung adenocarcinoma (LUAD) gene expression by RNAseq (There is no phenotype data here)
#https://xenabrowser.net/datapages/?dataset=TCGA.LUAD.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# Log2(x+1),  x is the RSEM value

###  Second Dataset Chosen : GDC TCGA Lung Adenocarcinoma (LUAD) Cohort :The FPKM Data
# https://xenabrowser.net/datapages/?dataset=TCGA-LUAD.htseq_fpkm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

# - -------------------------------------------------------------------------------------------------------------------------------------------------------
# Final Datasets Selected ---------------------------------------------------------------------------------------------------------------------------------



# Cohort GDC TCGA Lung Adenocarcinoma (LUAD)
# * High Throughput Sequencing Counts: TCGA-LUAD.htseq_counts.tsv 
# * The associated Phenotype Data: TCGA-LUAD.GDC_phenotype.tsv.gz
# * The FPKM Data: TCGA-LUAD.htseq_fpkm.ts

# NOTE: Download this data from :
# GDC TCGA Lung Adenocarcinoma (LUAD) > Cohort > RNA Seq Section found @
# https://xenabrowser.net


# Instructions  -------------------------------------------------------------
# Collapse data to a single file or convert  to .tsv



# LOAD THE DATA -------------------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)

# load the counts data for lung Adenocarcinoma
lung_count_data <- read.table("Data/_raw/TCGA-LUAD.htseq_counts.tsv", header = TRUE, row.names = 1)
lung_count_data 



# look at the structure of the data. 
dim(lung_count_data)
cat("There are 60,488 rows and 586 columns")


# load the phenotype data for lung Adenocarcinoma
lung_phenotype_data <- read_tsv("Data/_raw/TCGA-LUAD.GDC_phenotype.tsv.gz")
lung_phenotype_data


# load the fpkm data for lung Adenocarcinoma
lung_fpkm_data <- read_tsv("Data/_fpkm/TCGA-LUAD.htseq_fpkm.tsv")
lung_fpkm_data 

dim(lung_fpkm_data)

# save the data objects into one file. 
save(lung_count_data , lung_phenotype_data, lung_fpkm_data , file = "Data/lung_data.Rdata")



