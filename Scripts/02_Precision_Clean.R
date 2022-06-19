# Instructions  -------------------------------------------------------------------------------------------------------------------------------------------

# Remove invalid data, e.g. if you have amino acid sequence data,
# remove non-valid sequences containing X or other non-standard amino acid characters or fix columns,
# e.g. dates or when two labels are the same, but spelled differently


# Load libraries 
library(tidyverse)
library(stringr)


#  Count Data ---------------------------------------------------------------------------------------------------------------------------------------------

load("Data/lung_data.Rdata")
lung_count_data


# save the row names which are the gene names in ENSEMBLE and (edit them)
gene_names <-rownames(lung_count_data)

# convert the matrix into a tibble. 
count_data <- as_tibble(lung_count_data)

# combine the gene names and gene expression data. 
lung_count_data_wg <- cbind(gene_names, count_data) 
lung_count_data_wg

# change the column names in the data to replace the . with - using base R. 
colnames(lung_count_data_wg) <- gsub("\\.", "-", colnames(lung_count_data_wg))
lung_count_data_wg


# Note:https://stackoverflow.com/questions/68926419/how-to-edit-all-column-names-to-replace-a-certain-character-in-r

#lung_count_data_clean <- lung_count_data_wg %>% 
  # select all columns except the gene column
  #pivot_longer(cols = starts_with("TCGA"), names_to = c("samples")) %>% 
  # replace the . with a - ?????? check this. 
  #mutate(samples = str_replace_all("([\\.)", "-", samples)) %>% 
  # put them back into wide format. 
  #pivot_wider(names_from = samples)


# write out th final count data. 
lung_count_data_clean <- lung_count_data_wg

# Phenotype Data.  ----------------------------------------------------------------------------------------------------------------------------------------


# Instead of 125 columns in the phenotype data, we can choose 10 of interest. 
lung_phenotype_data_filt <- lung_phenotype_data %>% 
  dplyr::select(submitter_id.samples,  disease_type, sample_type.samples, sample_type_id.samples,
         age_at_initial_pathologic_diagnosis, 
         followup_treatment_success,
         person_neoplasm_cancer_status, 
         year_of_initial_pathologic_diagnosis, 
         gender.demographic,
         primary_diagnosis.diagnoses)


# Change some of the column names for the phenotype data. 
lung_pheno_data <- lung_phenotype_data_filt %>% 
  
  rename(samples = submitter_id.samples) %>% 
   
  # make most variables factors. 
   dplyr::mutate(disease_type = as.factor(disease_type),
         sample_type.samples = as.factor(sample_type.samples),
         person_neoplasm_cancer_status = as.factor(person_neoplasm_cancer_status),
         followup_treatment_success = as.factor(followup_treatment_success),
         gender.demographic = as.factor(gender.demographic),
         primary_diagnosis.diagnoses = as.factor(primary_diagnosis.diagnoses))

lung_pheno_data 

save(lung_pheno_data, lung_count_data_clean, lung_fpkm_data, file = "Data/lung_data_clean.Rdata")

