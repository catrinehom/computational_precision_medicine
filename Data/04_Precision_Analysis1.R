# Instructions --------------------------------------------------------------------------------------------------------------------------------------------
# Perform Analysis : Differential Expression using DESeq2
# Analysis number 1

# Load libraries 
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(stringr)
library(org.Hs.eg.db)
library(AnnotationDbi)

load("Data/lung_data_aug.Rdata")

# The count matrix and column data can typically be read into R from flat files using base R functions 
# such as read.csv or read.delim. For htseq-count files, see the dedicated input function below.
# With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet.

#write.csv(lung_count_data_ordered,"Path to export the DataFrame\\File Name.csv", row.names = FALSE)
#write.csv(lung_pheno_data_aug,"Path to export the DataFrame\\File Name.csv", row.names = FALSE)

  
# load the data. 
lung_count_data_ordered
# The expression data 
dim(lung_count_data_ordered)

# The phenotype data
dim(lung_pheno_data_aug)

lung_pheno_data_aug %>% 
  filter(sample_type.samples == "NA")

# Differential Expression ---------------------------------------------------------------------------------------------------------------------------------

# Step 1
dds <- DESeqDataSetFromMatrix(
  
  # round the count data because it needs to go in as integers, not floats. 
  countData = round(lung_count_data_ordered),
  
  colData = lung_pheno_data_aug,
  
  # the data is split based on the condition column.
  design = ~sample_type.samples
)

dds


# keep counts that are above or equal to 10 
# keep rows that have at least 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10

# save this as a dds object after sub-setting using keep. 
dds <- dds[keep, ]
dds

cat("10,000 genes just got deleted.")

# Set a factor level to compare between the normal sample and against the Glioblastoma level

# The first condition factor is normal so this has to be set as the reference.
dds$sample_type.samples <- relevel(dds$sample_type.samples, ref = "Solid_Tissue_Normal")
dds$sample_type.samples

# - -------------------------------------------------------------------------------------------------------------------------------------------------------

 
### The null hypothesis is that there is no differential expression across the two sample groups


# - -------------------------------------------------------------------------------------------------------------------------------------------------------


### Run the DESeq Function

# . Step 2
# The differential expression function.
dds <- DESeq(dds)
dds

#* using pre-existing size factors
#* estimating dispersions
#* found already estimated dispersions, replacing these
#* gene-wise dispersion estimates
#* mean-dispersion relationship
#* final dispersion estimates

# .Step 3, Visualise the results
res <- results(dds)
res


#* Recurrent.Tumor vs Solid.Tissue.Normal are compared
#* baseMean the average of the noromalsied counts for all samples.
#* The log2FoldChange of this gene in the GBM condition when compared with the normal condition.
#* The standard error estimates for the log2FoldChange.
#* The Wald test statistic
#* The p-value
#* The corrected p-values for multiple testing.

# Step 4, Look at the summary of the results with an alpha threshold 
summary(res, alpha = 0.05)

res

cat("There are 487 genes upregulateed and there are 2 genes downregulated")

# sort the summary list by p-value.
res2 <- data.frame(res) %>% 
  rownames_to_column("Gene") %>% as_tibble %>% 
  # The gene is more expressed in Treatment if log2FoldChange is more than 1. 
  #filter(padj < 0.05) %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  drop_na() %>% 
  arrange(padj) 
  #head(200)

res2 <- res2 %>% 
  # change the gene length to be 16 characters. 
  mutate(Gene = str_sub(Gene, 1, 15))


res2
# PLOTS ---------------------------------------------------------------------------------------------------------------------------------------------------

# Q: What does the log fold change really mean??

### MA Plot

ggplot(data = res2, 
       mapping = aes(x = baseMean, y = log2FoldChange, col = padj, label = Gene)) + 
  # make a scatter plot
  geom_point(alpha=0.5) + 
  
  ggrepel::geom_label_repel(max.overlaps = 30)+
  
  scale_x_log10() + 
  geom_hline(yintercept = 0.05, alpha = 2.00,color="darkorange2")+
  theme_classic()+
  labs(
    title = "MA Plot",
  )


### Volcano Plot
ggplot(
  data = res2,
  mapping = aes(x = log2FoldChange, y = -log10(pvalue),  col = padj,label = Gene)) +
  
  # scatter plot
  geom_point(size = 3, alpha = 0.7) +
  
  ggrepel::geom_label_repel(max.overlaps = 30)+
  
  # change the axes of the plots. 
 scale_y_continuous(limits = c(3, 5))+
  
  # Add a horizontal line. 
  geom_vline(xintercept = 0, alpha = 0.75, linetype = "dashed") +
  theme_minimal()+
  labs(
    title = "Volcano Plot",
    
  )

# .Step 5 Benjamini Hoschberg Multiple Correction
# sum(res$padj < 0.05, na.rm = TRUE)
# rownames(res)[res$padj < 0.05]


# Differential Expression Analysis Results ----------------------------------------------------------------------------------------------------------------
res2

# Note: Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. 



# Write out the results of differential expression analysis. 
save(res, file = "Data/lung_differential_expression_results_original.Rdata")

res2
save(res2, file = "Data/lung_differential_expression_results.Rdata")





# 2. Map the Ensemble Gene IDs into a SYMBOL Gene Id. 
ets <- AnnotationDbi::select(org.Hs.eg.db, 
                                       key= res2$Gene, 
                                       columns= "SYMBOL",keytype="ENSEMBL") %>% as_tibble()
res2

# Make an MA plot with genes instead of ensemble IDs.
inner_join(res2, ets, by=c("Gene"="ENSEMBL")) %>% 
  ggplot(
         mapping = aes(x = baseMean, y = log2FoldChange, col = padj, label = SYMBOL)) + 
  # make a scatter plot
  geom_point(alpha=0.5) + 
  
  ggrepel::geom_label_repel(max.overlaps = 30)+
  
  scale_x_log10() + 
  geom_hline(yintercept = 0.05, alpha = 2.00,color="darkorange2")+
  theme_classic()+
  labs(
    title = "MA Plot",
  )-> maplot

maplot

ggsave("Doc/images/maplot.png", maplot, width = 5, height = 5)


# Make a volcano plot with genes instead of ensemble IDS. 
inner_join(res2, ets, by=c("Gene"="ENSEMBL")) %>% 
  ggplot(
    mapping = aes(x = log2FoldChange, y = -log10(pvalue),  col = padj,label = SYMBOL)) +
  
  # scatter plot
  geom_point(size = 3, alpha = 0.7) +
  
  ggrepel::geom_label_repel(max.overlaps = 20)+
  
  # change the axes of the plots. 
  scale_y_continuous(limits = c(3, 5))+
  
  # Add a horizontal line. 
  geom_vline(xintercept = 0, alpha = 0.75, linetype = "dashed") +
  theme_minimal()+
  labs(
    title = "Volcano Plot"
  ) -> volc

volc

ggsave("Doc/images/volcano.png", volc, width = 5, height = 5)


# -  ------------------------------------------------------------------------------------------------------------------------------------------------------


# Source: differential expression explanation :https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html

