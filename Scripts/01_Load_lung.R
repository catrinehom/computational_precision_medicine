## Load packages
library(ggplot2)
library(GSVA)

## Load data
raw <- read.delim("Data/_raw/Lung_study/HiSeqV2-2", row.names=1)
lung_data <- round((2^raw)-1, digits = 0)
colnames(lung_data) <- gsub(pattern = "\\.", replacement = "-", x = colnames(lung_data))

#load pheno data
pheno <- read.delim("Data/_raw/Lung_study/TCGA.LUAD.sampleMap-LUAD_clinicalMatrix",row.names=1)
lung_pheno <- pheno
lung_pheno <- lung_pheno[rownames(lung_pheno) %in% colnames(lung_data),]
lung_pheno <- lung_pheno[match(colnames(lung_data), rownames(lung_pheno)),]
table(rownames(lung_pheno)==colnames(lung_data))

lung_data_cpm <- apply(lung_data, 2, function(x) x/sum(x)*1000000) 
df <- data.frame(gene = lung_data_cpm[5,])
ggplot(df, aes(x = gene)) +
  geom_density()

#PCA
pca <- prcomp(t(lung_data_cpm))
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
df$condition <- lung_pheno$Expression_Subtype
ggplot(df, aes(x = PC1, y = PC2, color=condition)) +
  geom_point()
