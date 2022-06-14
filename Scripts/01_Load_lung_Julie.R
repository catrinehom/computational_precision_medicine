## Load packages
library(ggplot2)
library(GSVA)
library(ComplexHeatmap)

## Load data
raw <- read.delim("Data/_raw/HiSeqV2-2", row.names=1)
lung_data <- round((2^raw)-1, digits = 0)
colnames(lung_data) <- gsub(pattern = "\\.", replacement = "-", x = colnames(lung_data))
colnames(raw) <- gsub(pattern = "\\.", replacement = "-", x = colnames(raw))


#load pheno data
pheno <- read.delim("Data/_raw/TCGA.LUAD.sampleMap-LUAD_clinicalMatrix",row.names=1)
lung_pheno <- pheno
lung_pheno <- lung_pheno[rownames(lung_pheno) %in% colnames(lung_data),]
lung_pheno <- lung_pheno[match(colnames(lung_data), rownames(lung_pheno)),]
table(rownames(lung_pheno)==colnames(lung_data))

lung_data <- lung_data[, lung_pheno$Expression_Subtype %in% c("Bronchioid", "Squamoid","Magnoid")]
lung_data <- lung_data[ rowSums(lung_data != 0) > 0,]
lung_pheno <- lung_pheno[lung_pheno$Expression_Subtype %in% c("Bronchioid", "Squamoid","Magnoid"),]

lung_data <- apply(lung_data, 2, function(x) x*1) 

lung_data_cpm <- apply(lung_data, 2, function(x) x/sum(x)*1000000) 
df <- data.frame(gene = lung_data_cpm[7,])
ggplot(df, aes(x = gene)) +
  geom_density()

#PCA
pca <- prcomp(t(lung_data_cpm), scale=TRUE)
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
df$condition <- lung_pheno$Expression_Subtype

ggplot(df, aes(x = PC1, y = PC2, color=condition)) +
  geom_point()

## Define a signature for each subtype 
# Using raw, run a Mann-Whitney U test for each gene, for samples in each subtype vs the rest of the subtypes combined
# For each comparison, save genes that are significantly differentially expressed (remember multiple testing correction! Use the function "p.adjust")
# For significantly differentially expressed genes, calculate the log2 fold change (log2(median(samples-in-subtype)/median(samples-in-rest-of-the-subtypes)))
# Pick the highest upregulated genes (either by a log2fc threshold like >1, or by some reasonable size - maybe 40 genes?)
# This is your signature for this subtype


probe_signatures_up <- list()

# Make a function to do a MWW U test
row_mww <- function(x) {
  subtype <- x[lung_pheno$Expression_Subtype == i]
  rest <- x[!lung_pheno$Expression_Subtype == i]
  res <- wilcox.test(subtype, rest)
  return(res$p.value)
}
# Make a function to calculate log2 fc
row_fc <- function(x) {
  subtype <- x[lung_pheno$Expression_Subtype == i]
  rest <- x[!lung_pheno$Expression_Subtype == i]
  res <- log2(median(subtype)/median(rest))
  return(res)
}
# Run the comparisons for each CIT class
for (i in unique(lung_pheno$Expression_Subtype)) {
  pvals <- apply(lung_data, 1, FUN = row_mww)
  padj <- p.adjust(pvals, method = "BY")
  sig_probes <- rownames(lung_data)[padj<0.05]
  subtype_sig <- lung_data[rownames(lung_data) %in% sig_probes,]
  log2fc <- apply(subtype_sig, 1, FUN = row_fc)
  probe_signatures_up[[i]] <- names(log2fc[order(log2fc, decreasing = TRUE)][1:50])
  print(i)
}


column_ha = HeatmapAnnotation(subtype = lung_pheno$Expression_Subtype)
Heatmap(scale(lung_data[rownames(lung_data) %in% unname(unlist(probe_signatures_up)),]), show_column_names = FALSE, show_row_names = FALSE, top_annotation = column_ha)

#ssGSEA

enrichment <- gsva(lung_data,probe_signatures_up,method = "ssgsea", ssgsea.norm = FALSE)
enrichment_subtypes <- rownames(enrichment)[apply(enrichment, 2, which.max)]

table(lung_pheno$Expression_Subtype == rownames(enrichment)[apply(enrichment, 2, which.max)])

## PCA for subgenes
sub_genes = unname(unlist(probe_signatures_up))
lung_data_sub <- lung_data[rownames(lung_data) %in% sub_genes, ]

## Run PCA
pca_sub <- prcomp(t(lung_data_sub), scale = TRUE)
df <- data.frame(PC1 = pca_sub$x[,1], PC2 = pca_sub$x[,2], subtype = lung_pheno$Expression_Subtype)
ggplot(df, aes(x = PC1, y = PC2, color = subtype)) +
  geom_point()

### Distance to centroid

# make empty prediction vector
pred_vector <- c()
# loop over samples (columns)
for(i in 1:ncol(lung_data_sub)) {
  # make a training matrix consisting of all samples but the one you leave out
  training <- CIT_sub[,-i]
  # remove that samples class from the class vector
  training_classes <- CIT_classes[-i]
  # make a test vector consisting of only that sample
  test <- CIT_sub[,i]
  # get the class of that sample
  test_class <- CIT_classes[i]
  
  # make an empty centroid matrix
  centroids <- NULL
  # loop over each of the six classes
  for (class in unique(CIT_classes)) {
    # for each of these classes, subset the training matrix to samples belonging to that class, and calculate the mean expression of each gene in the class
    class_centroid <- rowMeans(training[,training_classes==class])
    # add the mean vector to the centroids matrix
    centroids <- cbind(centroids, class_centroid)
  }
  # add colnames to the centroid matrix
  colnames(centroids) <- unique(CIT_classes)
  # calculate the distance of the test sample to the centroids
  d <- as.matrix(dist(t(cbind(centroids, test))))
  # assign the class of the closest centroid
  class_pred <- names(which.min(d[1:6,7]))
  # check if you got it right and make a logical vector
  pred_vector <- c(pred_vector, test_class==class_pred)
}
# check how many of the cases you got it right
table(pred_vector)
