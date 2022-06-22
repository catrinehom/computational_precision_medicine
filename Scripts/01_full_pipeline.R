##############################
# Load and clean data
##############################
# Load packages
library(ggplot2)
library(GSVA)
library(ComplexHeatmap)
if(!require('class')) {
  install.packages('class')
  library('class')
}
library(readr)
library(tidyverse)
library(rankdist)

# Load functions
source(file = "Scripts/99_project_functions.R")

## Load count data
lung_data <- read.delim("Data/_raw/TCGA-LUAD.htseq_fpkm.tsv", row.names=1)
lung_data <- round((2^lung_data)-1, digits = 0)
# Rename columns to fit pheno data
colnames(lung_data) <- gsub(pattern = "\\.", replacement = "-", x = colnames(lung_data))
colnames(lung_data) <- gsub(pattern = ".$", replacement = "", x = colnames(lung_data))

#load pheno data
pheno <- read.delim("Data/_raw/TCGA.LUAD.sampleMap-LUAD_clinicalMatrix",row.names=1)
lung_pheno <- pheno
lung_pheno <- lung_pheno[rownames(lung_pheno) %in% colnames(lung_data),]
lung_pheno <- lung_pheno[match(colnames(lung_data), rownames(lung_pheno)),]
remove(pheno) #Not used anymore

#Clean data 
table(rownames(lung_pheno)==colnames(lung_data))
lung_data <- lung_data[, lung_pheno$Expression_Subtype %in% c("Bronchioid", "Squamoid","Magnoid")]
lung_data <- lung_data[ rowSums(lung_data != 0) > 0,]
lung_pheno <- lung_pheno[lung_pheno$Expression_Subtype %in% c("Bronchioid", "Squamoid","Magnoid"),]

lung_data_remove_outliers <- lung_data[ , -which(names(lung_data) %in% c("TCGA-44-6147-01.1","TCGA-44-5645-01","TCGA-44-6146-01","TCGA-44-2662-01.1","TCGA-44-2666-01.1","TCGA-44-2656-01.1","TCGA-44-6775-01.1","TCGA-44-4112-01"))]
lung_data_remove_outliers <- lung_data_remove_outliers[ rowSums(lung_data_remove_outliers != 0) > 0,]

lung_pheno_remove_outliers <- lung_pheno[-which(row.names(lung_pheno) %in% c("TCGA-44-6147-01.1","TCGA-44-5645-01","TCGA-44-6146-01","TCGA-44-2662-01.1","TCGA-44-2666-01.1","TCGA-44-2656-01.1","TCGA-44-6775-01.1","TCGA-44-4112-01")),]

# Write data --------------------------------------------------------------
write_tsv(x = lung_pheno,
          file = "Data/01_lung_pheno.tsv")

write_tsv(x = as_tibble(lung_data),
          file = "Data/01_lung_data.tsv")

##############################
# Remove outliers
##############################

pca <- prcomp(t(lung_data), scale=TRUE)
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
df$Subtype <- lung_pheno$Expression_Subtype

ggplot(df, aes(x = PC1, y = PC2, color=Subtype)) +
  geom_point(size = 3, alpha = 1) + 
  theme(text = element_text(size = 20)) +
  xlab("PC1 (24.28%)") + ylab("PC2 (4.83%)")

ggsave("Results/plots/PCA_full.png")

# We can see 8 patient being weird
pca <- prcomp(t(lung_data_remove_outliers), scale=TRUE)
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
df$Subtype <- lung_pheno_remove_outliers$Expression_Subtype

ggplot(df, aes(x = PC1, y = PC2,color=Subtype)) +
  geom_point(size = 3, alpha = 1) +
  theme(text = element_text(size = 20)) +
  xlab("PC1 (5.81%)") + ylab("PC2 (5.40%)")

ggsave("Results/plots/PCA_outlier_removal.png")


lung_data <- lung_data_remove_outliers
lung_data <- apply(lung_data, 2, function(x) x*1) 
lung_pheno <- lung_pheno_remove_outliers

##############################
# Define Signatures
##############################
# Define how many number of genes to test
interval <- 2
no_genes <- seq(2,1000,by=interval)

# Define empty pred vectors
pred_gssea <- integer(length(no_genes))
pred_DTC <- integer(length(no_genes))
pred_knn <- integer(length(no_genes))

# Define empty gene signature list
signatures_all_up <- list()
signatures_all_down <- list()

# Run the comparisons for each subtype for upregulated
for (i in unique(lung_pheno$Expression_Subtype)) {
  pvals <- apply(lung_data, 1, FUN = row_mww)
  padj <- p.adjust(pvals, method = "BY")
  sig_probes <- rownames(lung_data)[padj<0.05]
  subtype_sig <- lung_data[rownames(lung_data) %in% sig_probes,]
  log2fc <- apply(subtype_sig, 1, FUN = row_fc)
  signatures_all_up[[i]] <- names(log2fc[order(log2fc, decreasing = TRUE)])
  signatures_all_down[[i]] <- names(log2fc[order(log2fc, decreasing = FALSE)])
}

# Rank count data
lung_data_unranked <- lung_data
lung_data <- apply(lung_data, 2, rank)

# Find best number of genes in the signature
for (u in no_genes){
  signatures_up <- signatures_all_up
  signatures_up[["Bronchioid"]] <- signatures_all_up[["Bronchioid"]][1:u]
  signatures_up[["Magnoid"]] <- signatures_all_up[["Magnoid"]][1:u]
  signatures_up[["Squamoid"]] <- signatures_all_up[["Squamoid"]][1:u]
  
  signatures_down <- signatures_all_down
  signatures_down[["Bronchioid"]] <- signatures_all_down[["Bronchioid"]][1:u]
  signatures_down[["Magnoid"]] <- signatures_all_down[["Magnoid"]][1:u]
  signatures_down[["Squamoid"]] <- signatures_all_down[["Squamoid"]][1:u]
  
  print(cat("Now calcalated number of genes: ", u))
  
  #Define subgenes (across all subtypes) and subset the full dataset
  sub_genes <- unname(c(unlist(signatures_up),unlist(signatures_down)))
  lung_data_sub <- lung_data[rownames(lung_data) %in% sub_genes, ]
  
  ##############################
  # Modeling
  ##############################
  
  ###############
  # SSGEA
  ###############
  # Find score for each up signature in each sample
  enrichment <- gsva(lung_data,signatures_up,method = "ssgsea", ssgsea.norm = FALSE)
  # Find score for each down signature in each sample
  detraction <- gsva(lung_data,signatures_down,method = "ssgsea", ssgsea.norm = FALSE)
  detraction <- detraction *-1
  
  total <- enrichment + detraction
  
  # Find the highest signature score for each sample
  enrichment_subtypes <- rownames(total)[apply(total, 2, which.max)]
  
  # answer_dtc is for all samples for a specific u
  answer_gssea <- table(lung_pheno$Expression_Subtype == rownames(enrichment)[apply(enrichment, 2, which.max)])
  
  # Remember this score for next loop
  pred_gssea[u/interval] <- answer_gssea["TRUE"]/(answer_gssea["FALSE"]+answer_gssea["TRUE"])
  
  ###############
  # KNN
  ###############
  # answer_knn is for all samples for a specific u
  answer_knn = integer(ncol(lung_data_sub))
  for(i in 1:ncol(lung_data_sub)){
    test_knn <- subset(lung_data_sub, select = c(i)) 
    train_knn <- subset(lung_data_sub, select = -c(i)) 
    
    cl <- subset(t(lung_pheno$Expression_Subtype), select = -c(i)) 
    cl <- factor(t(cl))
    train_knn <- t(train_knn)
    test_knn <- t(test_knn)
    
    subtype_knn <- as.character(knn(train_knn, test_knn, cl, k = 17, prob=TRUE))
    answer_knn[i] = subtype_knn == t(lung_pheno$Expression_Subtype[i])
    
  }
  # Remember this score for next loop
  pred_knn[u/interval] <- sum(answer_knn/ncol(lung_data_sub))
  
  ###############
  # DTC
  ###############
  # answer_dtc is for all samples for a specific u
  answer_dtc <- c()
  # loop over samples (columns)
  for(i in 1:ncol(lung_data_sub)) {
    # make a training matrix consisting of all samples but the one you leave out
    training <- lung_data_sub[,-i]
    # remove that samples class from the class vector
    training_classes <- lung_pheno$Expression_Subtype[-i]
    # make a test vector consisting of only that sample
    test <- lung_data_sub[,i]
    # get the class of that sample
    test_class <- lung_pheno$Expression_Subtype[i]
    
    # make an empty centroid matrix
    centroids <- NULL
    # loop over each of the three classes
    for (class in unique(lung_pheno$Expression_Subtype)) {
      # for each of these classes, subset the training matrix to samples belonging to that class, and calculate the mean expression of each gene in the class
      class_centroid <- rowMeans(training[,training_classes==class])
      # add the mean vector to the centroids matrix
      centroids <- cbind(centroids, class_centroid)
    }
    # add colnames to the centroid matrix
    colnames(centroids) <- unique(lung_pheno$Expression_Subtype)
    # calculate the distance of the test sample to the centroids
    d <- as.matrix(dist(t(cbind(centroids, test)),method="canberra"))
    # assign the class of the closest centroid
    class_pred <- names(which.min(d[1:3,4]))
    # check if you got it right and make a logical vector
    answer_dtc <- c(answer_dtc, test_class==class_pred)
  }
  # check how many of the cases you got it right
  answer_dtc <-table(answer_dtc)
  # Remember this score for next loop
  pred_DTC[u/interval] <- answer_dtc["TRUE"]/(answer_dtc["FALSE"]+answer_dtc["TRUE"])
}

##############################
# Prediction value of models
##############################
no_genes<-no_genes*2
#Combine pred score from the three models
prediction <-cbind(no_genes,pred_DTC,pred_knn,pred_gssea)
colnames(prediction) <- c("#genes","DTC","kNN","GSSEA")

#Plot the prediction scores
df <- data.frame(prediction)
ggplot(df, aes(x= no_genes, y = Prediction , color = Model)) +
  geom_point(aes(y=DTC,col="DTC"))+
  geom_point(aes(y=kNN,col="kNN"))+
  geom_point(aes(y=GSSEA,col="ssGSEA")) +
  theme(text = element_text(size = 20)) +
  xlab("Number of genes in signature") + ylab("Prediction score")
  

ggsave("Results/plots/01_model_predictions.png")

##############################
# Save prediction value of models
##############################

save(prediction, file="Results/prediction.Rdata")


#DTC: 0.8218182, 2000
#KNN:  0.8254545, 200
#SSGSEA:  0.7454545, 168

##############################
# Run confusion matrix on top results
##############################

###############
# SSGEA
###############
u = 168/2
signatures_up <- signatures_all_up
signatures_up[["Bronchioid"]] <- signatures_all_up[["Bronchioid"]][1:u]
signatures_up[["Magnoid"]] <- signatures_all_up[["Magnoid"]][1:u]
signatures_up[["Squamoid"]] <- signatures_all_up[["Squamoid"]][1:u]

signatures_down <- signatures_all_down
signatures_down[["Bronchioid"]] <- signatures_all_down[["Bronchioid"]][1:u]
signatures_down[["Magnoid"]] <- signatures_all_down[["Magnoid"]][1:u]
signatures_down[["Squamoid"]] <- signatures_all_down[["Squamoid"]][1:u]

# Find score for each up signature in each sample
enrichment <- gsva(lung_data,signatures_up,method = "ssgsea", ssgsea.norm = FALSE)
# Find score for each down signature in each sample
detraction <- gsva(lung_data,signatures_down,method = "ssgsea", ssgsea.norm = FALSE)
detraction <- detraction *-1

total <- enrichment + detraction

# Find the highest signature score for each sample
enrichment_subtypes <- rownames(total)[apply(total, 2, which.max)]

# answer_dtc is for all samples for a specific u
answer_gssea <- table(lung_pheno$Expression_Subtype == rownames(enrichment)[apply(enrichment, 2, which.max)])

confusionmatrix <- cbind(c(0,0,0),c(0,0,0),c(0,0,0))
colnames(confusionmatrix) <- c("Bronchioid","Magnoid","Squamoid")
rownames(confusionmatrix) <- c("Bronchioid","Magnoid","Squamoid")

for(i in 1:length(enrichment_subtypes)) {
confusionmatrix[enrichment_subtypes[i],lung_pheno$Expression_Subtype[i]] <-confusionmatrix[enrichment_subtypes[i],lung_pheno$Expression_Subtype[i]] +1 
}

precition_Bronchioid <- confusionmatrix["Bronchioid","Bronchioid"]/rowSums(confusionmatrix)[1]
precition_Magnoid <- confusionmatrix["Magnoid","Magnoid"]/rowSums(confusionmatrix)[2]
precition_Squamoid <- confusionmatrix["Squamoid","Squamoid"]/rowSums(confusionmatrix)[3]

recall_Bronchioid <- confusionmatrix["Bronchioid","Bronchioid"]/colSums(confusionmatrix)[1]
recall_Magnoid <- confusionmatrix["Magnoid","Magnoid"]/colSums(confusionmatrix)[2]
recall_Squamoid <- confusionmatrix["Squamoid","Squamoid"]/colSums(confusionmatrix)[3]

cbind(precition_Bronchioid,precition_Magnoid,precition_Squamoid)
cbind(recall_Bronchioid,recall_Magnoid,recall_Squamoid)

###############
# KNN
###############
u = 200/2
signatures_up <- signatures_all_up
signatures_up[["Bronchioid"]] <- signatures_all_up[["Bronchioid"]][1:u]
signatures_up[["Magnoid"]] <- signatures_all_up[["Magnoid"]][1:u]
signatures_up[["Squamoid"]] <- signatures_all_up[["Squamoid"]][1:u]

signatures_down <- signatures_all_down
signatures_down[["Bronchioid"]] <- signatures_all_down[["Bronchioid"]][1:u]
signatures_down[["Magnoid"]] <- signatures_all_down[["Magnoid"]][1:u]
signatures_down[["Squamoid"]] <- signatures_all_down[["Squamoid"]][1:u]

save(signatures_down, file="Results/signatures_down.Rdata")
save(signatures_down, file="Results/signatures_up.Rdata")

sub_genes <- unname(c(unlist(signatures_up),unlist(signatures_down)))
lung_data_sub <- lung_data[rownames(lung_data) %in% sub_genes, ]

save(lung_data_sub, file="Results/lung_data_sub.Rdata")

confusionmatrix <- cbind(c(0,0,0),c(0,0,0),c(0,0,0))
colnames(confusionmatrix) <- c("Bronchioid","Magnoid","Squamoid")
rownames(confusionmatrix) <- c("Bronchioid","Magnoid","Squamoid")

answer_knn = integer(ncol(lung_data_sub))
for(i in 1:ncol(lung_data_sub)){
  test_knn <- subset(lung_data_sub, select = c(i)) 
  train_knn <- subset(lung_data_sub, select = -c(i)) 
  
  cl <- subset(t(lung_pheno$Expression_Subtype), select = -c(i)) 
  cl <- factor(t(cl))
  train_knn <- t(train_knn)
  test_knn <- t(test_knn)
  
  subtype_knn <- as.character(knn(train_knn, test_knn, cl, k = 17, prob=TRUE))
  answer_knn[i] = subtype_knn == t(lung_pheno$Expression_Subtype[i])
  
  confusionmatrix[subtype_knn,lung_pheno$Expression_Subtype[i]] <-confusionmatrix[subtype_knn,lung_pheno$Expression_Subtype[i]] +1 
}

precition_Bronchioid <- confusionmatrix["Bronchioid","Bronchioid"]/rowSums(confusionmatrix)[1]
precition_Magnoid <- confusionmatrix["Magnoid","Magnoid"]/rowSums(confusionmatrix)[2]
precition_Squamoid <- confusionmatrix["Squamoid","Squamoid"]/rowSums(confusionmatrix)[3]

recall_Bronchioid <- confusionmatrix["Bronchioid","Bronchioid"]/colSums(confusionmatrix)[1]
recall_Magnoid <- confusionmatrix["Magnoid","Magnoid"]/colSums(confusionmatrix)[2]
recall_Squamoid <- confusionmatrix["Squamoid","Squamoid"]/colSums(confusionmatrix)[3]

cbind(precition_Bronchioid,precition_Magnoid,precition_Squamoid)
cbind(recall_Bronchioid,recall_Magnoid,recall_Squamoid)


# Plot subset 
pca <- prcomp(t(lung_data_sub), scale=TRUE)
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
df$Subtype <- lung_pheno$Expression_Subtype

ggplot(df, aes(x = PC1, y = PC2, color=Subtype)) +
  geom_point(size = 3, alpha = 1) +
  theme(text = element_text(size = 20)) +
  xlab("PC1 (18.39%)") + ylab("PC2 (10.08%)")


ggsave("Results/plots/PCA_subset.png")

# plot subset unranked
lung_data_sub_unranked <- lung_data_unranked[rownames(lung_data_unranked) %in% sub_genes, ]
pca <- prcomp(t(lung_data_sub_unranked), scale=TRUE)
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
df$condition <- lung_pheno$Expression_Subtype

ggplot(df, aes(x = PC1, y = PC2, color=condition)) +
  geom_point(size = 3, alpha = 1)

###############
# DTC
###############
u = 2000/2
signatures_up <- signatures_all_up
signatures_up[["Bronchioid"]] <- signatures_all_up[["Bronchioid"]][1:u]
signatures_up[["Magnoid"]] <- signatures_all_up[["Magnoid"]][1:u]
signatures_up[["Squamoid"]] <- signatures_all_up[["Squamoid"]][1:u]

signatures_down <- signatures_all_down
signatures_down[["Bronchioid"]] <- signatures_all_down[["Bronchioid"]][1:u]
signatures_down[["Magnoid"]] <- signatures_all_down[["Magnoid"]][1:u]
signatures_down[["Squamoid"]] <- signatures_all_down[["Squamoid"]][1:u]

sub_genes <- unname(c(unlist(signatures_up),unlist(signatures_down)))
lung_data_sub <- lung_data[rownames(lung_data) %in% sub_genes, ]

confusionmatrix <- cbind(c(0,0,0),c(0,0,0),c(0,0,0))
colnames(confusionmatrix) <- c("Bronchioid","Magnoid","Squamoid")
rownames(confusionmatrix) <- c("Bronchioid","Magnoid","Squamoid")


# answer_dtc is for all samples for a specific u
answer_dtc <- c()
# loop over samples (columns)
for(i in 1:ncol(lung_data_sub)) {
  # make a training matrix consisting of all samples but the one you leave out
  training <- lung_data_sub[,-i]
  # remove that samples class from the class vector
  training_classes <- lung_pheno$Expression_Subtype[-i]
  # make a test vector consisting of only that sample
  test <- lung_data_sub[,i]
  # get the class of that sample
  test_class <- lung_pheno$Expression_Subtype[i]
  
  # make an empty centroid matrix
  centroids <- NULL
  # loop over each of the three classes
  for (class in unique(lung_pheno$Expression_Subtype)) {
    # for each of these classes, subset the training matrix to samples belonging to that class, and calculate the mean expression of each gene in the class
    class_centroid <- rowMeans(training[,training_classes==class])
    # add the mean vector to the centroids matrix
    centroids <- cbind(centroids, class_centroid)
  }
  # add colnames to the centroid matrix
  colnames(centroids) <- unique(lung_pheno$Expression_Subtype)
  # calculate the distance of the test sample to the centroids
  d <- as.matrix(dist(t(cbind(centroids, test)),method="canberra"))
  # assign the class of the closest centroid
  class_pred <- names(which.min(d[1:3,4]))
  # check if you got it right and make a logical vector
  answer_dtc <- c(answer_dtc, test_class==class_pred)
  
  confusionmatrix[class_pred,lung_pheno$Expression_Subtype[i]] <-confusionmatrix[class_pred,lung_pheno$Expression_Subtype[i]] +1 
  
  
}


precition_Bronchioid <- confusionmatrix["Bronchioid","Bronchioid"]/rowSums(confusionmatrix)[1]
precition_Magnoid <- confusionmatrix["Magnoid","Magnoid"]/rowSums(confusionmatrix)[2]
precition_Squamoid <- confusionmatrix["Squamoid","Squamoid"]/rowSums(confusionmatrix)[3]

recall_Bronchioid <- confusionmatrix["Bronchioid","Bronchioid"]/colSums(confusionmatrix)[1]
recall_Magnoid <- confusionmatrix["Magnoid","Magnoid"]/colSums(confusionmatrix)[2]
recall_Squamoid <- confusionmatrix["Squamoid","Squamoid"]/colSums(confusionmatrix)[3]

cbind(precition_Bronchioid,precition_Magnoid,precition_Squamoid)
cbind(recall_Bronchioid,recall_Magnoid,recall_Squamoid)
