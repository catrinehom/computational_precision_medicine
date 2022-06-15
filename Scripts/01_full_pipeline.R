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

##############################
# Define Signatures
##############################
# Define how many number of genes to test
no_genes <- seq(5,1000,by=5)

# Define empty pred vectors
pred_gssea <- integer(length(no_genes))
pred_DTC <- integer(length(no_genes))
pred_knn <- integer(length(no_genes))

# Define empty gene signature list
signatures <- list()

# Run the comparisons for each subtype
for (i in unique(lung_pheno$Expression_Subtype)) {
  pvals <- apply(lung_data, 1, FUN = row_mww)
  padj <- p.adjust(pvals, method = "BY")
  sig_probes <- rownames(lung_data)[padj<0.05]
  subtype_sig <- lung_data[rownames(lung_data) %in% sig_probes,]
  log2fc <- apply(subtype_sig, 1, FUN = row_fc)
  signatures[[i]] <- names(log2fc[order(log2fc, decreasing = TRUE)])
}

# Find best number of genes in the signature
for (u in no_genes){
  signatures_up<- signatures 
  signatures_up[["Bronchioid"]] <- signatures_up[["Bronchioid"]][1:u]
  signatures_up[["Magnoid"]] <- signatures_up[["Magnoid"]][1:u]
  signatures_up[["Squamoid"]] <- signatures_up[["Squamoid"]][1:u]
  
  print(cat("Now calcalated number of genes: ", u))
  
  #Define subgenes (across all subtypes) and subset the full dataset
  sub_genes <- unname(unlist(signatures_up))
  lung_data_sub <- lung_data[rownames(lung_data) %in% sub_genes, ]
  
  ##############################
  # Modeling
  ##############################
  
  ###############
  # SSGEA
  ###############
  # Find score for each signature in each sample
  enrichment <- gsva(lung_data,signatures_up,method = "ssgsea", ssgsea.norm = FALSE)
  # Find the highest signature score for each sample
  enrichment_subtypes <- rownames(enrichment)[apply(enrichment, 2, which.max)]
  
  # answer_dtc is for all samples for a specific u
  answer_gssea <- table(lung_pheno$Expression_Subtype == rownames(enrichment)[apply(enrichment, 2, which.max)])
  
  # Remember this score for next loop
  pred_gssea[u/5] <- answer_gssea["TRUE"]/(answer_gssea["FALSE"]+answer_gssea["TRUE"])
  
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
    
    subtype_knn <- as.character(knn(train_knn, test_knn, cl, k = 10, prob=TRUE))
    answer_knn[i] = subtype_knn == t(lung_pheno$Expression_Subtype[i])
    
  }
  # Remember this score for next loop
  pred_knn[u/5] <- sum(answer_knn/ncol(lung_data_sub))
  
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
    d <- as.matrix(dist(t(cbind(centroids, test))))
    # assign the class of the closest centroid
    class_pred <- names(which.min(d[1:3,4]))
    # check if you got it right and make a logical vector
    answer_dtc <- c(answer_dtc, test_class==class_pred)
  }
  # check how many of the cases you got it right
  answer_dtc <-table(answer_dtc)
  # Remember this score for next loop
  pred_DTC[u/5] <- answer_dtc["TRUE"]/(answer_dtc["FALSE"]+answer_dtc["TRUE"])
}

##############################
# Prediction value of models
##############################

#Combine pred score from the three models
prediction <-cbind(no_genes,pred_DTC,pred_knn,pred_gssea)
colnames(prediction) <- c("#genes","DTC","kNN","GSSEA")

#Plot the prediction scores
df <- data.frame(prediction)
ggplot(df, aes(x= no_genes, y = Prediction , color = variable)) +
  geom_point(aes(y=pred_DTC,col="DTC"))+
  geom_point(aes(y=pred_knn,col="kNN"))+
  geom_point(aes(y=pred_gssea,col="GSSEA"))