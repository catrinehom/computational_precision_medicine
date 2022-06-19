# Instructions --------------------------------------------------------------------------------------------------------------------------------------------
# Perform Analysis: Principal Components Analysis

library(tidyverse)
library(broom)
library(factoextra)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggthemes)
library(ggrepel)
library(pheatmap)

#  load the data ------------------------------------------------------------------------------------------------------------------------------------------

load("Data/lung_data_aug.Rdata")

# Subset the expression matrix to only the differentially expressed genes. 
load("Data/lung_differential_expression_results_original.Rdata")


# Load all genes from DE
res <- data.frame(res)
res

# Filter for up/down genes. 
res <- res %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 4) %>% 
  arrange(padj) %>% 
  drop_na()
 
res

dim(res)

# Find only the row names of the genes. 
Gene <- rownames(res)
Gene 

# Change the column label. 
gene_rownames <- as_tibble(Gene)
gene_rownames

gene_rownames <- gene_rownames %>% 
  rename("value"="Gene")

gene_rownames

# Bind the rownames in a usable dataframe. 
signif_genes_from_DE <- gene_rownames %>% 
  cbind(as_tibble(res))



# Create the counts object --------------------------------------------------------------------------------------------------------------------------------
# Extract the count data for the significantly expressed genes. 
signif_genes_from_DE

lung_count_data_ordered
#rownames(lung_count_data_ordered)

lung_sig <- signif_genes_from_DE %>% 
  tibble::column_to_rownames('Gene')
lung_sig

rownames(lung_sig)

# match the rows in the count data to the rows in significantly expressed genes from DESeq2
genomic_idx <- base::match( rownames(lung_count_data_ordered), rownames(lung_sig))
genomic_idx

# Extract the count data for the significantly expressed genes. 
lu <- lung_count_data_ordered[genomic_idx, ]
lu


# Lung Gene Count Data ------------------------------------------------------------------------------------------------------------------------------------


# save the matrix with a good name (IMPORTANT OBJECT)
lung_expr_deg <- lu %>% 
  na.omit()

lung_expr_deg 
#as_tibble(lung_expr_deg)
dim(lung_expr_deg)


# The genes are the features across the top. (rotated)
lung_expr_deg_rot <- t(lung_expr_deg)
lung_expr_deg_rot

dim(lung_expr_deg_rot)
as_tibble(lung_expr_deg_rot)


# Change the Gene names   --------------------------------------------------------------------------------------------------------------

#1. get the row names 
rownames(lung_expr_deg)

#. 2. 

Gene_Rownames  <- rownames(lung_expr_deg) %>% 
  as_tibble() %>% 
  rename("value" = "Gene")

# Cut off the decimal point... 
lung_expr_deg_gn <- Gene_Rownames %>% 
  cbind(as_tibble(lung_expr_deg)) %>% 
  mutate(Gene = str_sub(Gene, 1, 15))

as_tibble(lung_expr_deg_gn)

# 2. Map the Ensemble Gene IDs into a SYMBOL Gene Id. 
ens_to_symbol <- AnnotationDbi::select(org.Hs.eg.db, 
                                       key= lung_expr_deg_gn$Gene, 
                                       columns= "SYMBOL",keytype="ENSEMBL")

ens_to_symbol <- as_tibble(ens_to_symbol)
ens_to_symbol

lung_expr_deg_joined <- inner_join(lung_expr_deg_gn, ens_to_symbol, by=c("Gene"="ENSEMBL")) %>% 
  relocate(SYMBOL)

# Create the new matrix of genes. 
lung_exp_deg_matx <- lung_expr_deg_joined %>% 
  dplyr::select(-Gene) %>% 
  drop_na()

lung_exp_deg_matx

# change the genes into rownames for unsupervised learning. 

lung_expr_DE <- column_to_rownames(lung_exp_deg_matx, var = "SYMBOL")
lung_expr_DE

as_tibble(lung_expr_DE)

# - -------------------------------------------------------------------------------------------------------------------------------------------------------





# normalisation..? 
# lung_expr_cpm <- apply(lung_expr_deg, 2 , function(x) x/sum(x)*100000)
# lung_expr_cpm


# Matching the Phenotype Data -----------------------------------------------------------------------------------------------------------------------------

dim(lung_pheno_data_aug)

rownames(lung_pheno_data_aug)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------





# gbm_expr_deg <- gbm_expr_cpm[rownames(gbm_expr_cpm) %in% rownames(res)[res$padj < 0.05], ]
# as_tibble(gbm_expr_deg)

# Principal Components Analysis ---------------------------------------------------------------------------------------------------------------------------

# 1. Perform PCA 
pca_deg <- prcomp(lung_expr_deg)
pca_deg

# 2. Check out how much variance is explained
summary(pca_deg)


# 3. put the summary data in a tidy fashion. 
pca_summary <- pca_deg %>%
  tidy(matrix = "eigenvalues") %>% 
  mutate(percentage = percent * 100)
pca_summary


# 4. Create a bar plot from this information
# Make a bar plot of the principal components. 
ggplot( data = pca_summary,
        
        mapping = aes(x = PC, y = percent)) +
  
  # make a bar plot
  geom_col(colour = "darkblue", fill = "#56B4E8", alpha = 0.3) +
  
  # Add a line and then some points. 
  geom_line()+
  
  geom_point(shape=21, color="black", fill="#69b3a2", size=2) +
  
  # Adjust the x axis and the y axis. 
  scale_x_continuous(breaks = 1:13, limits = c(0,20)) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.01))) +
  
  theme_classic(base_family = "Helvetica",
                base_size = 12) +
  
  # Add some labels. 
  labs(
    y = "Variance explained by each PC",
    x = "The Principal Component",
    title = "Principal Componnts Analysis",
    subtitle = "scree plot",
    caption = "Fig") 


# Create a table of PC1 and PC2 of the genes. 
df <- data.frame(PC1 = pca_deg$x[, 1], PC2 = pca_deg$x[, 2])
df

df_chan <- rownames_to_column(df, "Gene") %>% 
  mutate(Gene = str_sub(Gene, 1, 15))

etos <- AnnotationDbi::select(org.Hs.eg.db, 
                      key = df_chan$Gene , 
                      columns = "SYMBOL",keytype="ENSEMBL")
etos


df2 <- inner_join(df_chan, etos, 
                  by=c("Gene"="ENSEMBL")) %>% 
  relocate(SYMBOL)

as_tibble(df2)

# Make a plot of the principal components (without colour)
ggplot(data = df2, 
       mapping = aes(x = PC1, y = PC2, label = SYMBOL)) +
  geom_point(pch = 21, size = 4, alpha = 0.3, fill = "darkblue")+
  ggrepel::geom_label_repel(max.overlaps = 20)+
  ggthemes::scale_fill_fivethirtyeight()+
  theme_minimal()+
  labs(
    title = "Principal Component Analysis",
    x = "85.7 % Variance explained by PC1",
    y = "5.69 % Variance explained by PC2"
  )

# Hierarchical Clustering  --------------------------------------------------------------------------------------------------------------------------------

# Genes are on the columns here. 
gene_matrix <- dist(lung_expr_DE)

my_distances <- dist(gene_matrix, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

my_tree <- hclust(my_distances, method = "centroid")


# Install R Skittl brewer from dectools 
devtools::install_github('alyssafrazee/RSkittleBrewer')
mm_colors = RSkittleBrewer::RSkittleBrewer('M&M')

dend_plot <- fviz_dend(my_tree, 
          k = 6,
          color_labels_by_k = FALSE,
          kcolors = ggthemes::scale_colour_fivethirtyeight(alpha = 0.7),
          cex = 0.7,
          main = "Heirarchical Cluster Dendrogram",
          ylab = "The Euclidean Distance between the Genes", 
          subtitle = "Centroid Method, Unsupervised Machine Learning", 
          xlab = "Genes",
          horiz = FALSE,
          ggtheme = theme_clean())
 
dend_plot

# require("igraph")
# dend_plot <- fviz_dend(my_tree, k = 10,
#                        k_colors = "uchicago",
#                        color_labels_by_k = T,
#                        type = "rectangle", 
#                        cx = 0.7,
#                        lwd = 1.5,
#                        main = "Heirarchical Cluster Dendrogram - Average Method",
#                        ylab = "The Euclidean Distance between the Genes", 
#                        subtitle = "Unsupervised Machine Learning", 
#                        xlab = "Genes",
#                        horiz = FALSE,
#                        ggtheme = theme_clean())



ggsave("Doc/images/heirarchical.png", dend_plot, width = 5, height = 5)


# HEATMAP -----------------------------------------------------------------------------------------------------------------------------------------------

heat_comp <- pheatmap(lung_expr_DE, 
         color=viridis::viridis(17),
         fontsize_col = 1,
         fontsize_row = 11,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Heatmap",
         clustering_method = "centroid")
  


ggsave("Doc/images/heat.png",heat_comp, width = 5, height = 5)


# colours for plots: https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html

