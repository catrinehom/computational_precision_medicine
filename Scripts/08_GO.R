library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationHub)

#Define the list of background genes
backg<- as.character(row.names(lung_data))

#Define the list of genes of interest (Upregulated genes and down regulated genes)
u = 200/2
signatures_up <- signatures_all_up
signatures_up[["Bronchioid"]] <- signatures_all_up[["Bronchioid"]][1:u]
signatures_up[["Magnoid"]] <- signatures_all_up[["Magnoid"]][1:u]
signatures_up[["Squamoid"]] <- signatures_all_up[["Squamoid"]][1:u]

signatures_down <- signatures_all_down
signatures_down[["Bronchioid"]] <- signatures_all_down[["Bronchioid"]][1:u]
signatures_down[["Magnoid"]] <- signatures_all_down[["Magnoid"]][1:u]
signatures_down[["Squamoid"]] <- signatures_all_down[["Squamoid"]][1:u]

sub_genes_up <- unname(c(unlist(signatures_up)))                    
sub_genes_down <- unname(c(unlist(signatures_down))) 

save(sub_genes_up, file="Results/signatures_up_names.Rdata")
save(sub_genes_down, file="Results/signatures_down_names.Rdata")


sigUP<- sub_genes_up
sigDOWN<- sub_genes_down    

sigUP <- gsub(pattern = "\\.[0-9]*$", replacement = "", x = sigUP)
colnames(lung_data) <- gsub(pattern = ".$", replacement = "", x = colnames(lung_data))
                    
#ORA of gene ontology terms (biological processes)
#Upregulated - this command can take some time to run
ego <- enrichGO(gene = sigUP,
                universe = backg,
                keyType = "ENSEMBL", 
                OrgDb = org.Hs.eg.db,
                ont = "BP", #BP=biological process, write all to get all the GO
                pAdjustMethod = "BH", #bonferoni
                qvalueCutoff = 0.05, 
                pvalueCutoff = 0.05,
                readable = TRUE)
