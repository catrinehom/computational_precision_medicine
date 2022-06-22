library(fgsea)
library(reactome.db)
library(org.Hs.eg.db)
library(stringr)
library(dplyr)

#Load necessary data-------------------------------

load("Results/signatures_down.Rdata")
load("Results/signatures_up.Rdata")
lung_data <- read.delim("Data/_raw/TCGA-LUAD.htseq_fpkm.tsv", row.names=1) #Maybe find a faster way to load the initial gene set

magnoid_down <- signatures_down$Magnoid
squamoid_down <- signatures_down$Squamoid
bronchioid_down <- signatures_down$Bronchioid

magnoid_up <- signatures_up$Magnoid
squamoid_up <- signatures_up$Squamoid
bronchioid_up <- signatures_up$Bronchioid

universe <- as_tibble(rownames(lung_data))

remove(signatures_down, signatures_up, lung_data)

## Convert ENSEMBL to ENTREZID--------------------------------------------------

# Convert lists to tibbles
magnoid_down <- as_tibble(magnoid_down)
magnoid_down <- magnoid_down %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))

squamoid_down <- as_tibble(squamoid_down)
squamoid_down <- squamoid_down %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))

bronchioid_down <- as_tibble(bronchioid_down)
bronchioid_down <- bronchioid_down %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))

magnoid_up <- as_tibble(magnoid_up)
magnoid_up <- magnoid_up %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))

squamoid_up <- as_tibble(squamoid_up)
squamoid_up <- squamoid_up %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))

bronchioid_up <- as_tibble(bronchioid_up)
bronchioid_up <- bronchioid_up %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))

universe <- universe %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))
  
hs <- org.Hs.eg.db #Load conversion dataset

# Convert ENSEMBL to ENTREZ
magnoid_down <- AnnotationDbi::select(hs, 
                                 keys = magnoid_down$genes,
                                 columns = c("ENTREZID"),
                                 keytype = "ENSEMBL")

squamoid_down <- AnnotationDbi::select(hs,
                                  keys = squamoid_down$genes,
                                  columns = c("ENTREZID"),
                                  keytype = "ENSEMBL")

bronchioid_down <- AnnotationDbi::select(hs,
                                    keys = bronchioid_down$genes,
                                    columns = c("ENTREZID"),
                                    keytype = "ENSEMBL")

magnoid_up <- AnnotationDbi::select(hs, 
                                      keys = magnoid_up$genes,
                                      columns = c("ENTREZID"),
                                      keytype = "ENSEMBL")

squamoid_up <- AnnotationDbi::select(hs,
                                       keys = squamoid_up$genes,
                                       columns = c("ENTREZID"),
                                       keytype = "ENSEMBL")

bronchioid_up <- AnnotationDbi::select(hs,
                                         keys = bronchioid_up$genes,
                                         columns = c("ENTREZID"),
                                         keytype = "ENSEMBL")

universe <- AnnotationDbi::select(hs,
                                  keys = universe$genes,
                                  columns = c("ENTREZID"),
                                  keytype = "ENSEMBL")

# Check symbols---------------------------------------------------------
# magnoid_down <- AnnotationDbi::select(hs, 
#                                       keys = magnoid_down$ENSEMBL,
#                                       columns = c("SYMBOL"),
#                                       keytype = "ENSEMBL")
# 
# squamoid_down <- AnnotationDbi::select(hs,
#                                        keys = squamoid_down$ENSEMBL,
#                                        columns = c("SYMBOL"),
#                                        keytype = "ENSEMBL")
# 
# bronchioid_down <- AnnotationDbi::select(hs,
#                                          keys = bronchioid_down$ENSEMBL,
#                                          columns = c("SYMBOL"),
#                                          keytype = "ENSEMBL")
# 
# magnoid_up <- AnnotationDbi::select(hs, 
#                                     keys = magnoid_up$ENSEMBL,
#                                     columns = c("SYMBOL"),
#                                     keytype = "ENSEMBL")
# 
# squamoid_up <- AnnotationDbi::select(hs,
#                                      keys = squamoid_up$ENSEMBL,
#                                      columns = c("SYMBOL"),
#                                      keytype = "ENSEMBL")
# 
# bronchioid_up <- AnnotationDbi::select(hs,
#                                        keys = bronchioid_up$ENSEMBL,
#                                        columns = c("SYMBOL"),
#                                        keytype = "ENSEMBL")
# 
# universe <- AnnotationDbi::select(hs,
#                                   keys = universe$ENSEMBL,
#                                   columns = c("SYMBOL"),
#                                   keytype = "ENSEMBL")


remove(hs)

#Omit NAs----------------------------------------------------------------------
magnoid_down <- na.omit(magnoid_down)
squamoid_down <- na.omit(squamoid_down)
bronchioid_down <- na.omit(bronchioid_down)
magnoid_up <- na.omit(magnoid_up)
squamoid_up <- na.omit(squamoid_up)
bronchioid_up <- na.omit(bronchioid_up)
universe <- na.omit(universe)


## Pathways enrichment---------------------------------------------------
#Get pathways
pathways <- as.list(reactomeEXTID2PATHID)[universe$ENTREZID] #use genes ID to find Pathways ID
pathways <- as.list(reactomePATHID2EXTID)[unique(unlist(pathways))] #use pathways ID to get sets of genes

# Calculate pathways enrichment
magnoid_downRes <- fora(pathways, genes=magnoid_down$ENTREZID, universe=universe$ENTREZID)
magnoid_downRes <- magnoid_downRes %>% 
  mutate(pathway = as.list(reactomePATHID2NAME)[magnoid_downRes$pathway])

squamoid_downRes <- fora(pathways, genes=squamoid_down$ENTREZID, universe=universe$ENTREZID)
squamoid_downRes <- squamoid_downRes %>% 
  mutate(pathway = as.list(reactomePATHID2NAME)[squamoid_downRes$pathway])

bronchioid_downRes <- fora(pathways, genes=bronchioid_down$ENTREZID, universe=universe$ENTREZID)
bronchioid_downRes <- bronchioid_downRes %>% 
  mutate(pathway = as.list(reactomePATHID2NAME)[bronchioid_downRes$pathway])

magnoid_upRes <- fora(pathways, genes=magnoid_up$ENTREZID, universe=universe$ENTREZID)
magnoid_upRes <- magnoid_upRes %>% 
  mutate(pathway = as.list(reactomePATHID2NAME)[magnoid_upRes$pathway])

squamoid_upRes <- fora(pathways, genes=squamoid_up$ENTREZID, universe=universe$ENTREZID)
squamoid_upRes <- squamoid_upRes %>% 
  mutate(pathway = as.list(reactomePATHID2NAME)[squamoid_upRes$pathway])

bronchioid_upRes <- fora(pathways, genes=bronchioid_up$ENTREZID, universe=universe$ENTREZID)
bronchioid_upRes <- bronchioid_upRes %>% 
  mutate(pathway = as.list(reactomePATHID2NAME)[bronchioid_upRes$pathway])

save(magnoid_downRes, squamoid_downRes, bronchioid_downRes, magnoid_upRes, squamoid_upRes, bronchioid_upRes, file = "Results/03copy_pathway_enrichment.Rdata")

#Check genes in a pathway--------------------------------------------
# path1genes <- as_tibble(bronchioid_downRes$overlapGenes[[1]])
# 
# path1genes <- path1genes %>% 
#   rename("genes" = "value")
# 
# path1genes <- AnnotationDbi::select(hs,
#                                          keys = path1genes$genes,
#                                          columns = c("SYMBOL"),
#                                          keytype = "ENTREZID")

