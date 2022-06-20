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

magnoid <- append(magnoid_down, magnoid_up)
squamoid <- append(squamoid_down, squamoid_up)
bronchioid <- append(bronchioid_down, bronchioid_up)

universe <- as_tibble(rownames(lung_data))

remove(magnoid_down, magnoid_up, squamoid_down, squamoid_up, bronchioid_down, bronchioid_up, signatures_down, signatures_up, lung_data)

## Convert ENSEMBL to ENTREZID--------------------------------------------------

# Convert lists to tibbles
magnoid <- as_tibble(magnoid)
magnoid <- magnoid %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))

squamoid <- as_tibble(squamoid)
squamoid <- squamoid %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))

bronchioid <- as_tibble(bronchioid)
bronchioid <- bronchioid %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))

universe <- universe %>% 
  rename("genes" = "value") %>% 
  mutate(genes = str_sub(genes, 1, 15))
  
hs <- org.Hs.eg.db #Load conversion dataset

# Convert ENSEMBL to ENTREZ
magnoid <- AnnotationDbi::select(hs, 
                                 keys = magnoid$genes,
                                 columns = c("ENTREZID"),
                                 keytype = "ENSEMBL")

squamoid <- AnnotationDbi::select(hs,
                                  keys = squamoid$genes,
                                  columns = c("ENTREZID"),
                                  keytype = "ENSEMBL")

bronchioid <- AnnotationDbi::select(hs,
                                    keys = bronchioid$genes,
                                    columns = c("ENTREZID"),
                                    keytype = "ENSEMBL")

universe <- AnnotationDbi::select(hs,
                                  keys = universe$genes,
                                  columns = c("ENTREZID"),
                                  keytype = "ENSEMBL")
remove(hs)

#Omit NAs
magnoid <- na.omit(magnoid)
squamoid <- na.omit(squamoid)
bronchioid <- na.omit(bronchioid)
universe <- na.omit(universe)


## Pathways enrichment---------------------------------------------------
#Get pathways
pathways <- as.list(reactomeEXTID2PATHID)[universe$ENTREZID] #use genes ID to find Pathways ID
pathways <- as.list(reactomePATHID2EXTID)[unique(unlist(pathways))] #use pathways ID to get sets of genes

# Calculate pathways enrichment
magnoidRes <- fora(pathways, genes=magnoid$ENTREZID, universe=universe$ENTREZID)
magnoidRes <- magnoidRes %>% 
  mutate(pathway = as.list(reactomePATHID2NAME)[magnoidRes$pathway])

squamoidRes <- fora(pathways, genes=squamoid$ENTREZID, universe=universe$ENTREZID)
squamoidRes <- squamoidRes %>% 
  mutate(pathway = as.list(reactomePATHID2NAME)[squamoidRes$pathway])

bronchioidRes <- fora(pathways, genes=bronchioid$ENTREZID, universe=universe$ENTREZID)
bronchioidRes <- bronchioidRes %>% 
  mutate(pathway = as.list(reactomePATHID2NAME)[bronchioidRes$pathway])

save(magnoidRes, squamoidRes, bronchioidRes, file = "Results/03_pathway_enrichment.Rdata")



