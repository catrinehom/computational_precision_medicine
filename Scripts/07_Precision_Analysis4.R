# Instructions --------------------------------------------------------------------------------------------------------------------------------------------
# Perform Analysis : FGSEA

library(fgsea)
library(stringr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# load the data from differential expression testing. 
load("Data/lung_differential_expression_results_original.Rdata")

res <- data.frame(res)
res

Gene <- rownames(res)
Gene 

gene_rownames <- as_tibble(Gene) %>% 
  rename( "Gene" = "value")

gene_rownames

res <- gene_rownames %>% 
  cbind(as_tibble(res))

res

# Note: we might be interested in pathways / mechanisms that are altered and not just individual genes.
# fast-preranked gene set enrichment analysis (GSEA)
# The GSEA analysis is performed by:
  
# 1. ranking all genes in the data set based on their correlation to the chosen phenotype.
# 2. identifying the rank positions of all members of the gene set, and
# 3. calculating an enrichment score (ES) that represents the difference between the observed rankings 
# and that which would be expected assuming a random rank distribution.


# Note: Instructions for analysis from  https://sbc.shef.ac.uk/workshops/2019-01-14-rna-seq-r/rna-seq-gene-set-testing.nb.html
# Add the names of the genes as SYMBOLS instead of Ensemble Ids. 

# 1. This is the differential expression results data. # Get the ensemble names, up to 15 characters.

DEres <-  res %>% 
  mutate(Gene = str_sub(Gene, 1, 15))
  
DEres



# 2. Map the Ensemble Gene IDs into a SYMBOL Gene Id. 
ens_to_symbol <- AnnotationDbi::select(org.Hs.eg.db, 
                             key= DEres$Gene, 
                             columns= "SYMBOL",keytype="ENSEMBL")


ens_to_symbol <- as_tibble(ens_to_symbol)
ens_to_symbol


# 3. Join the columns. 
joinres <- inner_join(DEres, ens_to_symbol, by=c("Gene"="ENSEMBL"))
joinres


# 4.Find the mean statistic for each gene. 
# res3 <- joinres %>% 
#   dplyr::select(SYMBOL, stat) %>% 
#   
#   # REMOVE NA values in the genes.
#   na.omit() %>% 
#   distinct() %>% 
#   group_by(SYMBOL) %>% 
#   
#   # find the mean statistic. 
#   summarize(stat=mean(stat))
# 
# res3
# 
# ranks <- deframe(res3)
# head(ranks)
# #barplot(ranks)


# 4. Alternative method of above.... 
gseaInput <- joinres %>% 
  
  # REMOVE NA values in the genes. 
  filter( !is.na(SYMBOL)) %>% 
  
  # rank the genes based on the statistic, using the arrange function from tidyverse. 
  arrange(stat)


# 5. View the ranks as numbers. 
new_ranks <- pull(gseaInput,stat)
new_ranks

as_tibble(new_ranks)
# 6. Add gene names to the ranks
names(new_ranks) <- gseaInput$SYMBOL
new_ranks

head(new_ranks)

barplot(new_ranks) -> ranks_gen
ggsave("Doc/images/gene_ranks.png", ranks_gen , width = 5, height = 5)



# 8. We’re going to use the fgsea package for fast pre-ranked gene set enrichment analysis (GSEA)

# PATHWAYS ------------------------------------------------------------------------------------------------------------------------------------------------
# This gene set is downloaded from the website:

# Hypoxia Gene Set
# pathways_hallmark <- gmtPathways("Data/HALLMARK_HYPOXIA.v7.5.1.gmt")

# MSigDB gene sets as symbols
# pathways_hallmark <- gmtPathways("Data/msigdb.v7.5.1.symbols.gmt.txt")

# Oncogene signature gene sets? 
# pathways_hallmark <- gmtPathways("Data/c6.all.v7.5.1.symbols.gmt.txt")
# pathways_hallmark

# Gene Ontology Gene Sets
pathways_hallmark <- gmtPathways("Data/c5.go.bp.v7.5.1.symbols.gmt.txt")
pathways_hallmark

# Question: Which gene sets should be used? 
# Immunologic signature gene sets?


# View these gene sets and chromosomes. 
#pathways_hallmark %>% 
#  head() %>% 
#  lapply(head)


# fgsea Algorithm -----------------------------------------------------------------------------------------------------------------------------------------
# Run the fgsea program with 1000 permutations 

# kegg pathways::"data/msigdb/c2.cp.kegg.v6.2.symbols.gmt"), 
# go annotations:: "data/msigdb/c5.all.v6.2.symbols.gmt"

fgseaRes <- fgsea(pathways=pathways_hallmark, stats=new_ranks, nperm=10)
fgseaRes 

# Names of each pathway that was tested and the stats from doing the test. 
# Find the enriched pathways

fgseaResTidy <- as_tibble(fgseaRes)

dim(fgseaRes)



# Enrichment Scores ---------------------------------------------------------------------------------------------------------------------------------------

ggplot(data = fgseaResTidy, 
       mapping = aes(x = reorder(pathway, NES), y = NES)) +
       geom_col(aes(fill= padj < 0.05)) +
       coord_flip() +
  
       labs(x="Pathway", 
            y="Normalized Enrichment Score",
            title="Hallmark pathways NES from GSEA") + 
  theme_minimal()+
  theme(
    legend.position = "top",
    axis.text.y = element_text(angle =1, hjust=1, size = 6),
    plot.title=element_text(face="bold", size=20, colour = "darkblue"))


# Notes: https://stephenturner.github.io/deseq-to-fgsea/
# Source:http://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=HALLMARK_HYPOXIA&keywords=
# Source :Site:  https://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb


# Enrichment Plot -----------------------------------------------------------------------------------------------------------------------------------------

# The enrichment plot will show where the genes belonging to a particular gene set
# are towards the top or the bottom of the genelist, and how the enrichment score is calculated 
# across the dataset.

# Here we show the enrichment plot for the pathway with the most positive enrichment score

#plotEnrichment(pathways[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], new_ranks)













