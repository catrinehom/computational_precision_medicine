library(fgsea)
library(reactome.db)

load("Results/prediction.Rdata")

foraRes <- fora(examplePathways, genes=tail(names(exampleRanks), 200), universe=names(exampleRanks))

##Get Pathways

#0 use the Org.Hs.eg.db package to map gene Symbols to Ensemble to Entrez

#1 use your genes ID to find Pathways ID (reactomeEXTID2PATHID)
genesID2pathwaysID <- as.list(reactomeEXTID2PATHID)[c("1","10")] #TODO: substitute ["1"] with list of gene IDs

#2 use pathways ID to get the relative sets of genes reactomePATHID2EXTID
pathways <- as.list(reactomePATHID2EXTID)[unique(unlist(genesID2pathwaysID))]

#3 change pathways ID with names reactomePATHID2NAME
for (i in 1:length(pathways)){
  names(pathways)[i] <- as.list(reactomePATHID2NAME)[names(pathways)[i]]
}

# exampleRanks <- should I use prediction ? 
#   
# universe <- should I use the whole set of genes? I think so, convert rownames(lung_data) to entrez



foraRes <- fora(pathways, genes=tail(names(exampleRanks), 200), universe=names(exampleRanks))
