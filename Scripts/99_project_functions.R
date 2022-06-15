################################
# Functions for Wilcox
################################

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