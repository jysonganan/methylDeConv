### adapt the MCP-counter idea
### one versus all t-test, but only choose the up-regulated signatures for each cell type
source("refCompTableProbeSelection.R")
up_probes_oneVsAllttest_celltype <- function(ref_betamatrix, ref_phenotype, pv = 1e-8, MaxDMRs = 100){
  require(genefilter)
  
  ref_phenotype <- as.factor(ref_phenotype)
  
  tIndexes <- splitit(ref_phenotype)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(ref_betamatrix))
    x[i] <- 1
    return(rowttests(ref_betamatrix, factor(x)))
  })
  probeList <- lapply(tstatList, function(x) {
    y <- x[x[, "p.value"] < pv, ]
    yUp <- y[order(y[, "dm"], decreasing = FALSE), ]   ## cell type x higher than any other cell types
    c(rownames(yUp)[seq(MaxDMRs)])
  })
  return(probeList)
}

MCP_counter_score_within_celltype <- function(betamatrix, probes_celltype){
  require(psych)
  mat <- betamatrix[probes_celltype,]
  scores <- apply(mat, 2, function(x){return(geometric.mean(x,na.rm = TRUE))})
  return(scores)
}




#### or use oneVsAllttest (up-regulated) selected CpGs and perform ssGSEA directly on methylation data
ssGSEA_score_within_celltype <- function(betamatrix, probeList){
  require(GSVA)
  require(GSEABase)
  scores <- gsva(betamatrix, probeList, method = "ssgsea", ssgsea.norm = FALSE,parallel.sz = 4, parallel.type = 'SOCK')
  scores = scores - apply(scores,1,min)
  return(scores)
}








#### we can map the signatures of expression data as CpGs (xCell 489 signatures sets for 64 cell types)
#### use the mapped CpGs as signatures and then perform  MCP_counter directly on methylation data??

#### instead of collasping methylation data into gene expression data:
#### we can map the signatures of expression data as CpGs 
#### use the mapped CpGs as signatures and then perform ssGSEA directly on methylation data)




# ### ? how to extract signatures from xCell.data
# 
# 
# load("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/xCell.data.rda")
# signatures = xCell.data$signatures
# genes = xCell.data$genes
# 
# 
# ## build geneset
# geneIds <- geneIds(signatures[[1]]) # any character vector would do
# gs_new <- GeneSet(geneIds)
# setNames(gs_new) <- setName(signatures[[1]])
# 
# # build genesetCollection
# gsc <- GeneSetCollection(gs_new, signatures[[4]])




