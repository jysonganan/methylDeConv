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
  mat <- betamatrix[probes_celltype,]
  scores <- apply(mat, 2, function(x){return(geoMean(x,na.rm = TRUE))})
  return(scores)
}
