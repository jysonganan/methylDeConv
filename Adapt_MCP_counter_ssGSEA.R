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
  return(log(scores))
}




#### or use oneVsAllttest (up-regulated) selected CpGs and perform ssGSEA directly on methylation data
ssGSEA_score_within_celltype <- function(betamatrix, probeList){
  require(GSVA)
  require(GSEABase)
  scores <- gsva(betamatrix, probeList,mx.diff=TRUE, method = "ssgsea", ssgsea.norm = FALSE,parallel.sz = 4, parallel.type = 'SOCK')
  scores = scores - apply(scores,1,min)
  return(scores)
}


### implement ssGSEA (in ESTIMATE), instead of using the package GSVA

ssGSEAESTIMATE_score_within_celltype <- function(betamatrix, probeList){
  m <- betamatrix
  cpg.names <- rownames(betamatrix)
  sample.names <- colnames(betamatrix)
  Ns <- length(m[1, ]) # Number of sample
  Ng <- length(m[, 1]) # Number of cpgs
  
  ## Sample rank normalization
  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method="average")
  }
  m <- 10000*m/Ng   
  
  score.matrix <- matrix(0, nrow = length(probeList), ncol= Ns)
  
  for (i in 1:length(probeList)){
    signature.set <- probeList[[i]]
    cpg.overlap <- intersect(signature.set, cpg.names)
    if (length(cpg.overlap) == 0) { 
      score.matrix[i, ] <- rep(NA, Ns)
      next
    } else {
      ES.vector <- vector(length=Ns)
      ### enrichment score computation
      for (j in 1:Ns){
        cpg.list <- order(m[, j], decreasing=TRUE)            
        cpg.set2 <- match(cpg.overlap, cpg.names)
        correl.vector <- m[cpg.list, j]
        
        TAG <- sign(match(cpg.list, cpg.set2, nomatch=0))    # 1 (TAG) & 0 (no.TAG)
        no.TAG <- 1 - TAG 
        N <- length(cpg.list) 
        Nh <- length(cpg.set2) 
        Nm <-  N - Nh 
        correl.vector <- abs(correl.vector)^0.25
        sum.correl  <- sum(correl.vector[TAG == 1])
        P0 <- no.TAG / Nm
        F0 <- cumsum(P0)
        Pn <- TAG * correl.vector / sum.correl
        Fn <- cumsum(Pn)
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        
        if (max.ES > - min.ES) {
          arg.ES <- which.max(RES)
        } else {
          arg.ES <- which.min(RES)
        }
        ES <- sum(RES)
        EnrichmentScore <- list(ES=ES,
                                arg.ES=arg.ES,
                                RES=RES,
                                indicator=TAG)
        ES.vector[j] <- EnrichmentScore$ES
      }
      
      score.matrix[i, ] <- ES.vector
      
      }
      
    
  }

  score.data <- data.frame(score.matrix)
  names(score.data) <- sample.names
  row.names(score.data) <- names(probeList)
  return(score.data)
  
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

