load("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/xCell.data.rda")
library(GSVA)
library(GSEABase)
library(pracma)
library(stats)
library(MASS)


## suppose expr is already  # Transform the expression to rank: expr <- apply(expr, 2, rank)

xCellAnalysis <- function(expr, signatures=NULL, genes=NULL, spill=NULL, rnaseq=TRUE, file.name = NULL, scale=TRUE,
                          alpha = 0.5, save.raw = FALSE, parallel.sz = 4, parallel.type = 'SOCK',
                          cell.types.use = NULL) {
  if (is.null(signatures))
    signatures = xCell.data$signatures
  if (is.null(genes))
    genes = xCell.data$genes
  if (is.null(spill)) {
    if (rnaseq==TRUE) {
      spill = xCell.data$spill
    } else {
      spill = xCell.data$spill.array
    }
    # 10808 genes and 489 signatures
  }

  # Caulcate average ssGSEA scores for cell types
  if (is.null(file.name) || save.raw==FALSE) {
    fn <- NULL
  } else {
    fn <- paste0(file.name,'_RAW.txt')
  }

  
  
  if (!is.null(cell.types.use)) {
    A = intersect(cell.types.use,rownames(spill$K))
    if (length(A)<length(cell.types.use)) {
      return ('ERROR - not all cell types listed are available')
    }
  }




  scores <- rawEnrichmentAnalysis(expr,signatures,genes,fn, parallel.sz = parallel.sz, parallel.type = 'SOCK')

  # Transform scores from raw to percentages
  scores.transformed <- transformScores(scores, spill$fv, scale)

  # Adjust scores using the spill over compensation matrix
  if (is.null(file.name)) {
    fn <- NULL
  } else {
    fn <- file.name
  }

  if (is.null(cell.types.use)) {
    scores.adjusted <- spillOver(scores.transformed, spill$K, alpha,fn )
    scores.adjusted = microenvironmentScores(scores.adjusted)
  } else {
    scores.adjusted <- spillOver(scores.transformed[cell.types.use,], spill$K, alpha,fn )
  }
  return(scores.adjusted)
}












rawEnrichmentAnalysis <- function(expr, signatures, genes, file.name = NULL, parallel.sz = 4, parallel.type = 'SOCK') {

  # Reduce the expression dataset to contain only the required genes
  shared.genes <- intersect(rownames(expr), genes)
  #########################################
  # intersect 20621 genes with 10808 genes.
  ########################################
  print(paste("Num. of genes:", length(shared.genes)))
  expr <- expr[shared.genes, ]
  if (dim(expr)[1] < 5000) {
    print(paste("ERROR: not enough genes"))
    return - 1
  }

  

  # Run ssGSEA analysis for the ranked gene expression dataset
  scores <- gsva(expr, signatures, method = "ssgsea",
                       ssgsea.norm = FALSE,parallel.sz = parallel.sz,parallel.type = parallel.type)

  scores = scores - apply(scores,1,min)

  # Combine signatures for same cell types
  cell_types <- unlist(strsplit(rownames(scores), "%"))
  cell_types <- cell_types[seq(1, length(cell_types), 3)]
  agg <- aggregate(scores ~ cell_types, FUN = mean)
  rownames(agg) <- agg[, 1]
  scores <- agg[, -1]

  # Save raw scores
  if (!is.null(file.name)) {
    write.table(scores, file = file.name, sep = "\t",
                col.names = NA, quote = FALSE)
  }
  scores
}














transformScores <- function(scores, fit.vals, scale=TRUE,
                            fn = NULL) {
  rows <- rownames(scores)[rownames(scores) %in% rownames(fit.vals)]
  tscores <- scores[rows, ]
  minX <- apply(tscores, 1, min)
  A <- rownames(tscores)
  tscores <- (as.matrix(tscores) - minX)/5000
  tscores[tscores < 0] <- 0
  if (scale==FALSE) {
    fit.vals[A,3] = 1
  }

  tscores <- (tscores^fit.vals[A,2])/(fit.vals[A,3]*2)

  if (!is.null(fn)) {
    write.table(format(tscores, digits = 4), file = fn, sep = "\t",
                col.names = NA, quote = FALSE)
  }
  return(tscores)
}








spillOver <- function(transformedScores, K, alpha = 0.5, file.name = NULL) {
  K <- K * alpha
  diag(K) <- 1
  rows <- rownames(transformedScores)[rownames(transformedScores) %in%
                                        rownames(K)]
  scores <- apply(transformedScores[rows, ], 2, function(x) lsqlincon(K[rows,rows],
                                                                              x, lb = 0))

  scores[scores<0]=0
  rownames(scores) <- rows
  if (!is.null(file.name)) {
    scores = round(scores*10000)
    scores = scores/10000
    write.table(scores, file = file.name, sep = "\t",
                col.names = NA, quote = FALSE)
  }
  return(scores)
}








microenvironmentScores <- function(adjustedScores) {
  ImmuneScore = apply(adjustedScores[c('B-cells','CD4+ T-cells','CD8+ T-cells','DC','Eosinophils','Macrophages','Monocytes','Mast cells','Neutrophils','NK cells'),],2,sum)/1.5
  StromaScore = apply(adjustedScores[c('Adipocytes','Endothelial cells','Fibroblasts'),],2,sum)/2
  MicroenvironmentScore = ImmuneScore+StromaScore
  adjustedScores = rbind(adjustedScores,ImmuneScore,StromaScore,MicroenvironmentScore)
}






