.isMatrixBacked <- function(object){
  stopifnot(is(object, "SummarizedExperiment"))
  all(vapply(assays(object), is.matrix, logical(1L)))
}

.isMatrixBackedOrStop <- function(object, FUN){
  if (!.isMatrixBacked(object)){
    stop("'", FUN, "()' only supports matrix-backed minfi objects.",
         call. = FALSE)
  }
}

pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50,
                           #compositeCellType = compositeCellType,
                           probeSelect = probeSelect) {
  
## mSet is preprocessQuantile RGset.
  .isMatrixBackedOrStop(mSet)
  splitit <- function(x) {split(seq_along(x), x)}
  
  p <- getBeta(mSet)
  pd <- as.data.frame(colData(mSet))
  if (!is.null(cellTypes)) {
    if (!all(cellTypes %in% pd$CellType))
      stop("elements of argument 'cellTypes' is not part of ",
           "'mSet$CellType'")
    keep <- which(pd$CellType %in% cellTypes)
    pd <- pd[keep,]
    p <- p[,keep]
  }
  # NOTE: Make cell type a factor
  pd$CellType <- factor(pd$CellType, levels = cellTypes)
  ffComp <- rowFtests(p, pd$CellType)
  prof <- vapply(
    X = splitit(pd$CellType),
    FUN = function(j) rowMeans2(p, cols = j),
    FUN.VALUE = numeric(nrow(p)))
  r <- rowRanges(p)
  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2, -1, 0) + ncol(compTable)] <-
    c("low", "high", "range")
  tIndexes <- splitit(pd$CellType)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(p))
    x[i] <- 1
    return(rowttests(p, factor(x)))
  })
  
  
  
  if (probeSelect == "any") {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < 1e-8, ]
      yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
      c(rownames(yAny)[seq(numProbes * 2)])
    })
  } else {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < 1e-8, ]
      yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
      yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
      c(rownames(yUp)[seq_len(numProbes)],
        rownames(yDown)[seq_len(numProbes)])
    })
  }
  
  trainingProbes <- unique(unlist(probeList))
  p <- p[trainingProbes,]
  
  pMeans <- colMeans2(p)
  names(pMeans) <- pd$CellType
  
  # form <- as.formula(
  #   sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse = "+")))
  # phenoDF <- as.data.frame(model.matrix(~ pd$CellType - 1))
  # colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
  # if (ncol(phenoDF) == 2) {
  #   # Two group solution
  #   X <- as.matrix(phenoDF)
  #   coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
  # } else {
  #   # > 2 groups solution
  #   tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
  #   coefEsts <- tmp$coefEsts
  # }
  
  return(list(
    #coefEsts = coefEsts,
    compTable = compTable,
    #sampleMeans = pMeans
    trainingProbes = trainingProbes))
}

