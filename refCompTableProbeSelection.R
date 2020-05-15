
ref_phenotype <- {if(!is.null(ncol(CellLines.matrix))){
  Pheno1 <- c(rep("Cancer",ncol(CellLines.matrix)))
  ref_phenotype <- c(as.character(Pheno1), as.character(ref_phenotype))}
  else{
    ref_phenotype <- as.character(ref_phenotype)}
  return(ref_phenotype)}

ref_betamatrix <- cbind(CellLines.matrix, ref_betamatrix)


## tools functions
splitit <- function(x) {split(seq_along(x), x)}

#This function creates the pairs for the pairwise matrices
design.pairs <-
  function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    design
  }





### make compTable

ref_compTable <- function(ref_betamatrix, ref_phenotype){
  require(genefilter)

  ref_phenotype <- as.factor(ref_phenotype)
  
  ffComp <- rowFtests(ref_betamatrix, ref_phenotype)
  prof <- vapply(
    X = splitit(ref_phenotype),
    FUN = function(j) rowMeans2(ref_betamatrix, cols = j),
    FUN.VALUE = numeric(nrow(ref_betamatrix)))
  r <- rowRanges(ref_betamatrix)
  
  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", "high", "range")
  return(compTable)
}



### one versus all t-test: Pipeline default (minfi (estimateCellCounts), FlowSorted.Blood.450k)
ref_probe_selection_oneVsAllttest <- function(ref_betamatrix, ref_phenotype, probeSelect, pv = 1e-8, MaxDMRs = 100){
  require(genefilter)

  ref_phenotype <- as.factor(ref_phenotype)
  
  tIndexes <- splitit(ref_phenotype)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(ref_betamatrix))
    x[i] <- 1
    return(rowttests(ref_betamatrix, factor(x)))
  })
  
  if (probeSelect == "any") {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < pv, ]
      yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
      c(rownames(yAny)[seq(MaxDMRs)])
    })
  } else {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < pv, ]
      yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
      yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
      c(rownames(yUp)[seq_len(MaxDMRs/2)],
        rownames(yDown)[seq_len(MaxDMRs/2)])
    })
  }
  
  trainingProbes <- unique(unlist(probeList))
  # return(list(
  #   #coefEsts = coefEsts,
  #   compTable = compTable,
  #   #sampleMeans = pMeans
  #   trainingProbes = trainingProbes,
  #   tstatList = tstatList))
  return(trainingProbes)
}




ref_probe_selection_oneVsAllLimma <- function(ref_betamatrix, ref_phenotype, probeSelect, FDR = 1e-6, MaxDMRs = 100){
  require(genefilter)
  require(MKmisc)
  tIndexes <- splitit(as.factor(ref_phenotype))
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(ref_betamatrix))
    x[i] <- 1
    return(mod.t.test(ref_betamatrix, group = factor(x)))
  })
  
  if (probeSelect == "any") {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "adj.p.value"] < FDR, ]
      yAny <- y[order(abs(y[, "difference in means"]), decreasing = TRUE), ]
      c(rownames(yAny)[seq(MaxDMRs)])
    })
  } else {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "adj.p.value"] < FDR, ]
      yUp <- y[order(y[, "difference in means"], decreasing = TRUE), ]
      yDown <- y[order(y[, "difference in means"], decreasing = FALSE), ]
      c(rownames(yUp)[seq_len(MaxDMRs/2)],
        rownames(yDown)[seq_len(MaxDMRs/2)])
    })
  }
  
  trainingProbes <- unique(unlist(probeList))
  # return(list(
  #   #coefEsts = coefEsts,
  #   compTable = compTable,
  #   #sampleMeans = pMeans
  #   trainingProbes = trainingProbes,
  #   tstatList = tstatList))
  return(trainingProbes)
}





# design.pairs is the function from MethylCIBERSORT
### pairwise limma (pairwise moderated t-test)
ref_probe_selection_pairwiseLimma <- function(ref_betamatrix, ref_phenotype, FDR = 0.01, deltaBeta = 0.2, MaxDMRs = 100){
  
  require(dplyr)
  require(caret)
  require(glmnet)
  require(foreach)
  #require(NMF)
  require(doParallel)
  require(matrixStats)
  require(limma)
  
  ContrastMatrix <- design.pairs(levels(factor(ref_phenotype)))
  Des <- model.matrix(~0 + ref_phenotype)  #one-hot coding
  colnames(Des) <- rownames(ContrastMatrix)
  Fit <- lmFit(ref_betamatrix, Des)%>%
    contrasts.fit(., ContrastMatrix)%>%
    eBayes(.)
  
  FitList <- list()
  for(i in 1:ncol(ContrastMatrix)) {
    
    FitList[[i]] <- topTable(Fit, coef = i, number = nrow(ref_betamatrix))%>%
      mutate(ID = rownames(.))%>%
      filter(adj.P.Val < FDR)
    
    message(paste0(i, " done"))
  }
  
  ## pairwise delta mean/Beta estimates
  Transformed <- data.frame(t(ref_betamatrix))
  Split <- split(Transformed, ref_phenotype)
  Split <- lapply(Split, function(x) colMedians(data.matrix(x)))
  Split <- do.call(cbind, Split)
  rownames(Split) <- rownames(ref_betamatrix)
  
  dbList <- list()
  message("Getting Delta Beta estimates")
  for(i in 1:ncol(ContrastMatrix)) {
    
    dB <- with(data.frame(Split), eval(parse(text = colnames(ContrastMatrix)[[i]])))
    dB <- data.frame(dB = dB, ID = rownames(Split))
    dbList[[i]] <- dB
    message(paste0 (i, " done"))
  }
  #Filter by thresholds
  dbList <- lapply(dbList, function(x) filter(x, abs(dB) > deltaBeta))
  
  for(i in 1:length(FitList)) {
    A1 <- FitList[[i]]
    A1 <- filter(A1 , ID %in% dbList[[i]]$ID)
    A1 <- A1%>%.[rev(order(.$t)),]
    if(nrow(A1) > MaxDMRs) { A1 <-  A1[1:MaxDMRs,]                   }
    FitList[[i]] <- A1
  }
  
  Nonzeros <- lapply(FitList, function(x) dplyr::select(x,ID))
  Nonzeros <- do.call(rbind, Nonzeros)
  Nonzeros <- filter(Nonzeros, !duplicated(ID))
  select_probes <- Nonzeros$ID
  return(select_probes)
}










ref_probe_selection_pairwiseGlmnet <- function(ref_betamatrix, ref_phenotype, nCores = 4, reps.resamp = 20){
  require(dplyr)
  require(caret)
  require(glmnet)
  require(foreach)
  #require(NMF)
  require(doParallel)
  require(matrixStats)
  
  Features.CVparam<- trainControl(method = "boot632",number = reps.resamp,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)
  if(nCores > 1){
    registerDoParallel(makeCluster(nCores))
    message( "Parallelisation schema set up")}
  
  Pairs <- data.frame(t(combn(unique(ref_phenotype),2)), stringsAsFactors = F)
  ##Pairs <- filter(Pairs, !X1 == X2)
  
  FitList <- list()
  
  for(i in 1:nrow(Pairs)) {
    I1 <- ref_phenotype == Pairs[i,]$X1 | ref_phenotype == Pairs[i,]$X2   
    M1 <- ref_betamatrix[,I1]
    P1 <- as.character(ref_phenotype[I1])
    
    Model <- train(x = t(M1), y = factor(P1), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0.5,1),.lambda = seq(0,0.05,by=0.01)), metric = "Kappa")
    
    Nonzeros <- coef(Model$finalModel, s = Model$bestTune$lambda)
    Nonzeros <- as.matrix(Nonzeros)
    Nonzeros <- data.frame(ID = rownames(Nonzeros), Coef = as.numeric(Nonzeros[,1]))
    Nonzeros <- filter(Nonzeros, !Coef == 0)
    FitList[[i]] <- Nonzeros
    message(paste0("pair",i," done of ", nrow(Pairs)))
  }
  
  Nonzeros <- do.call(rbind, FitList)
  Nonzeros <- filter(Nonzeros, !duplicated(ID))
  
  select_probes <- Nonzeros$ID
  return(select_probes)
}









### multiclass glmnet
ref_probe_selection_multiclassGlmnet_cv <- function(ref_betamatrix, ref_phenotype, nCores = 4, reps.resamp = 10){
  require(dplyr)
  require(caret)
  require(glmnet)
  require(foreach)
  #require(NMF)
  require(doParallel)
  require(matrixStats)
  
  Features.CVparam<- trainControl(method = "cv",number = reps.resamp,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)
  if(nCores > 1){
    registerDoParallel(makeCluster(nCores))
    message( "Parallelisation schema set up")}
  
  Model <- train(x = t(ref_betamatrix), y = factor(ref_phenotype), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0.5,1),.lambda = seq(0,0.05,by=0.01)), metric = "Kappa")
  
  message("Retrieving Nonzero Coefficients")
  Nonzeros <-  coef(Model$finalModel, s = Model$bestTune$lambda)
  Nonzeros <- lapply(Nonzeros, function(x) data.frame(ID = rownames(x), Coef = as.numeric(x[,1])))
  Nonzeros <- lapply(Nonzeros, function(x) filter(x, !Coef == 0))
  Nonzeros <- do.call(rbind, Nonzeros)
  Nonzeros <- filter(Nonzeros, !duplicated(ID))
  
  select_probes <- Nonzeros$ID
  return(select_probes)
}
                     
 
ref_probe_selection_multiclassGlmnet <- function(ref_betamatrix, ref_phenotype, nCores = 4, reps.resamp = 20){
  require(dplyr)
  require(caret)
  require(glmnet)
  require(foreach)
  #require(NMF)
  require(doParallel)
  require(matrixStats)
  
  Features.CVparam<- trainControl(method = "boot632",number = reps.resamp,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)
  if(nCores > 1){
    registerDoParallel(makeCluster(nCores))
    message( "Parallelisation schema set up")}
  
  Model <- train(x = t(ref_betamatrix), y = factor(ref_phenotype), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0.5,1),.lambda = seq(0,0.05,by=0.01)), metric = "Kappa")
  
  message("Retrieving Nonzero Coefficients")
  Nonzeros <-  coef(Model$finalModel, s = Model$bestTune$lambda)
  Nonzeros <- lapply(Nonzeros, function(x) data.frame(ID = rownames(x), Coef = as.numeric(x[,1])))
  Nonzeros <- lapply(Nonzeros, function(x) filter(x, !Coef == 0))
  Nonzeros <- do.call(rbind, Nonzeros)
  Nonzeros <- filter(Nonzeros, !duplicated(ID))
  
  select_probes <- Nonzeros$ID
  #prediction_p <- predict(Model, test, type = "prob")
  return(list(select_probes, Model))
}

                 
