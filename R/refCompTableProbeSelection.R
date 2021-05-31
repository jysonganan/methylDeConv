## tools functions
splitit <- function(x) {split(seq_along(x), x)}

#This function creates the pairs for the pairwise matrices (MethylCIBERSORT)
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




#'Average reference profiles over each cell type
#'
#'Compute the average reference profiles for each cell type, using the reference matrix of purified samples
#'from different cell types (phenotype).
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@return The data maxtirx of average reference profiles (each column corresponds to each cell type),
#'except for the first column (Fstat) and the last three columns (low, high, range) which are the summary statistics.
#'@export


ref_compTable <- function(ref_betamatrix, ref_phenotype){
  ref_phenotype <- as.factor(ref_phenotype)

  ffComp <- genefilter::rowFtests(ref_betamatrix, ref_phenotype)
  prof <- vapply(
    X = splitit(ref_phenotype),
    FUN = function(j) matrixStats::rowMeans2(ref_betamatrix, cols = j),
    FUN.VALUE = numeric(nrow(ref_betamatrix)))
  r <- matrixStats::rowRanges(ref_betamatrix)

  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", "high", "range")
  return(compTable)
}




#'One-vs-All t test feature selection
#'
#'The One-vs-All t test feature selection based on the reference matrix.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@param probeSelect The selection can be "any" or "both". If "any", the function selects top probes regardless of
#'up-regulation or down-regulation. If "both", half of top probes are picked up from the up-regulated probes while
#'the other half of the top probes are picked up from the down-regulated probes.
#'@param pv The p-value threshold with default value as 1e-8
#'@param MaxDMRs The number of probes selected with default value as 100.
#'@return A vector of the selected probes.
#'@export


ref_probe_selection_oneVsAllttest <- function(ref_betamatrix, ref_phenotype, probeSelect, pv = 1e-8, MaxDMRs = 100){

  ref_phenotype <- as.factor(ref_phenotype)

  tIndexes <- splitit(ref_phenotype)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(ref_betamatrix))
    x[i] <- 1
    return(genefilter::rowttests(ref_betamatrix, factor(x)))
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
  return(trainingProbes)
}




#'One-vs-All moderated t test feature selection
#'
#'The moderated One-vs-All t test feature selection based on the reference matrix.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@param probeSelect The selection can be "any" or "both". If "any", the function selects top probes regardless of
#'up-regulation or down-regulation. If "both", half of top probes are picked up from the up-regulated probes while
#'the other half of the top probes are picked up from the down-regulated probes.
#'@param FDR The fdr threshold with default value as 1e-6
#'@param MaxDMRs The number of probes selected with default value as 100.
#'@return A vector of the selected probes.
#'@export

ref_probe_selection_oneVsAllLimma <- function(ref_betamatrix, ref_phenotype, probeSelect, FDR = 1e-6, MaxDMRs = 100){
  tIndexes <- splitit(as.factor(ref_phenotype))
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(ref_betamatrix))
    x[i] <- 1
    return(MKmisc::mod.t.test(ref_betamatrix, group = factor(x)))
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
  return(trainingProbes)
}





#'Pairwise moderated t test feature selection
#'
#'The pairwise moderated t test feature selection based on the reference matrix.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@param deltaBeta The threshold for the absolute mean difference. Default value is 0.2.
#'@param FDR The fdr threshold with default value as 0.01
#'@param MaxDMRs The number of probes selected with default value as 100.
#'@return A vector of the selected probes.
#'@export

ref_probe_selection_pairwiseLimma <- function(ref_betamatrix, ref_phenotype, FDR = 0.01, deltaBeta = 0.2, MaxDMRs = 100){

  require(dplyr)
  require(limma)

  #require(caret)
  #require(glmnet)
  #require(foreach)
  #require(doParallel)

  ContrastMatrix <- design.pairs(levels(factor(ref_phenotype)))
  Des <- stats::model.matrix(~0 + ref_phenotype)  #one-hot coding
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
  Split <- lapply(Split, function(x) matrixStats::colMedians(data.matrix(x)))
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
  dbList <- lapply(dbList, function(x) dplyr::filter(x, abs(dB) > deltaBeta))

  for(i in 1:length(FitList)) {
    A1 <- FitList[[i]]
    A1 <- dplyr::filter(A1 , ID %in% dbList[[i]]$ID)
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






#'Pairwise elastic net modeling feature selection
#'
#'Select the non-zero features from the pairwise elastic net modeling on the reference matrix.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@param nCores Number of parallel cores. Default value is 4.
#'@param reps.resamp The number of folds in k-fold cross-validation. Default value is 5.
#'@param reps.repeats The number of complete sets of folds to compute. Default value is 3.
#'@return A vector of the selected probes.
#'@export


ref_probe_selection_pairwiseGlmnet_cv <- function(ref_betamatrix, ref_phenotype, nCores = 4,
                                                  reps.resamp = 5, reps.repeats = 3){
  require(dplyr)
  require(caret)
  require(glmnet)
  #require(foreach)

  require(matrixStats)

  Features.CVparam<- trainControl(method="repeatedcv", repeats=reps.repeats, number = reps.resamp,verboseIter=TRUE,
                                  returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)
  if(nCores > 1){
    doParallel::registerDoParallel(parallel::makeCluster(nCores))
    message("Parallelisation schema set up")}

  Pairs <- data.frame(t(combn(unique(ref_phenotype),2)), stringsAsFactors = F)

  FitList <- list()

  for(i in 1:nrow(Pairs)) {
    I1 <- ref_phenotype == Pairs[i,]$X1 | ref_phenotype == Pairs[i,]$X2
    M1 <- ref_betamatrix[,I1]
    P1 <- as.character(ref_phenotype[I1])

    Model <- train(x = t(M1), y = factor(P1), trControl = Features.CVparam, method = "glmnet" ,
                   tuneGrid = expand.grid(.alpha=seq(0.1,1, by=0.1),.lambda = seq(0,1,by=0.01)), metric = "Kappa")

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
  select_probes <- as.character(select_probes)
  return(select_probes)
}



#'Multi-class elastic net modeling feature selection
#'
#'Select the non-zero features from the multi-class elastic net modeling on the reference matrix.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@param nCores Number of parallel cores. Default value is 4.
#'@param reps.resamp The number of folds in k-fold cross-validation. Default value is 5.
#'@param reps.repeats The number of complete sets of folds to compute. Default value is 3.
#'@return A vector of the selected probes.
#'@export

ref_probe_selection_multiclassGlmnet_cv <- function(ref_betamatrix, ref_phenotype, nCores = 4,
                                                    reps.resamp = 5, reps.repeats = 3){
  require(dplyr)
  require(caret)
  require(glmnet)
  #require(foreach)

  Features.CVparam<- trainControl(method="repeatedcv", repeats=reps.repeats,number = reps.resamp,
                                  verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)
  if(nCores > 1){
    doParallel::registerDoParallel(parallel::makeCluster(nCores))
    message( "Parallelisation schema set up")}

  Model <- train(x = t(ref_betamatrix), y = factor(ref_phenotype), trControl = Features.CVparam, method = "glmnet" ,
                 tuneGrid = expand.grid(.alpha=seq(0.1,1, by=0.1),.lambda = seq(0,1,by=0.01)), metric = "Kappa")

  message("Retrieving Nonzero Coefficients")
  Nonzeros <-  coef(Model$finalModel, s = Model$bestTune$lambda)
  Nonzeros <- lapply(Nonzeros, function(x) data.frame(ID = rownames(x), Coef = as.numeric(x[,1])))
  Nonzeros <- lapply(Nonzeros, function(x) filter(x, !Coef == 0))
  Nonzeros <- do.call(rbind, Nonzeros)
  Nonzeros <- filter(Nonzeros, !duplicated(ID))

  select_probes <- Nonzeros$ID
  select_probes <- as.character(select_probes)
  return(list(select_probes, Model))
}




#'Multi-class RF modeling feature selection
#'
#'Select features based on multi-class Random Forest modeling on the reference matrix of all probes.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@return A vector of the selected probes.
#'@export

ref_probe_selection_multiclassRF_cv <- function(ref_betamatrix, ref_phenotype, nCores = 4, reps.resamp = 5,
                                                reps.repeats = 3, tune_grid = 1:30){
  require(dplyr)
  require(caret)
  require(randomForest)
  #require(foreach)

  Features.CVparam<- trainControl(method="repeatedcv", repeats=reps.repeats,number = reps.resamp,
                                  verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)
  if(nCores > 1){
    doParallel::registerDoParallel(parallel::makeCluster(nCores))
    message( "Parallelisation schema set up")}

  Model <- train(x = t(ref_betamatrix), y = factor(ref_phenotype), trControl = Features.CVparam, method = "rf" ,
                 tuneGrid = expand.grid(.mtry = tune_grid), metric = "Kappa", importance = TRUE)

  return(Model)
}




#'Two-stage feature selection
#'
#'Select features in two stages: firstly, select top features from one-vs-all t test;
#'secondly, select the features with machine learning modeling.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@param preselect The number of top features per cell type selected from one-vs-all t tests. The default value is 300.
#'@param ml_model The machine learning model for feature selection in the second stage. The default value is "elastic net",
#'which correpsonds to selecting the non-zero features from multi-class elastic net modeling on the reference matrix. Otherwise,
#'if the parameter value is "RF", the model selection is based on the important variables learnt from multi-class Random forest modeling;
#'if the parameter value is "rfe", it selects features based on recursive feature elimination algorithm and
#'a Random Forest algorithm is used on each iteration to evaluate the model.
#'@return Model class.
#'@export

ref_probe_selection_twoStage <- function(ref_betamatrix, ref_phenotype, preselect = 300, ml_model = "elastic net"){
  probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,
                                                                probeSelect = "both", MaxDMRs = preselect)

  if (ml_model == "elastic net"){
    model <- ref_probe_selection_multiclassGlmnet_cv(ref_betamatrix[probes,], ref_phenotype)
  }
  else if (ml_model == "RF"){
    model <- ref_probe_selection_multiclassRF_cv(ref_betamatrix[probes,], ref_phenotype, reps.resamp = 5,
                                                 tune_grid = c(3,5,10,15,18,20,24,26,30,32,34,36,40,45,50,55,60))
  }
  else if (ml_model == "rfe"){
    control <- rfeControl(functions=rfFuncs, method="cv", number=5)
    model <- caret::rfe(t(ref_betamatrix[probes,]), factor(ref_phenotype), sizes=c(1:8), rfeControl=control)
  }
  return(model)
}





#'Feature selection based on variability
#'
#'Select the most variable probes across all samples in the reference matrix
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ranks The top number of most variable probes. Default value is 601:1200,
#'which corresponds the second 600 highest variable probes.
#'@return A vector of the selected probes.
#'@export

ref_probe_selection_HighVar <- function(ref_betamatrix, ranks = 601:1200){
  var_per_probe <- apply(ref_betamatrix,1,var)
  var_per_probe <- as.data.frame(var_per_probe)
  rownames(var_per_probe) <- rownames(ref_betamatrix)
  probes_sorted <- rownames(var_per_probe)[order(-var_per_probe[,1])]
  return(probes_sorted[ranks])
}






#'One-vs-All t test feature selection
#'
#'The One-vs-All t test feature selection based on the reference matrix.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@param pv The p-value threshold with default value as 1e-8
#'@param MaxDMRs The number of probes selected with default value as 100.
#'@return A vector of the selected probes.
#'@export

up_probes_oneVsAllttest_celltype <- function(ref_betamatrix, ref_phenotype, pv = 1e-8, MaxDMRs = 100){
  ref_phenotype <- as.factor(ref_phenotype)
  tIndexes <- splitit(ref_phenotype)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(ref_betamatrix))
    x[i] <- 1
    return(genefilter::rowttests(ref_betamatrix, factor(x)))
  })
  probeList <- lapply(tstatList, function(x) {
    y <- x[x[, "p.value"] < pv, ]
    yUp <- y[order(y[, "dm"], decreasing = FALSE), ]   ## cell type x higher than any other cell types
    c(rownames(yUp)[seq(MaxDMRs)])
  })
  return(probeList)
}
