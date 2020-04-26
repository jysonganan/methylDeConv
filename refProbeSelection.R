if(!is.null(ncol(CellLines.matrix))){
  Pheno1 <- c(rep("Cancer",ncol(CellLines.matrix)))
  ref_phenotype <- c(as.character(Pheno1), as.character(ref_phenotype))}
else{
  ref_phenotype <- as.character(ref_phenotype)}

ref_betamatrix <- cbind(CellLines.matrix, ref_betamatrix)



# design.pairs is the function from MethylCIBERSORT
### pairwise limma (pairwise moderated t-test)
ref_probe_selection_pairwiseLimma <- function(ref_betamatrix, ref_phenotype, FDR = 0.01, deltaBeta = 0.2, MaxDMRs = 100){
  
  require(dplyr)
  require(caret)
  require(glmnet)
  require(foreach)
  require(NMF)
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
  select_probes <- rownames(ref_betamatrix) %in% Nonzeros$ID
  return(select_probes)
}




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

                     
                     

                     
                     
 
ref_probe_selection_pairwiseGlmnet <- function(ref_betamatrix, ref_phenotype, nCores = 4, reps.resamp = 20){
  require(caret)
  require(glmnet)
  require(foreach)
  require(NMF)
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
  
  select_probes <- rownames(ref_betamatrix) %in% Nonzeros$ID
  return(select_probes)
}

                     
                     
                     
                     
                     
                     
                     
                     
### multiclass glmnet
ref_probe_selection_multiclassGlmnet <- function(ref_betamatrix, ref_phenotype, nCores = 4, reps.resamp = 20){
  require(caret)
  require(glmnet)
  require(foreach)
  require(NMF)
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

  select_probes <- rownames(ref_betamatrix) %in% Nonzeros$ID
  return(select_probes)
}
