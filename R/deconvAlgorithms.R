
projectCellType <- function(Y, coefCellType, contrastCellType = NULL, nonnegative = TRUE, lessThanOne = FALSE){
  library(quadprog)
  if (is.null(contrastCellType)) {
    Xmat <- coefCellType
  } else {
    Xmat <- tcrossprod(coefCellType, contrastCellType)
  }

  nCol <- dim(Xmat)[2]
  if (nCol == 2) {
    Dmat <- crossprod(Xmat)
    mixCoef <- t(apply(Y, 2, function(x) solve(Dmat, crossprod(Xmat, x))))
    colnames(mixCoef) <- colnames(Xmat)
    return(mixCoef)
  } else {
    nSubj <- dim(Y)[2]

    mixCoef <- matrix(0, nSubj, nCol)
    rownames(mixCoef) <- colnames(Y)
    colnames(mixCoef) <- colnames(Xmat)

    if (nonnegative) {
      if (lessThanOne) {
        Amat <- cbind(rep(-1, nCol), diag(nCol))
        b0vec <- c(-1, rep(0, nCol))
      } else {
        Amat <- diag(nCol)
        b0vec <- rep(0, nCol)
      }
      for (i in seq_len(nSubj)) {
        obs <- which(!is.na(Y[,i]))
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve.QP(
          Dmat = Dmat,
          dvec = crossprod(Xmat[obs,], Y[obs,i]),
          Amat = Amat,
          bvec = b0vec)$sol
      }
    } else {
      for (i in seq_len(nSubj)) {
        obs <- which(!is.na(Y[,i]))
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
      }
    }
    mixCoef
  }
}

#'Houseman algorithm
#'
#'Houseman deconvolution algorithm
#'@param betamatrix Beta value matrix of methylation array for mixture samples.
#'@param compTable Average reference profiles over each cell type.
#'@param probes_select Probes selected.
#'@return The matrix of deconvolution results with rows as mixture samples and columns as cell types.
#'@export

Houseman_project <- function(betamatrix, compTable, probes_select){
  probes <- intersect(rownames(betamatrix), probes_select)
  res <- projectCellType(betamatrix[probes,], as.matrix(compTable[probes,]))
  return(res)
}


#'RPC algorithm
#'
#'RPC deconvolution algorithm
#'@param betamatrix Beta value matrix of methylation array for mixture samples.
#'@param compTable Average reference profiles over each cell type.
#'@param probes_select Probes selected.
#'@return The matrix of deconvolution results with rows as mixture samples and columns as cell types.
#'@export

RPC <- function(betamatrix, compTable, probes_select){
  probes <- intersect(rownames(betamatrix), probes_select)
  res <- EpiDISH::epidish(betamatrix, as.matrix(compTable[probes,]), method = "RPC")$estF
  return(res)
}


#'CBS algorithm
#'
#'CBS deconvolution algorithm
#'@param betamatrix Beta value matrix of methylation array for mixture samples.
#'@param compTable Average reference profiles over each cell type.
#'@param probes_select Probes selected.
#'@return The matrix of deconvolution results with rows as mixture samples and columns as cell types.
#'@export

CBS <- function(betamatrix, compTable, probes_select){
  probes <- intersect(rownames(betamatrix), probes_select)
  res <- EpiDISH::epidish(betamatrix, as.matrix(compTable[probes,]), method = "CBS")$estF
  return(res)
}



#'MethylResolver algorithm
#'
#'MethylResolver deconvolution algorithm
#'@param methylMix Beta value matrix of methylation array for mixture samples.
#'@param methylSig Average reference profiles over each cell type.
#'@param alpha Tuning parameters. Default range 0.5–0.9, which corresponds to fitting a trimmed least square regression
#'to 50–90%, of the cpgs. Get the best alpha value based on RMSE between original and reconstructed mixture, and select
#'best alpha for each mixture sample separately
#'@param probes_select Probes selected.
#'@return The matrix of deconvolution results with rows as mixture samples and columns as cell types.
#'@export

MethylResolver <- function(methylMix, methylSig,
                           alpha = seq(0.5, 0.9, by = 0.05)) {
  i <- NULL
  methylSig <- as.data.frame(methylSig)
  methylMix <- data.matrix(varhandle::unfactor(methylMix))
  overlappingCpGs = intersect(rownames(methylSig),rownames(methylMix))
  methylSig = methylSig[overlappingCpGs,]
  methylMix = methylMix[overlappingCpGs,]


  ltsModel = foreach::foreach(i=1:ncol(methylMix), .combine = 'rbind',.packages = c("robustbase","Metrics")) %do% {
    regressionFormula = as.formula(paste0("methylMix[,i] ~ ",paste(colnames(methylSig),sep="",collapse=" + ")))

    if(length(alpha)>1){
      alphaRMSEs = c()
      for(alphaVal in alpha){
        deconvoluteSample <- robustbase::ltsReg(regressionFormula, data = methylSig, alpha = alphaVal)
        bestCpGs = deconvoluteSample$best
        deconvoluteSample <- deconvoluteSample$coefficients[2:length(deconvoluteSample$coefficients)]
        deconvoluteSample[which(deconvoluteSample<0)]<-0
        if(sum(deconvoluteSample)==0){
          deconvoluteSample[1:length(deconvoluteSample)] = rep((1/length(deconvoluteSample)),length(deconvoluteSample))
        }
        deconvoluteSample <- deconvoluteSample/sum(deconvoluteSample)
        mHat = deconvoluteSample%*%t(data.matrix(methylSig))
        rmse2 <- Metrics::rmse(methylMix[,i],mHat)
        alphaRMSEs = c(alphaRMSEs,rmse2)
      }
      alphaBest = alpha[which(alphaRMSEs == min(alphaRMSEs))]
      if (length(alphaBest) > 1){alphaBest = alphaBest[1]}
    }else{
      alphaBest = alpha
    }


    print(i)
    print('alphaBest selected')
    print(alphaRMSEs)
    print(alphaBest)

    #the actual model
    regressionFormula = as.formula(paste0("methylMix[,i] ~ ",paste(colnames(methylSig),sep="",collapse=" + ")))
    deconvoluteSample <- robustbase::ltsReg(regressionFormula, data = methylSig, alpha = alphaBest)

    #get the optimal cpgs that are used
    bestCpGs = deconvoluteSample$best
    deconvoluteSample <- deconvoluteSample$coefficients[2:length(deconvoluteSample$coefficients)]
    deconvoluteSample[which(deconvoluteSample<0)]<-0
    if(sum(deconvoluteSample)==0){
      deconvoluteSample[1:length(deconvoluteSample)] = rep((1/length(deconvoluteSample)),length(deconvoluteSample))
    }
    deconvoluteSample <- deconvoluteSample/sum(deconvoluteSample)

    #metrics for deconvolution accuracy
    mHat = deconvoluteSample%*%t(data.matrix(methylSig))
    mHat2 = mHat[bestCpGs]
    pearson1 = cor(methylMix[bestCpGs,i],as.vector(mHat2))
    rmse1 <- Metrics::rmse(methylMix[bestCpGs,i],mHat2)
    pearson2 = cor(methylMix[,i],as.vector(mHat))
    rmse2 <- Metrics::rmse(methylMix[,i],mHat)

    #return each sample deconvolution
    deconvoluteOne = data.frame(matrix(c(rmse1,pearson1,rmse2,pearson2,deconvoluteSample),nrow = 1))
    deconvoluteOne
  }

  print('completed')
  #add row and column names
  rownames(ltsModel) <- colnames(methylMix)
  colnames(ltsModel) <- c("RMSE1","R1","RMSE2","R2",colnames(methylSig))


  ignoreMetrics = which(colnames(ltsModel) %in% c("RMSE1","R1","RMSE2","R2"))
  Fractions = ltsModel[,-ignoreMetrics]
  ltsModel = cbind(ltsModel,Fractions)

  cat("\nCompleted LTS Deconvolution For This Mixture...\n")
  return(list(ltsModel, bestCpGs))
}











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

#'Enrichment score-based algorithm
#'
#'CBS deconvolution algorithm
#'@param betamatrix Beta value matrix of methylation array for mixture samples.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@param method The enrichment score-based deconvolution algorithm such as "MCP-counter", "ssGSEA" and "ESTIMATE",
#'which are adapted from gene-expression based algorithm. Default value: "MCP-counter".
#'@param pv The p-value threshold with default value as 1e-8
#'@param MaxDMRs The number of probes selected with default value as 100.
#'@return The matrix of deconvolution results with rows as mixture samples and columns as cell types.
#'@export


enrichment_score <- function(betamatrix, ref_betamatrix, ref_phenotype, method = "MCP-counter", pv = 1e-8, MaxDMRs = 100){
  probes_select_up <- up_probes_oneVsAllttest_celltype(ref_betamatrix, ref_phenotype, pv = pv, MaxDMRs = MaxDMRs)

  if (method == "MCP-counter"){
    mat <- betamatrix[probes_select_up,]
    scores <- apply(mat, 2, function(x){return(mean(x,na.rm = TRUE))})
    return(log(scores))
  }
  else if (method == "ssGSEA"){
    scores <- GSVA::gsva(betamatrix, probes_select_up, mx.diff=TRUE, method = "ssgsea",
                   ssgsea.norm = FALSE,parallel.sz = 4, parallel.type = 'SOCK')
    scores = scores - apply(scores,1,min)
    return(scores)
  }
  else if (method == "ESTIMATE"){
    scores <- ssGSEAESTIMATE_score_within_celltype(betamatrix, probes_select_up)
    return(scores)
  }
}



