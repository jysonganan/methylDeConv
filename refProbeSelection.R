# design.pairs is the function from MethylCIBERSORT
### pairwise limma (pairwise moderated t-test)
ref_probe_selection_pairwiseLimma <- function(ref_betamatrix, ref_phenotype, FDR = 0.01, deltaBeta = 0.2, MaxDMRs = 100){

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
