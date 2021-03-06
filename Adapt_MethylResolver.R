# MethylResolver(methylMix = MethylMix)
# MethylResolver(methylMix = MethylMix, methylSig = MethylSig, purityModel = RFmodel)

set.seed(5)

## alpha: from 0.5–0.9 in corresponding to fitting a regression to 50–90%, of the cpgs.
## select best alpha for each mixture sample separately
#Least Trimmed Squares Regression (LTS Regression)
MethylResolver <- function(methylMix = NULL, methylSig = MethylSig, betaPrime = TRUE, #outputPath = "./",
                           #outputName = "MethylResolver", #doPar = FALSE, numCores = 1,
                           alpha = seq(0.5,0.9,by = 0.05), absolute = TRUE, purityModel = RFmodel) {
    i <- NULL 
  
    if(is.null(methylMix)){
    cat("Please provide a methylation mixture file")
  } else{
      cat("Beginning LTS Deconvolution For This Mixture...\n")
      
      #Format input matrices
      methylSig <- as.data.frame(methylSig)
      methylMix <- data.matrix(varhandle::unfactor(methylMix))
      overlappingCpGs = intersect(rownames(methylSig),rownames(methylMix))
      methylSig = methylSig[overlappingCpGs,]
      methylMix = methylMix[overlappingCpGs,]
      
      #If methylation mixture is supplied as Beta values, convert to Beta'
      if(betaPrime==FALSE){
        methylMix = 1 - methylMix
      }
      
      ltsModel = foreach::foreach(i=1:ncol(methylMix), .combine = 'rbind',.packages = c("robustbase","Metrics")) %do% {
        #the actual model
        regressionFormula = as.formula(paste0("methylMix[,i] ~ ",paste(colnames(methylSig),sep="",collapse=" + ")))
        #check if a specific alpha value is supplied, otherwise use a grid search
        if(length(alpha)>1){
          #get the best alpha value based on RMSE between original and reconstructed mixture
          alphaRMSEs = c()
          for(alphaVal in alpha){
            deconvoluteSample <- robustbase::ltsReg(regressionFormula, data = methylSig, alpha = alphaVal)
            #get the optimal cpgs that are used
            bestCpGs = deconvoluteSample$best
            #get the coefficients (the cell type percentages)
            deconvoluteSample <- deconvoluteSample$coefficients[2:length(deconvoluteSample$coefficients)]
            #change any negative values to 0
            deconvoluteSample[which(deconvoluteSample<0)]<-0
            #if all coefficients are zero, set them all to the same value
            if(sum(deconvoluteSample)==0){
              deconvoluteSample[1:length(deconvoluteSample)] = rep((1/length(deconvoluteSample)),length(deconvoluteSample))
            }
            deconvoluteSample <- deconvoluteSample/sum(deconvoluteSample)
            #evaluate rmse between original and reconstructed mixture
            mHat = deconvoluteSample%*%t(data.matrix(methylSig))
            rmse2 <- Metrics::rmse(methylMix[,i],mHat)
            alphaRMSEs = c(alphaRMSEs,rmse2)
          }
          #choose the best empirical alpha value
          alphaBest = alpha[which(alphaRMSEs == min(alphaRMSEs))]
          if (length(alphaBest) > 1){alphaBest = alphaBest[1]}
        }else{
          alphaBest = alpha
        }
        
        ##
        print(i)
        print('alphaBest selected')
        print(alphaRMSEs)
        print(alphaBest)
            
        #the actual model
        regressionFormula = as.formula(paste0("methylMix[,i] ~ ",paste(colnames(methylSig),sep="",collapse=" + ")))
        deconvoluteSample <- robustbase::ltsReg(regressionFormula, data = methylSig, alpha = alphaBest)
        
        #get the optimal cpgs that are used
        bestCpGs = deconvoluteSample$best
        #get the coefficients (the cell type percentages)
        deconvoluteSample <- deconvoluteSample$coefficients[2:length(deconvoluteSample$coefficients)]
        #change any negative values to 0
        deconvoluteSample[which(deconvoluteSample<0)]<-0
        #if all coefficients are zero, set them all to the same value
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
      
      
      #See if we should predict absolute fractions
      if(absolute == TRUE){
        #Predict sample purity and calculate absolute fractions
        ignoreMetrics = which(colnames(ltsModel) %in% c("RMSE1","R1","RMSE2","R2"))
        if(class(purityModel)[1] == "randomForest.formula"){
          rfModel = purityModel
        }
        else{
          cat("Provided purity model is not correct. Please provide a random forest model trained to
              predict purity...")
        }
        purityPrediction = predict(rfModel,ltsModel)
        absoluteFractions = ltsModel[,-ignoreMetrics]*(1-purityPrediction)
        colnames(absoluteFractions) = paste0("abs_",colnames(absoluteFractions))
        absoluteFractions$Purity = purityPrediction
        ltsModel = cbind(ltsModel,absoluteFractions)
      }
      
            #write to file
      #write.table(ltsModel,file=paste0(outputPath,outputName,".txt"),quote = F, sep="\t")
      
      cat("\nCompleted LTS Deconvolution For This Mixture...\n")
      return(list(ltsModel, bestCpGs))
      
   }
}
     
