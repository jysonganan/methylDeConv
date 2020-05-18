### feature deconvolution performance

library(FlowSorted.Blood.450k)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]

##### 1. change top number of probes oneVsallttest/Limma and pairwise Limma, probeSelect = any/both

source("refCompTableProbeSelection.R")
probes_oneVsAllttest_100_any <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "any", MaxDMRs = 100)
probes_oneVsAllttest_150 <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both", MaxDMRs = 150)
probes_oneVsAllttest_150_any <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "any", MaxDMRs = 150)
probes_oneVsAllttest_200 <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both", MaxDMRs = 200)
probes_oneVsAllttest_200_any <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "any", MaxDMRs =200)

probes_oneVsAllLimma_100_any <- ref_probe_selection_oneVsAllLimma(ref_betamatrix, ref_phenotype, probeSelect = "any", MaxDMRs = 100)
probes_oneVsAllLimma_150 <- ref_probe_selection_oneVsAllLimma(ref_betamatrix, ref_phenotype, probeSelect = "both", MaxDMRs = 150)
probes_oneVsAllLimma_150_any <- ref_probe_selection_oneVsAllLimma(ref_betamatrix, ref_phenotype, probeSelect = "any", MaxDMRs = 150)
probes_oneVsAllLimma_200 <- ref_probe_selection_oneVsAllLimma(ref_betamatrix, ref_phenotype, probeSelect = "both", MaxDMRs = 200)
probes_oneVsAllLimma_200_any <- ref_probe_selection_oneVsAllLimma(ref_betamatrix, ref_phenotype, probeSelect = "any", MaxDMRs = 200)

## pairwise Limma: probe_select actually correpsond to any.
probes_pairwiseLimma_150 <- ref_probe_selection_pairwiseLimma(ref_betamatrix, ref_phenotype,MaxDMRs = 150) 
probes_pairwiseLimma_200 <- ref_probe_selection_pairwiseLimma(ref_betamatrix, ref_phenotype,MaxDMRs = 200) 
save("probes_oneVsAllttest_100_any", "probes_oneVsAllttest_150", "probes_oneVsAllttest_150_any","probes_oneVsAllttest_200",
     "probes_oneVsAllttest_200_any", "probes_oneVsAllLimma_100_any", "probes_oneVsAllLimma_150", "probes_oneVsAllLimma_150_any",
     "probes_oneVsAllLimma_200","probes_oneVsAllLimma_200_any", "probes_pairwiseLimma_150", "probes_pairwiseLimma_200",
     file = "Flow450kprobesChangeNoBothAny.RData")



##### 2. for pairwise glmnet and multiclass glmnet, change the cross-validation parameters as 5-folds CV (repeatedcv, repeat three times)
#set.seed(5)
source("refCompTableProbeSelection.R")
probes_pairwiseGlmnet_cv <- ref_probe_selection_pairwiseGlmnet_cv(ref_betamatrix, ref_phenotype, reps.resamp = 5)
probes_multiclassGlmnet_cv <- ref_probe_selection_multiclassGlmnet_cv(ref_betamatrix, ref_phenotype, reps.resamp = 5)
save("probes_pairwiseGlmnet_cv", "probes_multiclassGlmnet_cv", file = "Flow450kGlmnetAllprobesCV.RData")



#### 2.1. more tunning? 
#tuneGrid = expand.grid(.alpha=seq(0.1,1, by=0.1),.lambda=seq(0,1,by=0.01)),

ref_probe_selection_multiclassGlmnet_cv <- function(ref_betamatrix, ref_phenotype, nCores = 4, reps.resamp = 10, reps.repeats = 3){
  require(dplyr)
  require(caret)
  require(glmnet)
  require(foreach)
  #require(NMF)
  require(doParallel)
  require(matrixStats)
  
  Features.CVparam<- trainControl(method="repeatedcv", repeats=reps.repeats,number = reps.resamp,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)
  if(nCores > 1){
    registerDoParallel(makeCluster(nCores))
    message( "Parallelisation schema set up")}
  
  Model <- train(x = t(ref_betamatrix), y = factor(ref_phenotype), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=seq(0.1,1, by=0.1),.lambda =seq(0,1,by=0.01)), metric = "Kappa")
  
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

probes_multiclassGlmnet_cv_moregrid <- ref_probe_selection_multiclassGlmnet_cv(ref_betamatrix, ref_phenotype, reps.resamp = 5)
save("probes_multiclassGlmnet_cv_moregrid", file = "Flow450kGlmnetAllprobesMoreGrid.RData")



### 2.2 RF on all probes

ref_probe_selection_multiclassRF_cv <- function(ref_betamatrix, ref_phenotype, nCores = 4, reps.resamp = 10, reps.repeats = 3, tune_grid = 1:30){
  require(dplyr)
  require(caret)
  require(randomForest)
  require(foreach)
  require(doParallel)
  require(matrixStats)
  
  Features.CVparam<- trainControl(method="repeatedcv", repeats=reps.repeats,number = reps.resamp,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)
  if(nCores > 1){
    registerDoParallel(makeCluster(nCores))
    message( "Parallelisation schema set up")}
 
  Model <- train(x = t(ref_betamatrix), y = factor(ref_phenotype), trControl = Features.CVparam, method = "rf" , tuneGrid = expand.grid(.mtry = tune_grid), metric = "Kappa", importance = TRUE)
  
  return(Model)
}

model_multiclassRF_cv <- ref_probe_selection_multiclassRF_cv(ref_betamatrix, ref_phenotype, reps.resamp = 5, tune_grid = 300:600)
save("model_multiclassRF_cv", file = "model_multiclassRF_cv.RData")










### 3.  more probe preselection: most variable probes ; RF varImp

ref_probe_selection_HighVar <- function(ref_betamatrix, ranks = 601:1200){
  var_per_probe <- apply(ref_betamatrix,1,var)
  var_per_probe <- as.data.frame(var_per_probe)
  rownames(var_per_probe) <- rownames(ref_betamatrix)
  probes_sorted <- rownames(var_per_probe)[order(-var_per_probe[,1])]
  return(probes_sorted[ranks])
}

probes_HighVar_1_600<- ref_probe_selection_HighVar(ref_betamatrix, ranks = 1:600)
probes_HighVar_601_1200 <- ref_probe_selection_HighVar(ref_betamatrix)
probes_HighVar_1201_1800<- ref_probe_selection_HighVar(ref_betamatrix, ranks = 1201:1800)
save("probes_HighVar_1_600", "probes_HighVar_601_1200", "probes_HighVar_1201_1800", file = "Flow450kProbesHighVar.RData")




### automatic feature seletion: recursive feature elimination
## A Random Forest algorithm is used on each iteration to evaluate the model. 
# The algorithm is configured to explore all possible subsets of the attributes. 
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(ref_betamatrix, factor(ref_phenotype), sizes=c(1:8), rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))











##### 4. apply different machine learning methods in multiclass classification
### start with pre-selected probes, not all probes e.g. probes_oneVsAllttest

## 4.1 stacking
require(dplyr)
require(caret)
require(glmnet)
require(foreach)
#require(NMF)
require(doParallel)
require(matrixStats)

nCores = 4
probes <- probes_oneVsAllttest


if(nCores > 1){
  registerDoParallel(makeCluster(nCores))
  message( "Parallelisation schema set up")}


grid.xgboost <- expand.grid(.nrounds = c(40, 50, 60),
                            .eta = c(0.2, 0.3, 0.4),                
                            .gamma = c(0, 1),
                            .max_depth = c(2, 3, 4),
                            .colsample_bytree = c(0.8),                
                            .subsample = c(1), 
                            .min_child_weight = c(1))

grid.rf <- expand.grid(.mtry = 3:6)

Features.CVparam<- trainControl(method = "boot632",number = reps.resamp,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)


Model <- train(x = t(ref_betamatrix), y = factor(ref_phenotype), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0.5,1),.lambda = seq(0,0.05,by=0.01)), metric = "Kappa")


