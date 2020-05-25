### feature deconvolution performance
data_type = "Flow450k"
library(FlowSorted.Blood.450k)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]



library(minfi)
library(GEOquery)
geoMat <- getGEO("GSE77797")
rgSet <- read.metharray.exp("GSE77797/idat")

pD.all <- pData(geoMat[[1]])
pD <- pD.all[, c("title", "geo_accession", "b cell (%):ch1", "cd4+ t cell (%):ch1", "cd8+ t cell (%):ch1",
                 "granulocyte (%):ch1", "monocyte (%):ch1", "natural killer cell (%):ch1")]
sampleNames(rgSet) <- substr(sampleNames(rgSet), 1, 10)
pD <- pD[sampleNames(rgSet),]
pD <- as(pD, "DataFrame")
pData(rgSet) <- pD
grSet <- preprocessNoob(rgSet, dyeMethod = "single")

benchmark_betaMat <- getBeta(grSet)

pD <- pD.all[, c("title", "geo_accession", "b cell (%):ch1", "cd4+ t cell (%):ch1", "cd8+ t cell (%):ch1",
                       "granulocyte (%):ch1", "monocyte (%):ch1", "natural killer cell (%):ch1")]
pD <- pD[sampleNames(rgSet),]
facs <- pD[,3:8]
colnames(facs) <- c("Bcell", "CD4T", "CD8T","Gran","Mono","NK")
benchmark_trueProp <- facs
benchmark_trueProp <- as.data.frame(benchmark_trueProp)
for (i in 1:6){
  benchmark_trueProp[,i] <- as.numeric(as.character(benchmark_trueProp[,i]))
}
benchmark_trueProp <- benchmark_trueProp/100



##### 1. five default probe selection
set.seed(5)
source("refCompTableProbeSelection.R")
probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both")
probes_oneVsAllLimma <- ref_probe_selection_oneVsAllLimma(ref_betamatrix, ref_phenotype, probeSelect = "both")
probes_pairwiseLimma <-  ref_probe_selection_pairwiseLimma(ref_betamatrix, ref_phenotype)
probes_pairwiseGlmnet <- ref_probe_selection_pairwiseGlmnet_cv(ref_betamatrix, ref_phenotype)
probes_multiclassGlmnet <- ref_probe_selection_multiclassGlmnet_cv(ref_betamatrix, ref_phenotype)
save("probes_oneVsAllttest", "probes_oneVsAllLimma", "probes_pairwiseLimma", 
     "probes_pairwiseGlmnet", "probes_multiclassGlmnet", file = paste0(data_type,"Probesdefault.RData"))





##### 2. change top number of probes oneVsallttest/Limma and pairwise Limma, probeSelect = any/both
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
     file = paste0(data_type, "probesChangeNoBothAny.RData"))





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
save("probes_HighVar_1_600", "probes_HighVar_601_1200", "probes_HighVar_1201_1800", file = paste0(data_type, "ProbesHighVar.RData"))



### 4.  automatic feature seletion: recursive feature elimination  --- Better version of RF varImp
## A Random Forest algorithm is used on each iteration to evaluate the model. 
# The algorithm is configured to explore all possible subsets of the attributes. 
library(caret)
set.seed(5)
control <- rfeControl(functions=rfFuncs, method="cv", number=5)
# run the RFE algorithm
results <- rfe(t(ref_betamatrix), factor(ref_phenotype), sizes=c(1:8), rfeControl=control)
save("results",  file = paste0(data_type, "RfeResults.RData"))
# summarize the results
#print(results)
# list the chosen features
#predictors(results)
# plot the results




#### 5. RF on all probes
set.seed(5)
ref_probe_selection_multiclassRF_cv <- function(ref_betamatrix, ref_phenotype, nCores = 4, reps.resamp = 5, reps.repeats = 3, tune_grid = 1:30){
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

model_multiclassRF_cv <- ref_probe_selection_multiclassRF_cv(ref_betamatrix, ref_phenotype, reps.resamp = 5, 
                                                             tune_grid = c(300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000))
save("model_multiclassRF_cv", file = paste0(data_type, "AllProbes_multiclassRF.RData"))






##### 6. apply different machine learning methods in multiclass classification
### start with pre-selected probes, not all probes e.g. probes_oneVsAllttest, RFE top600, probes_multiclassGlmnet


### 6.1 RF on top 600, 1200, 1800 (one vs all ttest probes)
set.seed(5)
source("refCompTableProbeSelection.R")
probes_oneVsAllttest_100 <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 100)
probes_oneVsAllttest_200 <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 200)
probes_oneVsAllttest_300 <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)


ref_probe_selection_multiclassRF_cv <- function(ref_betamatrix, ref_phenotype, nCores = 4, reps.resamp = 5, reps.repeats = 3, tune_grid = 1:30){
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

model_multiclassRF_cv_1 <- ref_probe_selection_multiclassRF_cv(ref_betamatrix[probes_oneVsAllttest_100,], ref_phenotype, reps.resamp = 5, 
                                                             tune_grid = c(3,5,10,15,18,20,24,26,30,32,34,36,40,45,50,55,60))
model_multiclassRF_cv_2 <- ref_probe_selection_multiclassRF_cv(ref_betamatrix[probes_oneVsAllttest_200,], ref_phenotype, reps.resamp = 5, 
                                                               tune_grid = c(3,5,10,15,18,20,24,26,30,32,34,36,40,45,50,55,60))
model_multiclassRF_cv_3 <- ref_probe_selection_multiclassRF_cv(ref_betamatrix[probes_oneVsAllttest_300,], ref_phenotype, reps.resamp = 5, 
                                                               tune_grid = c(3,5,10,15,18,20,24,26,30,32,34,36,40,45,50,55,60))


save("model_multiclassRF_cv_1","model_multiclassRF_cv_2", "model_multiclassRF_cv_3", 
     file = paste0(data_type, "ProbePreselect_multiclassRF.RData"))







##6.2 ML  on top 1800 (one vs all ttest probes) kNN, svm, glmnet

library(dplyr)
library(caret)
library(doParallel)
library(matrixStats)
library(plyr)
library(recipes)
library(adabag)

set.seed(5)
nCores = 4

source("refCompTableProbeSelection.R")
probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)
ProbePreselect_multiclassGlmnet <- ref_probe_selection_multiclassGlmnet_cv(ref_betamatrix[probes,], ref_phenotype)

mlpred_preselectProbes_benchmark_corr <- function(probes_preselect, benchmark_betaMat, ref_betamatrix, ref_phenotype, benchmark_trueProp){
  ## fit a KNN model on the reference data with the probes selected by oneVsAllttest/onevsAllLimma
  ## obtain the class predicted probabilities of the mixture (e.g. benchmark data)
  
  
  print('start knn!')
  
  knnFit <- train(x = t(ref_betamatrix[probes_preselect,]), y = as.factor(ref_phenotype),method = "knn", 
                  preProc = c("center", "scale"), tuneGrid = expand.grid(k = c(2,3,4,5,6,7,8)),
                  trControl = trainControl(method = "cv", classProbs = TRUE))
  knnpredProb <- predict(knnFit, newdata = t(benchmark_betaMat[probes_preselect,]), type = "prob") %>% 
    mutate('class'=names(.)[apply(., 1, which.max)])
  
  corr <- rep(NA, ncol(benchmark_betaMat))
  for (i in 1:ncol(benchmark_betaMat)){
    corr[i] <-cor(as.numeric(knnpredProb[i,1:ncol(benchmark_trueProp)]),as.numeric(as.character(benchmark_trueProp[i,])),method = "spearman")
  }
  print(mean(corr))
  corr <- rep(NA, ncol(benchmark_trueProp))
  for (i in 1:ncol(benchmark_trueProp)){
    corr[i] <-cor(as.numeric(knnpredProb[,i]),as.numeric(as.character(benchmark_trueProp[,i])),method = "spearman")
  }
  print(corr)
  print('complete knn!')
  
  
  
  
  print('start svmRBF!')
  
  svmFit <- train(x = t(ref_betamatrix[probes_preselect,]), y = as.factor(ref_phenotype),method = "svmRadial", 
                  preProc = c("center", "scale"), tuneLength = 10,
                  trControl = trainControl(method = "cv", classProbs = TRUE))
  svm_predProb <- predict(svmFit, newdata = t(benchmark_betaMat[probes_preselect,]), type = "prob") %>% 
    mutate('class'=names(.)[apply(., 1, which.max)])
  
  corr <- rep(NA, ncol(benchmark_betaMat))
  for (i in 1:ncol(benchmark_betaMat)){
    corr[i] <-cor(as.numeric(svm_predProb[i,1:ncol(benchmark_trueProp)]),as.numeric(as.character(benchmark_trueProp[i,])),method = "spearman")
  }
  print(mean(corr))
  corr <- rep(NA, ncol(benchmark_trueProp))
  for (i in 1:ncol(benchmark_trueProp)){
    corr[i] <-cor(as.numeric(svm_predProb[,i]),as.numeric(as.character(benchmark_trueProp[,i])),method = "spearman")
  }
  print(corr)
  print('complete svmRBF!')

}


mlpred_preselectProbes_benchmark_corr(probes, benchmark_betaMat, ref_betamatrix, ref_phenotype, benchmark_trueProp)
  
save("ProbePreselect_multiclassGlmnet", file = paste0(data_type, "ProbePreselect_multiclassGlmnet.RData"))










#### analysis on results
source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)
probeSelect_deconv_benchmark_corr(probes_oneVsAllttest,benchmark_betaMat,compTable[,3:8],benchmark_trueProp)









































### 4.2 ML (stacking) on top 1800 (one vs all ttest probes)
library(dplyr)
library(caret)
library(glmnet)
library(randomForest)
library(gbm)
library(xgboost)
library(caretEnsemble)
library(doParallel)
library(matrixStats)

set.seed(5)
nCores = 4

library(FlowSorted.Blood.450k)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]
source("refCompTableProbeSelection.R")
probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)

algorithmList <- c('glmnet','rf') #'gbm', 'xgbTree'
control <- trainControl(method="repeatedcv", number=5, repeats=3, savePredictions=TRUE, classProbs=TRUE)
grid.glmnet <- expand.grid(.alpha=seq(0.5,1),.lambda =seq(0,1,by=0.01))
grid.rf <- expand.grid(.mtry = c(3,5,10,15,18,20,24,26,30,32,34,36,40,45,50,55,60))


grid.xgboost <- expand.grid(.nrounds = c(50, 75, 100),
                            .eta = c(0.2, 0.3, 0.4),                
                            .gamma = c(0, 1),
                            .max_depth = c(2, 3, 4),
                            .colsample_bytree = c(0.8),                
                            .subsample = c(1), 
                            .min_child_weight = c(1))

grid.gbm <- expand.grid(n.trees = c(100, 200, 250),
                        interaction.depth = c(1, 4, 6),
                        shrinkage = c(0.05,0.1),
                        n.minobsinnode = 10)

models <- caretList(x = t(ref_betamatrix[probes,]), y = factor(ref_phenotype), trControl=control, metric = "Kappa",
                    methodList=algorithmList, tuneList = list(glmnet = caretModelSpec(method = "glmnet", tuneGrid = grid.glmnet),
                                                              rf = caretModelSpec(method = "rf", tuneGrid = grid.rf)
                                                              #gbm = caretModelSpec(method = "gbm", tuneGrid = grid.gbm)
                                                              #xgbTree = caretModelSpec(method="xgbTree", tuneGrid = grid.xgboost)
                                                              #nnet=caretModelSpec(method="nnet", trace=FALSE, tuneLength=1)
                                                              ))

print("finish caretlist")
# stack using glmnet

stackControl <- trainControl(method="repeatedcv", number=5, repeats=3, savePredictions=TRUE, classProbs=TRUE)

stack.glm <- caretStack(models, tuneLength = 10, method="glm", metric="Accuracy", trControl=stackControl)
### We can also use more sophisticated algorithms to combine predictions in an effort 
## to tease out when best to use the different methods. 
## In this case, we can use the random forest algorithm to combine the predictions

# stack using random forest
stack.rf <- caretStack(models, tuneLength = 10, method="rf", metric="Acurracy", trControl=stackControl)  ## also method = "gbm"


## caretList uses lm,rpart and glm to fit x1 and x2 to y (i.e., y~x1+x2)
## This gives three predictions : plm, prpart, pglm
## Then caretStack creates the model ensemble and provides a global prediction with rf such as y~plm+prpart+pglm
save(stack.glm, stack.rf, file = "Flow450kStackingRes.RData")


### SVM, knn requires preprocessing (scale, center)
### RF and other tree-based methods needn't.












## 4.3  ML  on top 1800 (one vs all ttest probes)
library(dplyr)
library(caret)
library(doParallel)
library(matrixStats)
library(plyr)
library(recipes)
library(adabag)

set.seed(5)
nCores = 4

library(FlowSorted.Blood.450k)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]
source("refCompTableProbeSelection.R")
probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)
load("tmp.RData")




## fit a KNN model on the reference data with the probes selected by oneVsAllttest/onevsAllLimma
## obtain the class predicted probabilities of the mixture (e.g. benchmark data)
print('start knn!')

knnFit <- train(x = t(ref_betamatrix[probes,]), y = as.factor(ref_phenotype),method = "knn", 
                preProc = c("center", "scale"), tuneGrid = expand.grid(k = c(2,3,4,5,6,7,8)),
                trControl = trainControl(method = "cv", classProbs = TRUE))
knnpredProb <- predict(knnFit, newdata = t(betaMat_77797[probes,]), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])

corr <- rep(NA, 18)
for (i in 1:18){
  corr[i] <-cor(as.numeric(knnpredProb[i,1:6]),as.numeric(as.character(facs_77797_prop[i,])),method = "spearman")
}
print(mean(corr))
corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(as.numeric(knnpredProb[,i]),as.numeric(as.character(facs_77797_prop[,i])),method = "spearman")
}
print(corr)
print('complete knn!')




print('start svmRBF!')

svmFit <- train(x = t(ref_betamatrix[probes,]), y = as.factor(ref_phenotype),method = "svmRadial", 
                preProc = c("center", "scale"), tuneLength = 10,
                trControl = trainControl(method = "cv", classProbs = TRUE))
svm_predProb <- predict(svmFit, newdata = t(betaMat_77797[probes,]), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])

corr <- rep(NA, 18)
for (i in 1:18){
  corr[i] <-cor(as.numeric(svm_predProb[i,1:6]),as.numeric(as.character(facs_77797_prop[i,])),method = "spearman")
}
print(mean(corr))
corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(as.numeric(svm_predProb[,i]),as.numeric(as.character(facs_77797_prop[,i])),method = "spearman")
}
print(corr)
print('complete svmRBF!')



print('start Bagged AdaBoost!')

AdaBagFit <- train(x = t(ref_betamatrix[probes,]), y = as.factor(ref_phenotype),method = "AdaBag", 
                tuneLength = 10,
                trControl = trainControl(method = "cv", classProbs = TRUE))
AdaBag_predProb <- predict(AdaBagFit, newdata = t(betaMat_77797[probes,]), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])

corr <- rep(NA, 18)
for (i in 1:18){
  corr[i] <-cor(as.numeric(AdaBag_predProb[i,1:6]),as.numeric(as.character(facs_77797_prop[i,])),method = "spearman")
}
print(mean(corr))
corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(as.numeric(AdaBag_predProb[,i]),as.numeric(as.character(facs_77797_prop[,i])),method = "spearman")
}
print(corr)
print('complete  Bagged AdaBoost!')




print('start bayesglm!')

bayesglmFit <- train(x = t(ref_betamatrix[probes,]), y = as.factor(ref_phenotype),method = "bayesglm", 
                   preProc = c("center", "scale"), tuneLength = 10,
                   trControl = trainControl(method = "cv", classProbs = TRUE))
bayesglm_predProb <- predict(bayesglmFit, newdata = t(betaMat_77797[probes,]), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])

corr <- rep(NA, 18)
for (i in 1:18){
  corr[i] <-cor(as.numeric(bayesglm_predProb[i,1:6]),as.numeric(as.character(facs_77797_prop[i,])),method = "spearman")
}
print(mean(corr))
corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(as.numeric(bayesglm_predProb[,i]),as.numeric(as.character(facs_77797_prop[,i])),method = "spearman")
}
print(corr)
print('complete  bayesglm!')


