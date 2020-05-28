#################################### 
###### EPIC blood
####################################
## make FlowSorted.Blood.EPIC.RData
# library(ExperimentHub)  
# hub <- ExperimentHub()  
# query(hub, "FlowSorted.Blood.EPIC")  
# FlowSorted.Blood.EPIC <- hub[["EH1136"]]  
# FlowSorted.Blood.EPIC  

#Bcell  CD4T  CD8T  Mono   Neu    NK 
#6     7     6     6     6     6 
data_type = "FlowEPIC"
library(GEOquery)
library(minfi)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
load("FlowSorted.Blood.EPIC.RData")
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.EPIC, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.EPIC))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]
#450 CpGs of six immune cell subtypes: neutrophils, B cells, monocytes, NK cells, CD4+ T cells, and CD 8+ T cells
# Consists of 37 magnetic sorted blood cell references and 12 artificial mixture samples.

#################################### 
###### EPIC Epithelial
####################################
#save("betaMat_122126","phenotype_122126", file = "ref_122126_EPICEpithelial.RData")
data_type = "FlowEPIC_Epithelial"
library(GEOquery)
library(minfi)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
load("FlowSorted.Blood.EPIC.RData")
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.EPIC, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.EPIC))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]

load("ref_122126_EPICEpithelial.RData")

# cfDNA         cfDNA In vitro mix 
# 58                          5 
# **Colon epithelial cells           Cortical neurons 
# 3                          2 
# Hepatocytes               In vitro mix 
# 2                          9 
# Leukocytes      **Lung epithelial cells 
# 1                          3 
# **Pancreatic acinar cells      Pancreatic beta cells 
# 2                          1 
# ** Pancreatic duct cells Vascular endothelial cells 
# 2                          2 
set.seed(3)
betaMat_122126_sub <- cbind(betaMat_122126[,phenotype_122126%in%c("Colon epithelial cells","Lung epithelial cells",
                                                                 "Pancreatic acinar cells","Pancreatic duct cells")],
                            betaMat_122126[,sample(which(phenotype_122126 == "cfDNA"),10,replace = FALSE)])
phenotype_122126_sub <- c(rep("Epithelial",10), rep("cfDNA", 10))

ref_betamatrix <- cbind(ref_betamatrix, betaMat_122126_sub)
ref_phenotype <- c(ref_phenotype, phenotype_122126_sub)






##### 1. five default probe selection
set.seed(2)
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
set.seed(2)
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
set.seed(2)
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
set.seed(2)
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

set.seed(2)
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


mlpred_preselectProbes_benchmark_corr(probes, benchmark_betamatrix, ref_betamatrix, ref_phenotype, benchmark_trueprop)

save("ProbePreselect_multiclassGlmnet", file = paste0(data_type, "ProbePreselect_multiclassGlmnet.RData"))


  
  











#### EPIC benchmark with true proportions
####################################

###!!! 12 mixture with known proportions can be used as benchmark
load("FlowSorted.Blood.EPIC.RData")
annot <- as.data.frame(colData(FlowSorted.Blood.EPIC))
benchmark <- which(annot$CellType == "MIX")
tmp <- getBeta(preprocessNoob(FlowSorted.Blood.EPIC, dyeMethod = "single"))
benchmark_betamatrix <- tmp[,rownames(annot)[benchmark]]
benchmark_trueprop <- annot[benchmark, c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]



### GSE112618 6 samples of known proportions  ## maybe problematic as the sum wasn't 1.
geoMat <- getGEO("GSE112618")
pD.all <- pData(geoMat[[1]])
pD <- cbind(as.numeric(pD.all[,"bcell proportion:ch1"]), as.numeric(pD.all[,"cd4t proportion:ch1"]), 
            as.numeric(pD.all[,"cd8t proportion:ch1"]), as.numeric(pD.all[,"monocytes proportion:ch1"]),
            as.numeric(pD.all[,"neutrophils proportion:ch1"]),as.numeric(pD.all[,"nk proportion:ch1"]))
pD <- as.data.frame(pD)
colnames(pD) <- c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")
rownames(pD) <- rownames(pD.all)

getGEOSuppFiles("GSE112618")
untar("GSE112618/GSE112618_RAW.tar", exdir = "GSE112618/idat")
head(list.files("GSE112618/idat", pattern = "idat"))

idatFiles <- list.files("GSE112618/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

rgSet_112618 <- read.metharray.exp("GSE112618/idat",force = TRUE)

benchmark_betamatrix_1 <- getBeta(preprocessNoob(rgSet_112618, dyeMethod = "single"))
benchmark_trueprop_1 <- pD








#### analysis
source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)
probes <- probes_HighVar_1_600
probes <- probes_HighVar_601_1200
probes <- probes_HighVar_1201_1800
probeSelect_deconv_benchmark_corr(probes,benchmark_betamatrix,compTable[,3:8],benchmark_trueprop)

probeSelect_deconv_benchmark_corr(probes,benchmark_betamatrix_1,compTable[,3:8],benchmark_trueprop_1)


### multiglmnet on preselected 1800 onevsall ttest
probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)
library(dplyr)
multiGlmnet_predProb <- predict(ProbePreselect_multiclassGlmnet[[2]], newdata = t(benchmark_betamatrix[probes,]), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])
rownames(multiGlmnet_predProb) <- colnames(benchmark_betamatrix)

corr <- rep(NA, ncol(benchmark_betamatrix))
for (i in 1:ncol(benchmark_betamatrix)){
  corr[i] <-cor(as.numeric(multiGlmnet_predProb[i,1:ncol(benchmark_trueprop)]),as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
}
print(mean(corr))
corr <- rep(NA, ncol(benchmark_trueprop))
for (i in 1:ncol(benchmark_trueprop)){
  corr[i] <-cor(as.numeric(multiGlmnet_predProb[,i]),as.numeric(as.character(benchmark_trueprop[,i])),method = "spearman")
}
print(corr)







### multiglmnet on all probes
library(dplyr)
multiGlmnet_predProb <- predict(probes_multiclassGlmnet[[2]], newdata = t(benchmark_betamatrix), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])
rownames(multiGlmnet_predProb) <- colnames(benchmark_betamatrix)

corr <- rep(NA, ncol(benchmark_betamatrix))
for (i in 1:ncol(benchmark_betamatrix)){
  corr[i] <-cor(as.numeric(multiGlmnet_predProb[i,1:ncol(benchmark_trueprop)]),as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
}
print(mean(corr))
corr <- rep(NA, ncol(benchmark_trueprop))
for (i in 1:ncol(benchmark_trueprop)){
  corr[i] <-cor(as.numeric(multiGlmnet_predProb[,i]),as.numeric(as.character(benchmark_trueprop[,i])),method = "spearman")
}
print(corr)




### high var, combine T cells
benchmark_trueprop_merge <- cbind(benchmark_trueprop[,1], benchmark_trueprop[,2]+benchmark_trueprop[,3],
                                  benchmark_trueprop[,4], benchmark_trueprop[,5], benchmark_trueprop[,6])
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)
probes <- probes_HighVar_1_600
library(EpiDISH)
Houseman_res <- projectCellType(benchmark_betamatrix[probes,],as.matrix(compTable[probes,3:8]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes,3:8]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes,3:8]), method = "CBS")$estF

Houseman_res_merge <- cbind(Houseman_res[,1], Houseman_res[,2]+Houseman_res[,3],
                            Houseman_res[,4],Houseman_res[,5],Houseman_res[,6])
RPC_res_merge <- cbind(RPC_res[,1], RPC_res[,2]+RPC_res[,3],
                       RPC_res[,4], RPC_res[,5], RPC_res[,6])

CBS_res_merge <- cbind(CBS_res[,1], CBS_res[,2]+CBS_res[,3],
                            CBS_res[,4], CBS_res[,5], CBS_res[,6])

corr <- rep(NA, ncol(benchmark_betamatrix))
for (i in 1:ncol(benchmark_betamatrix)){
  corr[i] <-cor(Houseman_res_merge[i,],as.numeric(as.character(benchmark_trueprop_merge[i,])),method = "spearman")
}
print(mean(corr))

corr <- rep(NA, 5)
for (i in 1:5){
  corr[i] <-cor(Houseman_res_merge[,i],as.numeric(as.character(benchmark_trueprop_merge[,i])),method = "spearman")
}
print(corr)



#### EPIC epithelial reference performance analysis
probes <- probes_HighVar_1201_1800

library(caret)
library(randomForest)
importance <- varImp(model_multiclassRF_cv_3, scale=FALSE) 
probes <- rownames(importance$importance)[1:600]

probes <- ProbePreselect_multiclassGlmnet[[1]][-1]
# probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)
# library(dplyr)
# pred_prop <- predict(ProbePreselect_multiclassGlmnet[[2]], newdata = t(benchmark_betamatrix[probes,]), type = "prob") %>%
#   mutate('class'=names(.)[apply(., 1, which.max)])
# rownames(pred_prop) <- colnames(benchmark_betamatrix)



library(EpiDISH)
Houseman_res <- projectCellType(benchmark_betamatrix[probes,],as.matrix(compTable[probes,3:10]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes,3:10]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes,3:10]), method = "CBS")$estF

pred_prop <-CBS_res


corr <- rep(NA, ncol(benchmark_betamatrix))
for (i in 1:ncol(benchmark_betamatrix)){
  corr[i] <-cor(as.numeric(pred_prop[i,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]),
                as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
}
print(mean(corr))


corr <- rep(NA, ncol(benchmark_trueprop))
for (i in 1:ncol(benchmark_trueprop)){
  corr[i] <-cor(as.numeric(pred_prop[,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")][,i]),
                as.numeric(as.character(benchmark_trueprop[,i])),method = "spearman")
}
print(corr)








#### in addition, EPIC epithelial reference have another reference data: (from other cfDNAs)
### check whether the highest predicted cell type are cfDNAs.

load("ref_122126_EPICEpithelial.RData")

# cfDNA         cfDNA In vitro mix 
# 58                          5 
# **Colon epithelial cells           Cortical neurons 
# 3                          2 
# Hepatocytes               In vitro mix 
# 2                          9 
# Leukocytes      **Lung epithelial cells 
# 1                          3 
# **Pancreatic acinar cells      Pancreatic beta cells 
# 2                          1 
# ** Pancreatic duct cells Vascular endothelial cells 
# 2                          2 
set.seed(3)
betaMat_122126_cfDNA <- betaMat_122126[,which(phenotype_122126 == "cfDNA")]
forref <- colnames(betaMat_122126)[sample(which(phenotype_122126 == "cfDNA"),10,replace = FALSE)]
benchmark_betamatrix <- betaMat_122126_cfDNA[,!colnames(betaMat_122126_cfDNA)%in%forref]
### 48 samples of cfDNA



probes <- probes_HighVar_1_600

library(caret)
library(randomForest)
importance <- varImp(model_multiclassRF_cv_3, scale=FALSE) 
probes <- rownames(importance$importance)[1:600]

probes <- ProbePreselect_multiclassGlmnet[[1]][-1]
# probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)
# library(dplyr)
# pred_prop <- predict(ProbePreselect_multiclassGlmnet[[2]], newdata = t(benchmark_betamatrix[probes,]), type = "prob") %>%
#   mutate('class'=names(.)[apply(., 1, which.max)])
# rownames(pred_prop) <- colnames(benchmark_betamatrix)



library(EpiDISH)
Houseman_res <- projectCellType(benchmark_betamatrix[probes,],as.matrix(compTable[probes,3:10]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes,3:10]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes,3:10]), method = "CBS")$estF

table(apply(Houseman_res,1,which.max))
table(apply(RPC_res,1,which.max))
table(apply(CBS_res,1,which.max))




write.csv(RNAMatrix[rownames(BetaMatrix_average),], "/Users/junesong/Downloads/RNAMatrix_average.csv")
write.csv(dat_Quantile, "/Users/junesong/Desktop/methyl_Kuan/dat_Quantile.csv")
write.csv(resdata1, file = "DESeq2-results-with-normalized_A1vsD1.csv")
write.csv(resdata1, file = "DESeq2-results-with-normalized_A1vsC1.csv")
write.csv(RNAMatrix[rownames(BetaMatrix_average),], "/Users/junesong/Downloads/RNAMatrix_average.csv")
write.csv(BetaMatrix_average,"/Users/junesong/Downloads/BetaMatrix_average.csv")






###### EPIC EWAS

## EPIC data of 150 normal saliva samples (24 female + 126 male) 
#168 samples, 24 from females and 144 from males
library(GEOquery)
library(minfi)
getGEOSuppFiles("GSE111631")
untar("GSE111631/GSE111631_RAW.tar", exdir = "GSE111631/idat")
head(list.files("GSE111631/idat", pattern = "idat"))

idatFiles <- list.files("GSE111631/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)


## GSE116298  GBM

# 3-4 samples per tumor from 12 adult glioblastoma patients (38 samples), 
# and from 3 adult meningioma patients (9 samples) as a homogeneous reference giving a total of 47 samples.
library(GEOquery)
library(minfi)
getGEOSuppFiles("GSE116298")
untar("GSE116298/GSE116298_RAW.tar", exdir = "GSE116298/idat")
head(list.files("GSE116298/idat", pattern = "idat"))

idatFiles <- list.files("GSE116298/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)
