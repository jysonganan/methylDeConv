### 
## Reformat FlowSorted.450k; with preprocessNoob
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
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)
probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both", MaxDMRs = 100)
probes_oneVsAllLimma <- ref_probe_selection_oneVsAllLimma(ref_betamatrix, ref_phenotype, probeSelect = "both")
## the above two methods: intersection 538/600 probes
probes_pairwiseLimma <- ref_probe_selection_pairwiseLimma(ref_betamatrix, ref_phenotype)  #681 probes selected
probes_pairwiseGlmnet <- ref_probe_selection_pairwiseGlmnet(ref_betamatrix, ref_phenotype)
probes_multiclassGlmnet <- ref_probe_selection_multiclassGlmnet(ref_betamatrix, ref_phenotype)



### compare the deconvolution performance of probes
load("/Users/junesong/Desktop/causal inference/Flow450kProbeCompTable.RData")

## on 450k benchmark dataset 77797
library(minfi)
library(GEOquery)
geoMat_77797<- getGEO("GSE77797")
pD.all <- pData(geoMat_77797[[1]])
rgSet_77797 <- read.metharray.exp("GSE77797/idat")

pD <- pD.all[, c("title", "geo_accession", "b cell (%):ch1", "cd4+ t cell (%):ch1", "cd8+ t cell (%):ch1",
                 "granulocyte (%):ch1", "monocyte (%):ch1", "natural killer cell (%):ch1")]
## add phenotype data
sampleNames(rgSet_77797) <- substr(sampleNames(rgSet_77797), 1, 10)
pD <- pD[sampleNames(rgSet_77797),]
pD <- as(pD, "DataFrame")
pData(rgSet_77797) <- pD

grSet_77797 <- preprocessNoob(rgSet_77797, dyeMethod = "single")
betaMat_77797 <- getBeta(grSet_77797)




pD_77797 <- pD.all[, c("title", "geo_accession", "b cell (%):ch1", "cd4+ t cell (%):ch1", "cd8+ t cell (%):ch1",
                       "granulocyte (%):ch1", "monocyte (%):ch1", "natural killer cell (%):ch1")]
pD_77797 <- pD_77797[sampleNames(rgSet_77797),]
facs_77797 <- pD_77797[,3:8]
colnames(facs_77797) <- c("Bcell", "CD4T", "CD8T","Gran","Mono","NK")
#facs_77797_prop <- cbind(facs_77797[,"CD8T"], facs_77797[,"CD4T"], facs_77797[,"NK"],
#                         facs_77797[,"Bcell"],facs_77797[,"Mono"],facs_77797[,"Gran"])
facs_77797_prop <- facs_77797
facs_77797_prop <- as.data.frame(facs_77797_prop)
for (i in 1:6){
  facs_77797_prop[,i] <- as.numeric(as.character(facs_77797_prop[,i]))
}
facs_77797_prop <- facs_77797_prop/100




library(EpiDISH)
## Houseman
source("projectCellType.R")
Houseman_res1 <- projectCellType(betaMat_77797[probes_oneVsAllttest,],as.matrix(compTable[probes_oneVsAllttest,3:8]))
Houseman_res2 <- projectCellType(betaMat_77797[probes_oneVsAllLimma,],as.matrix(compTable[probes_oneVsAllLimma,3:8]))
Houseman_res3 <- projectCellType(betaMat_77797[probes_pairwiseLimma,],as.matrix(compTable[probes_pairwiseLimma,3:8]))
Houseman_res4 <- projectCellType(betaMat_77797[probes_pairwiseGlmnet,],as.matrix(compTable[probes_pairwiseGlmnet,3:8]))
Houseman_res5 <- projectCellType(betaMat_77797[probes_multiclassGlmnet,],as.matrix(compTable[probes_multiclassGlmnet,3:8]))

RPC_res1 <- epidish(betaMat_77797, as.matrix(compTable[probes_oneVsAllttest,3:8]), method = "RPC")$estF
RPC_res2 <- epidish(betaMat_77797, as.matrix(compTable[probes_oneVsAllLimma,3:8]), method = "RPC")$estF
RPC_res3 <- epidish(betaMat_77797, as.matrix(compTable[probes_pairwiseLimma,3:8]), method = "RPC")$estF
RPC_res4 <- epidish(betaMat_77797, as.matrix(compTable[probes_pairwiseGlmnet,3:8]), method = "RPC")$estF
RPC_res5 <- epidish(betaMat_77797, as.matrix(compTable[probes_multiclassGlmnet,3:8]), method = "RPC")$estF

CBS_res1 <- epidish(betaMat_77797, as.matrix(compTable[probes_oneVsAllttest,3:8]), method = "CBS")$estF
CBS_res2 <- epidish(betaMat_77797, as.matrix(compTable[probes_oneVsAllLimma,3:8]), method = "CBS")$estF
CBS_res3 <- epidish(betaMat_77797, as.matrix(compTable[probes_pairwiseLimma,3:8]), method = "CBS")$estF
CBS_res4 <- epidish(betaMat_77797, as.matrix(compTable[probes_pairwiseGlmnet,3:8]), method = "CBS")$estF
CBS_res5 <- epidish(betaMat_77797, as.matrix(compTable[probes_multiclassGlmnet,3:8]), method = "CBS")$estF

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(Houseman_res5[,i],as.numeric(as.character(facs_77797_prop[,i])),method = "spearman")
}
corr



#### prediction score of multiple class elastic net:
length(probes_oneVsAllttest)
length(probes_oneVsAllLimma)
length(probes_pairwiseLimma)
length(probes_pairwiseGlmnet)
length(probes_multiclassGlmnet)

probe12 <- intersect(probes_oneVsAllttest,probes_oneVsAllLimma)
probe123 <- intersect(probe12, probes_pairwiseLimma)
probe1234 <- intersect(probe12, intersect(probes_pairwiseLimma, probes_pairwiseGlmnet))
probe12345 <- intersect(probe1234, probes_multiclassGlmnet)
probe123_all <- unique(c(probes_oneVsAllLimma, probes_oneVsAllttest, probes_pairwiseLimma))
probe12_all <- unique(c(probes_oneVsAllLimma, probes_oneVsAllttest))
res1 <- projectCellType(betaMat_77797[probe12_all,],as.matrix(compTable[probe12_all,3:8]))
res2 <- epidish(betaMat_77797, as.matrix(compTable[probe12_all,3:8]), method = "RPC")$estF
res3 <- epidish(betaMat_77797, as.matrix(compTable[probe12_all,3:8]), method = "CBS")$estF

corr <- rep(NA, 18)
for (i in 1:18){
  corr[i] <-cor(res3[i,],as.numeric(as.character(facs_77797_prop[i,])),method = "spearman")
}
mean(corr)





####
# http://rstudio-pubs-static.s3.amazonaws.com/251240_12a8ecea8e144fada41120ddcf52b116.html#tuning-glmnet-models

# multinomial logistic regression r predicted probabilities

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








