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
load("/Users/junesong/Desktop/causal inference/Flow450kprobesChangeNoBothAny.RData")
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








## using the model fitted in  ref_probe_selection_multiclassGlmnet 
## to obtain the class predicted probabilities of the mixture (e.g. benchmark data)
load("/Users/junesong/Desktop/causal inference/Flow450kProbeMultiPred.RData")  ## default multi glmnet
load("/Users/junesong/Desktop/causal inference/Flow450kGlmnetAllprobesCV.RData") ## pairwise, multiclass glmnet (repeatedcv,5 fold, rep 3) 
### conclusion: Same probes selected in glmnet on whole data, even given different cross-validation methods.
library(caret)
ggplot(probes_multiclassGlmnet[[2]])
ggplot(probes_multiclassGlmnet_cv[[2]])
ggplot(probes_pairwiseGlmnet_cv[[2]])
## a little differences... almost all cross-validation accuracy == 1 even changing regularization parameters.

library(dplyr)
multiGlmnet_predProb <- predict(probes_multiclassGlmnet[[2]], newdata = t(betaMat_77797), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])
rownames(multiGlmnet_predProb) <- colnames(betaMat_77797)

multiGlmnet_cv_predProb <- predict(probes_multiclassGlmnet_cv[[2]], newdata = t(betaMat_77797), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])
rownames(multiGlmnet_cv_predProb) <- colnames(betaMat_77797)

corr <- rep(NA, 18)
for (i in 1:18){
  corr[i] <-cor(as.numeric(multiGlmnet_predProb[i,1:6]),as.numeric(as.character(facs_77797_prop[i,])),method = "spearman")
}
mean(corr)

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(as.numeric(multiGlmnet_predProb[,i]),as.numeric(as.character(facs_77797_prop[,i])),method = "spearman")
}
corr

corr <- rep(NA, 18)
for (i in 1:18){
  corr[i] <-cor(as.numeric(multiGlmnet_cv_predProb[i,1:6]),as.numeric(as.character(facs_77797_prop[i,])),method = "spearman")
}
mean(corr)

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(as.numeric(multiGlmnet_cv_predProb[,i]),as.numeric(as.character(facs_77797_prop[,i])),method = "spearman")
}
corr
#### conclusion: Same prediction performance in glmnet on whole data, even given different cross-validation methods.


##### check with more grids
load("/Users/junesong/Desktop/causal inference/Flow450kGlmnetAllprobesMoreGrid.RData")
length(intersect(probes_multiclassGlmnet_cv_moregrid[[1]],probes_multiclassGlmnet[[1]]))
### it selected more probes than the default, however it contains the previous probes in default (910 in 2495)
multiGlmnet_predProb <- predict(probes_multiclassGlmnet_cv_moregrid[[2]], newdata = t(betaMat_77797), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])
rownames(multiGlmnet_predProb) <- colnames(betaMat_77797)

corr <- rep(NA, 18)
for (i in 1:18){
  corr[i] <-cor(as.numeric(multiGlmnet_predProb[i,1:6]),as.numeric(as.character(facs_77797_prop[i,])),method = "spearman")
}
mean(corr)

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(as.numeric(multiGlmnet_predProb[,i]),as.numeric(as.character(facs_77797_prop[,i])),method = "spearman")
}
corr

#### also better deconvolution performance with better tuning grid
ggplot(probes_multiclassGlmnet_cv_moregrid[[2]])


## pairwise glmnet with better tuning
load("/Users/junesong/Desktop/causal inference/Flow450kGlmnetPairwiseMoreGrid.RData")
length(intersect(probes_pairwiseGlmnet_cv_moregrid,as.character(probes_pairwiseGlmnet)))  ## 9047/17855




### high var probes
load("/Users/junesong/Desktop/causal inference/Flow450kProbesHighVar.RData")
length(intersect(probes_HighVar_1_600, probes_oneVsAllttest))
length(intersect(probes_HighVar_601_1200, probes_oneVsAllttest))
length(intersect(probes_HighVar_1201_1800, probes_oneVsAllttest))




# estimate variable importance
importance <- varImp(model, scale=FALSE) #For most classification models, 
#each predictor will have a separate variable importance for each class 
#(the exceptions are classification trees, bagged trees and boosted trees).
# summarize importance
print(importance)
# plot importance
plot(importance)

varImp(modelFit)

library("gbm")
library(caretEnsemble)
### stacking (bagging/boosting)
# Example of Stacking algorithms
# create submodels
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions=TRUE, classProbs=TRUE)
algorithmList <- c('lda', 'rpart', 'glm', 'knn', 'svmRadial')
set.seed(seed)
models <- caretList(Class~., data=dataset, trControl=control, methodList=algorithmList)
results <- resamples(models)  ## Number of resamples: 30 
summary(results)
dotplot(results)
## When we combine the predictions of different models using stacking, 
## it is desirable that the predictions made by the sub-models have low correlation. 
## If the predictions for the sub-models were highly corrected (>0.75) 
## then they would be making the same or very similar predictions most of the time 
## reducing the benefit of combining the predictions.

# correlation between results
modelCor(results)
splom(results)

# stack using glm
stackControl <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions=TRUE, classProbs=TRUE)
set.seed(seed)
stack.glm <- caretStack(models, method="glm", metric="Accuracy", trControl=stackControl)
print(stack.glm)
predict(stack.glm, newx, newy, type = "prob")
### We can also use more sophisticated algorithms to combine predictions in an effort 
## to tease out when best to use the different methods. 
## In this case, we can use the random forest algorithm to combine the predictions

# stack using random forest
set.seed(seed)
stack.rf <- caretStack(models, method="rf", metric="Accuracy", trControl=stackControl).  ## also method = "gbm"
print(stack.rf)

## caretList uses lm,rpart and glm to fit x1 and x2 to y (i.e., y~x1+x2)
## This gives three predictions : plm, prpart, pglm
## Then caretStack creates the model ensemble and provides a global prediction with rf such as y~plm+prpart+pglm


grid.xgboost <- expand.grid(.nrounds = c(40, 50, 60),
                            .eta = c(0.2, 0.3, 0.4),                
                            .gamma = c(0, 1),
                            .max_depth = c(2, 3, 4),
                            .colsample_bytree = c(0.8),                
                            .subsample = c(1), 
                            .min_child_weight = c(1))

grid.rf <- expand.grid(.mtry = 3:6)

ctrl <- trainControl(method="cv",
                     number=5,
                     returnResamp = "final",
                     savePredictions = "final",
                     classProbs = TRUE,
                     selectionFunction = "oneSE",
                     verboseIter = TRUE,
                     summaryFunction = twoClassSummary)


model_list <- caretList(Class ~.,
                        data = Sonar,
                        trControl = ctrl,
                        tuneList = list(
                          xgbTree = caretModelSpec(method="xgbTree",
                                                   tuneGrid = grid.xgboost),
                          rf = caretModelSpec(method = "rf",
                                              tuneGrid = grid.rf))
)


models_stack <- caretStack(
  model_list,
  tuneLength = 10,
  method ="glmnet",
  metric = "ROC",
  trControl = ctrl
)















library(EpiDISH)
## Houseman
source("projectCellType.R")
Houseman_res1 <- projectCellType(betaMat_77797[probes_oneVsAllttest,],as.matrix(compTable[probes_oneVsAllttest,3:8]))
Houseman_res2 <- projectCellType(betaMat_77797[probes_oneVsAllLimma_200_any,],as.matrix(compTable[probes_oneVsAllLimma_200_any,3:8]))
Houseman_res3 <- projectCellType(betaMat_77797[probes_pairwiseLimma_200,],as.matrix(compTable[probes_pairwiseLimma_200,3:8]))
Houseman_res4 <- projectCellType(betaMat_77797[probes_pairwiseGlmnet,],as.matrix(compTable[probes_pairwiseGlmnet,3:8]))
Houseman_res5 <- projectCellType(betaMat_77797[probes_multiclassGlmnet[[1]],],as.matrix(compTable[probes_multiclassGlmnet[[1]],3:8]))

RPC_res1 <- epidish(betaMat_77797, as.matrix(compTable[probes_oneVsAllttest,3:8]), method = "RPC")$estF
RPC_res2 <- epidish(betaMat_77797, as.matrix(compTable[probes_oneVsAllLimma_200_any,3:8]), method = "RPC")$estF
RPC_res3 <- epidish(betaMat_77797, as.matrix(compTable[probes_pairwiseLimma_200,3:8]), method = "RPC")$estF
RPC_res4 <- epidish(betaMat_77797, as.matrix(compTable[probes_pairwiseGlmnet,3:8]), method = "RPC")$estF
RPC_res5 <- epidish(betaMat_77797, as.matrix(compTable[probes_multiclassGlmnet[[1]],3:8]), method = "RPC")$estF

CBS_res1 <- epidish(betaMat_77797, as.matrix(compTable[probes_oneVsAllttest,3:8]), method = "CBS")$estF
CBS_res2 <- epidish(betaMat_77797, as.matrix(compTable[probes_oneVsAllLimma_200_any,3:8]), method = "CBS")$estF
CBS_res3 <- epidish(betaMat_77797, as.matrix(compTable[probes_pairwiseLimma_200,3:8]), method = "CBS")$estF
CBS_res4 <- epidish(betaMat_77797, as.matrix(compTable[probes_pairwiseGlmnet,3:8]), method = "CBS")$estF
CBS_res5 <- epidish(betaMat_77797, as.matrix(compTable[probes_multiclassGlmnet,3:8]), method = "CBS")$estF


probes <- c(probes_HighVar_1_600,probes_HighVar_601_1200, probes_HighVar_1201_1800)
probes <- probes_HighVar_1201_1800

probes <- probes_multiclassGlmnet_cv_moregrid[[1]][-1]
#probes <- as.character(probes_multiclassGlmnet[[1]])[-1]

#probes <- as.character(probes_pairwiseGlmnet)[-1]
probes <- probes_pairwiseGlmnet_cv_moregrid[-1]
#load("/Users/junesong/Desktop/causal inference/RfeResults600.RData")
load("/Users/junesong/Desktop/causal inference/RfeResult.RData")  ## best 485512
probes <- results$optVariables[1:600]

Houseman_res <- projectCellType(betaMat_77797[probes,],as.matrix(compTable[probes,3:8]))
RPC_res <- epidish(betaMat_77797, as.matrix(compTable[probes,3:8]), method = "RPC")$estF
CBS_res <- epidish(betaMat_77797, as.matrix(compTable[probes,3:8]), method = "CBS")$estF

corr <- rep(NA, 18)
for (i in 1:18){
    corr[i] <-cor(CBS_res[i,],as.numeric(as.character(facs_77797_prop[i,])),method = "spearman")
}
mean(corr)



corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(Houseman_res[,i],as.numeric(as.character(facs_77797_prop[,i])),method = "spearman")
}
corr

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(area1 = 600, area2 = 600,  cross.area = 1, 
                   category = c("oneVsAllttest", "HighVar:1st top"), 
                   lty = "blank", fill = c("blue", "mediumorchid")) # area2 7916 area3 7099

grid.newpage()
draw.pairwise.venn(area1 = 600, area2 = 600,  cross.area = 4, 
                   category = c("oneVsAllttest", "HighVar: 2nd top"), 
                   lty = "blank", fill = c("blue", "mediumorchid")) # area2 7916 area3 7099


grid.newpage()
draw.pairwise.venn(area1 = 600, area2 = 600,  cross.area = 23, 
                   category = c("oneVsAllttest", "HighVar: 3rd top"), 
                   lty = "blank", fill = c("blue", "mediumorchid")) # area2 7916 area3 7099
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




