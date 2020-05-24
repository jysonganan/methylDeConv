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


### RFE recursive RF feature selection

#load("/Users/junesong/Desktop/causal inference/RfeResults600.RData")
load("/Users/junesong/Desktop/causal inference/RfeResult.RData")  ## best 485512
probes <- results$optVariables[1:600]


### estimate variable importance
## task 2.2
## RF on all probes ## tune_grid = c(300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000)
# on server
load("model_multiclassRF_cv_smallergrid.RData")

library(dplyr)
multiRF_predProb <- predict(model_multiclassRF_cv_4, newdata = t(betaMat_77797[rownames(importance_4$importance),]), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])
rownames(multiRF_predProb) <- colnames(betaMat_77797)

corr <- rep(NA, 18)
for (i in 1:18){
  corr[i] <-cor(as.numeric(multiRF_predProb[i,1:6]),as.numeric(as.character(facs_77797_prop[i,])),method = "spearman")
}
mean(corr)

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(as.numeric(multiRF_predProb[,i]),as.numeric(as.character(facs_77797_prop[,i])),method = "spearman")
}
corr

importance <- varImp(model_multiclassRF_cv, scale=FALSE) 
#For most classification models, 
#each predictor will have a separate variable importance for each class 
#(the exceptions are classification trees, bagged trees and boosted trees).
# summarize importance
probes <- rownames(importance$importance)[1:600]


















#load("/Users/junesong/Desktop/causal inference/RfeResults600.RData")
load("/Users/junesong/Desktop/causal inference/RfeResult.RData")  ## best 485512
probes <- results$optVariables[1:600]






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

# #
# # smoking EWAS example
# setwd("/Users/junesong/Desktop/causal inference")
# input_methyl <- read.table("GSE50660_matrix_processed.txt", header = T, sep = "\t")
# library(EpiDISH)
# 
# print('Reference matrix set as whole blood by default!')
# dat <- input_methyl[,-1]
# rownames(dat) <- input_methyl[,1]
# cell_Prop <- epidish(beta.m = dat, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
# boxplot(cell_Prop)
# 
# write.table(cell_Prop, "GSE50660_cell_prop_smoking_Epidish.txt", sep = "\t", col.names = T, quote = F)
# 
# data(centDHSbloodDMC.m)
# library(RefFreeEWAS)
# 
# 
# common_probe=intersect(row.names(centDHSbloodDMC.m),row.names(dat))
# cell_Prop_2 = projectMix(dat[common_probe,], centDHSbloodDMC.m[common_probe,])
# 
# ## boxplot of the estimated cell proportions
# par(mfrow=c(1,2))
# boxplot(cell_Prop, main = "Epidish RPC")
# boxplot(cell_Prop_2, main = "Houseman")
# 
# ## 
# plot(cell_Prop[,1],cell_Prop_2[,1])
# abline(1,1, col = "red")
# 
# library(stats)
# for(i in 1:7){
#   tmp <- cbind(cell_Prop[,i], cell_Prop_2[,i])
#   colnames(tmp) <- c("epidishRPC", "houseman")
#   tmp <- as.data.frame(tmp)
#   reg <- lm(houseman~epidishRPC, data = tmp)
#   plot(tmp, main = colnames(cell_Prop_2)[i])
#   abline(reg, col = "blue")
#   legend("topleft", bty="n", legend=paste("R2 is", format(summary(reg)$adj.r.squared)), col = "red", cex = 0.6)
# }

