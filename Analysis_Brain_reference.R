########################################################################  
###### 450k brain  new ref (29+31 NeuN+ (Neuron) vs 29+31 NeuN-(Glia))
######################################################################## 
data_type = "Brain450k"
library(FlowSorted.DLPFC.450k)
#library(minfi)
GRset_frontCortex <-  preprocessNoob(FlowSorted.DLPFC.450k, dyeMethod = "single")
CellLines.matrix = NULL
cellTypes = c("NeuN_neg","NeuN_pos")
ref_betamatrix <- getBeta(GRset_frontCortex)
ref_phenotype <- as.data.frame(colData(FlowSorted.DLPFC.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]
#### GSE66351
library(GEOquery)
library(minfi)
# getGEOSuppFiles("GSE66351")
# untar("GSE66351/GSE66351_RAW.tar", exdir = "GSE66351/idat")
# head(list.files("GSE66351/idat", pattern = "idat"))
# 
# idatFiles <- list.files("GSE66351/idat", pattern = "idat.gz$", full = TRUE)
# sapply(idatFiles, gunzip, overwrite = TRUE)

rgSet <- read.metharray.exp("GSE66351/idat",force = TRUE)
betaMat <- getBeta(preprocessNoob(rgSet,dyeMethod = "single"))
geoMat <- getGEO("GSE66351")
pD <- pData(geoMat[[1]])

keep_samples <- rownames(pD[pD[,"cell type:ch1"]%in%c("Glia","Neuron"),])
phenotype <- pD[keep_samples, "cell type:ch1"]
colnames(betaMat) <- substr(colnames(betaMat),1,10)
betamatrix <- betaMat[,keep_samples]

ref_betamatrix <- cbind(ref_betamatrix, betamatrix)
ref_phenotype <- c(ref_phenotype, phenotype)

## combine neuron&NeuN+   glia&NeuN-
for (i in 1:length(ref_phenotype)){
  if (ref_phenotype[i] == "Glia"){
    ref_phenotype[i] <- "NeuN_neg"
  }
  if (ref_phenotype[i] == "Neuron"){
    ref_phenotype[i] <- "NeuN_pos"
  }
}






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


##### 2. glmnet on 1800 preselected probes
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

save("ProbePreselect_multiclassGlmnet", file = paste0(data_type, "ProbePreselect_multiclassGlmnet.RData"))



#################################### 
###### The new 450k brain performance on the prediction of purified neurons 


#GSE98203. 88 samples
#Genomic DNA was isolated from FACS-sorted human brain neuronal nuclei.
####################################
source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)

library(GEOquery)
library(minfi)
# getGEOSuppFiles("GSE98203")
# untar("GSE98203/GSE98203_RAW.tar", exdir = "GSE98203/idat")
# head(list.files("GSE98203/idat", pattern = "idat"))
# 
# idatFiles <- list.files("GSE98203/idat", pattern = "idat.gz$", full = TRUE)
# sapply(idatFiles, gunzip, overwrite = TRUE)

rgSet <- read.metharray.exp("GSE98203/idat",force = TRUE)
betaMat <- getBeta(preprocessNoob(rgSet,dyeMethod = "single"))
geoMat <- getGEO("GSE98203")
pD <- pData(geoMat[[1]])

set.seed(2)
source("refCompTableProbeSelection.R")
probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both")
probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both", MaxDMRs = 300)
probes_preselectGlmnet <- ref_probe_selection_pairwiseGlmnet_cv(ref_betamatrix[probes,], ref_phenotype)

library(EpiDISH)
source("projectCellType.R")
probes_select <- probes_oneVsAllttest
probes_select <- probes_preselectGlmnet[-1]
Houseman_res <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:4]))
RPC_res <- epidish(betaMat, as.matrix(compTable[probes_select,3:4]), method = "RPC")$estF
CBS_res <- epidish(betaMat, as.matrix(compTable[probes_select,3:4]), method = "CBS")$estF



#################################### 
###### EPIC Brain
####################################
data_type = "BrainEPIC"
library(GEOquery)
library(minfi)
geoMat <- getGEO("GSE111165")

pD.all <- pData(geoMat[[2]])
pD_ref <- pD.all[pD.all[,"tissue:ch1"]%in%c("brain_neg","brain_pos"),]

rgSet<- read.metharray.exp("GSE111165/idat",force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
colnames(betaMat) <- substr(colnames(betaMat),1,10)
ref_betamatrix <- betaMat[,rownames(pD_ref)]
ref_phenotype <- c(rep("NeuN_neg",12),rep("NeuN_pos",5))
CellLines.matrix = NULL
cellTypes = c("NeuN_neg","NeuN_pos")










#################################### #################################### 
#################################### ####################################  
#### matched 6 samples of whole brain (deconv on 450k vs deconv on EPIC)
#################################### 
## check Analysis_matched450kEPIC.R
#################################### 
## 1. 450k using "Brain450k" reference
### ... omitted here
rgSet<- read.metharray.exp("GSE111165/idat/450k",force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
colnames(betaMat) <- substr(colnames(betaMat),1,10)
betaMat <- betaMat[,matched_samples_450k]


source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)

set.seed(2)
source("refCompTableProbeSelection.R")
probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both")
probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both", MaxDMRs = 300)
probes_preselectGlmnet <- ref_probe_selection_pairwiseGlmnet_cv(ref_betamatrix[probes,], ref_phenotype)
library(EpiDISH)
source("projectCellType.R")
probes_select <- probes_oneVsAllttest
Houseman_res_Brain450k <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:4]))
RPC_res_Brain450k <- epidish(betaMat, as.matrix(compTable[probes_select,3:4]), method = "RPC")$estF
CBS_res_Brain450k <- epidish(betaMat, as.matrix(compTable[probes_select,3:4]), method = "CBS")$estF

probes_select <- probes_preselectGlmnet[-1]
Houseman_res_glmnetpreselect_Brain450k <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:4]))
RPC_res_glmnetpreselect_Brain450k <- epidish(betaMat, as.matrix(compTable[probes_select,3:4]), method = "RPC")$estF
CBS_res_glmnetpreselect_Brain450k <- epidish(betaMat, as.matrix(compTable[probes_select,3:4]), method = "CBS")$estF
save("Houseman_res_Brain450k","RPC_res_Brain450k","CBS_res_Brain450k","Houseman_res_glmnetpreselect_Brain450k",
     "RPC_res_glmnetpreselect_Brain450k","CBS_res_glmnetpreselect_Brain450k", file = paste0(tissue, "MatchBrain450kRes.RData"))


#################################### 
## 2. EPIC data using "BrainEPIC" reference
#################################### 
### ... omitted here
rgSet<- read.metharray.exp("GSE111165/idat",force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
colnames(betaMat) <- substr(colnames(betaMat),1,10)
betaMat <- betaMat[,matched_samples_EPIC]

source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)


set.seed(2)
source("refCompTableProbeSelection.R")
probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both")
probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both", MaxDMRs = 300)
probes_preselectGlmnet <- ref_probe_selection_pairwiseGlmnet_cv(ref_betamatrix[probes,], ref_phenotype)
library(EpiDISH)
source("projectCellType.R")
probes_select <- probes_oneVsAllttest
Houseman_res_BrainEPIC <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:4]))
RPC_res_BrainEPIC <- epidish(betaMat, as.matrix(compTable[probes_select,3:4]), method = "RPC")$estF
CBS_res_BrainEPIC <- epidish(betaMat, as.matrix(compTable[probes_select,3:4]), method = "CBS")$estF

probes_select <- probes_preselectGlmnet[-1]
Houseman_res_glmnetpreselect_BrainEPIC <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:4]))
RPC_res_glmnetpreselect_BrainEPIC <- epidish(betaMat, as.matrix(compTable[probes_select,3:4]), method = "RPC")$estF
CBS_res_glmnetpreselect_BrainEPIC <- epidish(betaMat, as.matrix(compTable[probes_select,3:4]), method = "CBS")$estF
save("Houseman_res_BrainEPIC","RPC_res_BrainEPIC","CBS_res_BrainEPIC","Houseman_res_glmnetpreselect_BrainEPIC",
     "RPC_res_glmnetpreselect_BrainEPIC","CBS_res_glmnetpreselect_BrainEPIC", file = paste0(tissue, "MatchBrainEPICRes.RData"))


#### check whether different brain regions have differences
