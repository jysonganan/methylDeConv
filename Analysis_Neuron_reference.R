#################################### 
###### 450k Neuron
####################################
data_type = "FlowBrain450k"
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
###### EPIC Neuron
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
