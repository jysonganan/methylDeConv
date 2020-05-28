## matched 450k and EPIC samples GSE111165
# Additionally, fluorescence activated cell sorting (FACS) was performed to separate cells positive for a neuronal marker,
#resulting in five neuronal positive samples (brain_pos) and 
# twelve neuronal negative ones (brain_neg) with sufficient DNA quantity for analysis with the EPIC array.
# > table(pD.all[,"tissue:ch1"])
# 
# blood   brain_neg   brain_pos brain_whole      buccal      saliva 
# 21          12           5          21          21          21 
# > 
#   > pD.all <- pData(geoMat[[1]])
# > table(pD.all[,"tissue:ch1"])
# 
# blood brain_whole      saliva 
# 12          22          12 

## pD.all[,"subject:ch1"]
## pD.all[,"tissue:ch1"]
library(GEOquery)
library(minfi)
# getGEOSuppFiles("GSE111165")
# untar("GSE111165/GSE111165_RAW.tar", exdir = "GSE111165/idat")
# head(list.files("GSE111165/idat", pattern = "idat"))
# 
# idatFiles <- list.files("GSE111165/idat", pattern = "idat.gz$", full = TRUE)
# sapply(idatFiles, gunzip, overwrite = TRUE)

geoMat <- getGEO("GSE111165")


tissue = "saliva"


if (tissue == "blood"){    ## 6 matched samples subject:101,102,103,104,105,106
  pD.450k <- pData(geoMat[[1]])
  pD.EPIC <- pData(geoMat[[2]])
  common_subject_id <- intersect(pD.450k[pD.450k[,"tissue:ch1"]=="blood","subject:ch1"], 
                                 pD.EPIC[pD.EPIC[,"tissue:ch1"]=="blood","subject:ch1"])
  matched_samples_450k <- rownames(pD.450k[pD.450k[,"tissue:ch1"]=="blood"&pD.450k[,"subject:ch1"]%in%common_subject_id,])
  matched_samples_EPIC <- rownames(pD.EPIC[pD.EPIC[,"tissue:ch1"]=="blood"&pD.EPIC[,"subject:ch1"]%in%common_subject_id,])
}

if (tissue == "saliva"){ ## 6 matched samples subject:101,102,103,104,105,106
  pD.450k <- pData(geoMat[[1]])
  pD.EPIC <- pData(geoMat[[2]])
  common_subject_id <- intersect(pD.450k[pD.450k[,"tissue:ch1"]=="saliva","subject:ch1"], 
                                 pD.EPIC[pD.EPIC[,"tissue:ch1"]=="saliva","subject:ch1"])
  matched_samples_450k <- rownames(pD.450k[pD.450k[,"tissue:ch1"]=="saliva"&pD.450k[,"subject:ch1"]%in%common_subject_id,])
  matched_samples_EPIC <- rownames(pD.EPIC[pD.EPIC[,"tissue:ch1"]=="saliva"&pD.EPIC[,"subject:ch1"]%in%common_subject_id,])
}

if (tissue == "brain_whole"){ ## matched samples 14 vs 6 (there are replicates with same subject id in 450k) ## make brain region
  ## also matched; then matched samples: 6 subject:101,102,103,104,105,106
  pD.450k <- pData(geoMat[[1]])
  pD.EPIC <- pData(geoMat[[2]])
  common_subject_id <- intersect(pD.450k[pD.450k[,"tissue:ch1"]=="brain_whole","subject:ch1"], 
                                 pD.EPIC[pD.EPIC[,"tissue:ch1"]=="brain_whole","subject:ch1"])
  matched_samples_EPIC <- rownames(pD.EPIC[pD.EPIC[,"tissue:ch1"]=="brain_whole"&pD.EPIC[,"subject:ch1"]%in%common_subject_id,])
  matched_samples_450k <- rep(NA, length(matched_samples_EPIC))
  for (i in 1:length(matched_samples_EPIC)){
    subject <- pD.EPIC[matched_samples_EPIC, "subject:ch1"][i]
    brain_region <- as.character(pD.EPIC[matched_samples_EPIC, "characteristics_ch1.5"])[i]
    matched_samples_450k[i] <- rownames(pD.450k[pD.450k[,"subject:ch1"] == subject & substr(as.character(pD.450k[,"characteristics_ch1.5"]),1,18) == substr(brain_region,1,18),])
}}





## GSE111165 
## for 450k   46 samples
rgSet<- read.metharray.exp("GSE111165/idat/450k",force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
colnames(betaMat) <- substr(colnames(betaMat),1,10)
betaMat <- betaMat[,matched_samples_450k]


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

load("Flow450kProbesdefault.RData")
probes_select <- probes_oneVsAllttest

library(EpiDISH)
source("projectCellType.R")
Houseman_res <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:8]))
RPC_res <- epidish(betaMat, as.matrix(compTable[probes_select,3:8]), method = "RPC")$estF
CBS_res <- epidish(betaMat, as.matrix(compTable[probes_select,3:8]), method = "CBS")$estF

load("Flow450kProbePreselect_multiclassGlmnet.RData")
#### the deconv + Glmnet on preselected 1800 probes -> select ~922 probes in glmnet
probes_select <- ProbePreselect_multiclassGlmnet[[1]][-1]
Houseman_res_glmnetpreselect <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:8]))
RPC_res_glmnetpreselect <- epidish(betaMat, as.matrix(compTable[probes_select,3:8]), method = "RPC")$estF
CBS_res_glmnetpreselect <- epidish(betaMat, as.matrix(compTable[probes_select,3:8]), method = "CBS")$estF

probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)
library(dplyr)
multiGlmnet_predProb <- predict(ProbePreselect_multiclassGlmnet[[2]], newdata = t(betaMat[probes,]), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])
rownames(multiGlmnet_predProb) <- colnames(betaMat)


save("Houseman_res","RPC_res","CBS_res","Houseman_res_glmnetpreselect","RPC_res_glmnetpreselect",
     "CBS_res_glmnetpreselect", "multiGlmnet_predProb", file = paste0(tissue,"Match450kRes.RData"))







## for EPIC 101 samples   EPIC-epithelial 
rgSet<- read.metharray.exp("GSE111165/idat",force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
colnames(betaMat) <- substr(colnames(betaMat),1,10)
betaMat <- betaMat[,matched_samples_EPIC]



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
set.seed(3)
betaMat_122126_sub <- cbind(betaMat_122126[,phenotype_122126%in%c("Colon epithelial cells","Lung epithelial cells",
                                                                  "Pancreatic acinar cells","Pancreatic duct cells")],
                            betaMat_122126[,sample(which(phenotype_122126 == "cfDNA"),10,replace = FALSE)])
phenotype_122126_sub <- c(rep("Epithelial",10), rep("cfDNA", 10))
ref_betamatrix <- cbind(ref_betamatrix, betaMat_122126_sub)
ref_phenotype <- c(ref_phenotype, phenotype_122126_sub)
source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)



load("FlowEPIC_EpithelialProbesdefault.RData")
probes_select <- probes_oneVsAllttest
library(EpiDISH)
source("projectCellType.R")
Houseman_res_epicEpithelial <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:10]))
RPC_res_epicEpithelial <- epidish(betaMat, as.matrix(compTable[probes_select,3:10]), method = "RPC")$estF
CBS_res_epicEpithelial <- epidish(betaMat, as.matrix(compTable[probes_select,3:10]), method = "CBS")$estF

load("FlowEPIC_EpithelialProbePreselect_multiclassGlmnet.RData")
#### the deconv + Glmnet on preselected 1800 probes -> select ~922 probes in glmnet
probes_select <- ProbePreselect_multiclassGlmnet[[1]][-1]
Houseman_res_glmnetpreselect_epicEpithelial <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:10]))
RPC_res_glmnetpreselect_epicEpithelial <- epidish(betaMat, as.matrix(compTable[probes_select,3:10]), method = "RPC")$estF
CBS_res_glmnetpreselect_epicEpithelial <- epidish(betaMat, as.matrix(compTable[probes_select,3:10]), method = "CBS")$estF

probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)
library(dplyr)
multiGlmnet_predProb_epicEpithelial <- predict(ProbePreselect_multiclassGlmnet[[2]], newdata = t(betaMat[probes,]), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])
rownames(multiGlmnet_predProb_epicEpithelial) <- colnames(betaMat)

save("Houseman_res_epicEpithelial","RPC_res_epicEpithelial","CBS_res_epicEpithelial","Houseman_res_glmnetpreselect_epicEpithelial","RPC_res_glmnetpreselect_epicEpithelial",
     "CBS_res_glmnetpreselect_epicEpithelial", "multiGlmnet_predProb_epicEpithelial", file = paste0(tissue, "MatchEPICEpithelialRes.RData"))























## for EPIC 101 samples   FlowEPIC
rgSet<- read.metharray.exp("GSE111165/idat",force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
colnames(betaMat) <- substr(colnames(betaMat),1,10)
betaMat <- betaMat[,matched_samples_EPIC]



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
source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)



load("FlowEPICProbesdefault.RData")
probes_select <- probes_oneVsAllttest
library(EpiDISH)
source("projectCellType.R")
Houseman_res_epic <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:8]))
RPC_res_epic <- epidish(betaMat, as.matrix(compTable[probes_select,3:8]), method = "RPC")$estF
CBS_res_epic <- epidish(betaMat, as.matrix(compTable[probes_select,3:8]), method = "CBS")$estF

load("FlowEPICProbePreselect_multiclassGlmnet.RData")
#### the deconv + Glmnet on preselected 1800 probes -> select ~922 probes in glmnet
probes_select <- ProbePreselect_multiclassGlmnet[[1]][-1]
Houseman_res_glmnetpreselect_epic <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:8]))
RPC_res_glmnetpreselect_epic <- epidish(betaMat, as.matrix(compTable[probes_select,3:8]), method = "RPC")$estF
CBS_res_glmnetpreselect_epic <- epidish(betaMat, as.matrix(compTable[probes_select,3:8]), method = "CBS")$estF

probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)
library(dplyr)
multiGlmnet_predProb_epic <- predict(ProbePreselect_multiclassGlmnet[[2]], newdata = t(betaMat[probes,]), type = "prob") %>% 
  mutate('class'=names(.)[apply(., 1, which.max)])
rownames(multiGlmnet_predProb_epic) <- colnames(betaMat)

save("Houseman_res_epic","RPC_res_epic","CBS_res_epic","Houseman_res_glmnetpreselect_epic","RPC_res_glmnetpreselect_epic",
     "CBS_res_glmnetpreselect_epic", "multiGlmnet_predProb_epic", file = paste0(tissue, "MatchEPICRes.RData"))
