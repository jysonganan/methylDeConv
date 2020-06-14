### hierarchical deconvolution


#################################### 
###### EPIC Epithelial (no cfDNA)
####################################
#save("betaMat_122126","phenotype_122126", file = "ref_122126_EPICEpithelial.RData")
data_type = "FlowEPIC_Epithelial_nocfDNA"
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
betaMat_122126_sub <- betaMat_122126[,phenotype_122126%in%c("Colon epithelial cells","Lung epithelial cells",
                                                            "Pancreatic acinar cells","Pancreatic duct cells")]
phenotype_122126_sub <- rep("Epithelial",10)

ref_betamatrix <- cbind(ref_betamatrix, betaMat_122126_sub)
ref_phenotype <- c(ref_phenotype, phenotype_122126_sub)
#################################### 


ref_betamatrix_sub <- ref_betamatrix[,which(ref_phenotype %in% c("CD4T", "CD8T"))]
ref_phenotype_sub <- ref_phenotype[which(ref_phenotype %in% c("CD4T", "CD8T"))]

ref_phenotype[which(ref_phenotype %in% c("CD4T", "CD8T"))] <- "Tcell"


source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)
compTable_sub <- ref_compTable(ref_betamatrix_sub, ref_phenotype_sub)

## step 1 feature selection
set.seed(2)
nCores = 4
probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both")
probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both", MaxDMRs = 300)
ProbePreselect_multiclassGlmnet <- ref_probe_selection_multiclassGlmnet_cv(ref_betamatrix[probes,], ref_phenotype)
probes_select <- ProbePreselect_multiclassGlmnet[[1]][-1]

## step 2 feature selection
set.seed(2)
nCores = 4
probes_oneVsAllttest_sub <- ref_probe_selection_oneVsAllttest(ref_betamatrix_sub, ref_phenotype_sub, probeSelect = "both")
probes_sub <- ref_probe_selection_oneVsAllttest(ref_betamatrix_sub, ref_phenotype_sub, probeSelect = "both", MaxDMRs = 300)
ProbePreselect_multiclassGlmnet_sub <- ref_probe_selection_pairwiseGlmnet_cv(ref_betamatrix_sub[probes_sub,], ref_phenotype_sub)
probes_select_sub <- ProbePreselect_multiclassGlmnet_sub[-1]


###### input benchmark data

library(EpiDISH)
source("projectCellType.R")
## step 1 deconvolution
## onevsall
Houseman_res <- projectCellType(benchmark_betamatrix[probes_oneVsAllttest,],as.matrix(compTable[probes_oneVsAllttest,3:8]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_oneVsAllttest,3:8]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_oneVsAllttest,3:8]), method = "CBS")$estF

## preselect + glmnet
Houseman_res_glmnet <- projectCellType(benchmark_betamatrix[probes_select,],as.matrix(compTable[probes_select,3:8]))
RPC_res_glmnet <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:8]), method = "RPC")$estF
CBS_res_glmnet <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:8]), method = "CBS")$estF



## step 2 deconvolution
## onevsall
Houseman_res_sub <- projectCellType(benchmark_betamatrix[probes_oneVsAllttest_sub,],as.matrix(compTable_sub[probes_oneVsAllttest_sub,3:4]))
RPC_res_sub <- epidish(benchmark_betamatrix, as.matrix(compTable_sub[probes_oneVsAllttest_sub,3:4]), method = "RPC")$estF
CBS_res_sub <- epidish(benchmark_betamatrix, as.matrix(compTable_sub[probes_oneVsAllttest_sub,3:4]), method = "CBS")$estF

## preselect + glmnet
Houseman_res_glmnet_sub <- projectCellType(benchmark_betamatrix[probes_select_sub,],as.matrix(compTable_sub[probes_select_sub,3:4]))
RPC_res_glmnet_sub <- epidish(benchmark_betamatrix, as.matrix(compTable_sub[probes_select_sub,3:4]), method = "RPC")$estF
CBS_res_glmnet_sub <- epidish(benchmark_betamatrix, as.matrix(compTable_sub[probes_select_sub,3:4]), method = "CBS")$estF



### generate the final deconvolution results and check the correlation with true proportions
finalDeconv <- function(step1Deconv, step2Deconv){
  CD4Tres <- step1Deconv[,"Tcell"] * step2Deconv[,"CD4T"]
  CD8Tres <- step1Deconv[,"Tcell"] * step2Deconv[,"CD8T"]
  res <- cbind(step1Deconv[,"Bcell"], CD4Tres, CD8Tres, step1Deconv[,"Epithelial"], 
               step1Deconv[,"Mono"], step1Deconv[,"Neu"], step1Deconv[,"NK"])
 
  
  corr <- rep(NA, ncol(benchmark_betamatrix))
  for (i in 1:ncol(benchmark_betamatrix)){
    corr[i] <-cor(as.numeric(res[i,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]),
                  as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
  }
  print(mean(corr))
  
 
  corr <- rep(NA, ncol(benchmark_trueprop))
  for (i in 1:ncol(benchmark_trueprop)){
    corr[i] <-cor(as.numeric(res[,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")][,i]),
                  as.numeric(as.character(benchmark_trueprop[,i])),method = "spearman")
  }
  print(corr)
  
  }

