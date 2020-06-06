#### Analysis_EPIC_reference with scores
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




##########################################################################################
#### EPIC benchmark with true proportions
####################################

###!!! 12 mixture with known proportions can be used as benchmark
##########################################################################################
load("FlowSorted.Blood.EPIC.RData")
annot <- as.data.frame(colData(FlowSorted.Blood.EPIC))
benchmark <- which(annot$CellType == "MIX")
tmp <- getBeta(preprocessNoob(FlowSorted.Blood.EPIC, dyeMethod = "single"))
benchmark_betamatrix <- tmp[,rownames(annot)[benchmark]]
benchmark_trueprop <- annot[benchmark, c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]



##############################################
## MCP-counter scores, within each cell type, compared across samples
##############################################
source("adpat-MCP-counter.R")
probes <- up_probes_oneVsAllttest_celltype(ref_betamatrix, ref_phenotype, pv = 1e-8, MaxDMRs = 100)

scores_Bcells <- MCP-counter_score_within_celltype(benchmark_betamatrix, probes[["Bcells"]])
scores_CD4T <- MCP-counter_score_within_celltype(benchmark_betamatrix, probes[["CD4T"]])
scores_CD8T <- MCP-counter_score_within_celltype(benchmark_betamatrix, probes[["CD8T"]])
scores_Mono <- MCP-counter_score_within_celltype(benchmark_betamatrix, probes[["Mono"]])
scores_Neu <- MCP-counter_score_within_celltype(benchmark_betamatrix, probes[["Neu"]])
scores_NK <- MCP-counter_score_within_celltype(benchmark_betamatrix, probes[["NK"]])


