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
source("Adpat_MCP_counter_ssGSEA.R")
probes <- up_probes_oneVsAllttest_celltype(ref_betamatrix, ref_phenotype, pv = 1e-8, MaxDMRs = 100)

scores_Bcell <- MCP_counter_score_within_celltype(benchmark_betamatrix, probes[["Bcell"]])
scores_CD4T <- MCP_counter_score_within_celltype(benchmark_betamatrix, probes[["CD4T"]])
scores_CD8T <- MCP_counter_score_within_celltype(benchmark_betamatrix, probes[["CD8T"]])
scores_Mono <- MCP_counter_score_within_celltype(benchmark_betamatrix, probes[["Mono"]])
scores_Neu <- MCP_counter_score_within_celltype(benchmark_betamatrix, probes[["Neu"]])
scores_NK <- MCP_counter_score_within_celltype(benchmark_betamatrix, probes[["NK"]])

cor(scores_Bcell, benchmark_trueprop[,"Bcell"], method = "spearman")
cor(scores_CD4T, benchmark_trueprop[,"CD4T"], method = "spearman")
cor(scores_CD8T, benchmark_trueprop[,"CD8T"], method = "spearman")
cor(scores_Mono, benchmark_trueprop[,"Mono"], method = "spearman")
cor(scores_Neu, benchmark_trueprop[,"Neu"], method = "spearman")
cor(scores_NK, benchmark_trueprop[,"NK"], method = "spearman")



scores1 <- ssGSEA_score_within_celltype(benchmark_betamatrix, probes)
# benchmark_betamatrix_rank <- apply(benchmark_betamatrix,2,rank)   ## rank or not rank . same
# scores1 <- ssGSEA_score_within_celltype(benchmark_betamatrix_rank, probes)

cor(scores1[1,], benchmark_trueprop[,"Bcell"], method = "spearman")
cor(scores1[2,], benchmark_trueprop[,"CD4T"], method = "spearman")
cor(scores1[3,], benchmark_trueprop[,"CD8T"], method = "spearman")
cor(scores1[4,], benchmark_trueprop[,"Mono"], method = "spearman")
cor(scores1[5,], benchmark_trueprop[,"Neu"], method = "spearman")
cor(scores1[6,], benchmark_trueprop[,"NK"], method = "spearman")






##############################################
## MCP-counter scores, within each cell type, compared across samples (probes selected based on 450k data)
##############################################
library(FlowSorted.Blood.450k)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]


source("Adpat_MCP_counter_ssGSEA.R")
probes <- up_probes_oneVsAllttest_celltype(ref_betamatrix, ref_phenotype, pv = 1e-8, MaxDMRs = 100)

scores_Bcell <- MCP_counter_score_within_celltype(benchmark_betamatrix, intersect(rownames(benchmark_betamatrix),probes[["Bcell"]]))
scores_CD4T <- MCP_counter_score_within_celltype(benchmark_betamatrix, intersect(rownames(benchmark_betamatrix),probes[["CD4T"]]))
scores_CD8T <- MCP_counter_score_within_celltype(benchmark_betamatrix, intersect(rownames(benchmark_betamatrix),probes[["CD8T"]]))
scores_Mono <- MCP_counter_score_within_celltype(benchmark_betamatrix, intersect(rownames(benchmark_betamatrix),probes[["Mono"]]))
scores_Neu <- MCP_counter_score_within_celltype(benchmark_betamatrix, intersect(rownames(benchmark_betamatrix),probes[["Gran"]]))
scores_NK <- MCP_counter_score_within_celltype(benchmark_betamatrix, intersect(rownames(benchmark_betamatrix),probes[["NK"]]))

cor(scores_Bcell, benchmark_trueprop[,"Bcell"], method = "spearman")
cor(scores_CD4T, benchmark_trueprop[,"CD4T"], method = "spearman")
cor(scores_CD8T, benchmark_trueprop[,"CD8T"], method = "spearman")
cor(scores_Mono, benchmark_trueprop[,"Mono"], method = "spearman")
cor(scores_Neu, benchmark_trueprop[,"Neu"], method = "spearman")
cor(scores_NK, benchmark_trueprop[,"NK"], method = "spearman")




scores1 <- ssGSEA_score_within_celltype(benchmark_betamatrix, probes)
# benchmark_betamatrix_rank <- apply(benchmark_betamatrix,2,rank)   ## rank or not rank . same
# scores1 <- ssGSEA_score_within_celltype(benchmark_betamatrix_rank, probes)

cor(scores1[1,], benchmark_trueprop[,"Bcell"], method = "spearman")
cor(scores1[2,], benchmark_trueprop[,"CD4T"], method = "spearman")
cor(scores1[3,], benchmark_trueprop[,"CD8T"], method = "spearman")
cor(scores1[5,], benchmark_trueprop[,"Mono"], method = "spearman")
cor(scores1[4,], benchmark_trueprop[,"Neu"], method = "spearman")
cor(scores1[6,], benchmark_trueprop[,"NK"], method = "spearman")

