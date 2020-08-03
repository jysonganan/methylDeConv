## simple mixture

purified_datasets <- as.matrix(compTable[,3:9])



library(MCMCpack)

## 1. no epithelial
set.seed(3)
proportions_sim <- rdirichlet(600, c(1,1,1,1,1,1))
proportions_sim <- cbind(proportions_sim, 0)
proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]


mixture_sim_mat <- matrix(NA, nrow(purified_datasets), 600)
for (i in 1:nrow(purified_datasets)){
  mixture_sim_mat[i,] <- purified_datasets[i,] %*% t(proportions_sim)
}

rownames(mixture_sim_mat) <- rownames(purified_datasets)
#mixture_sim_mat <- purified_datasets %*% t(proportions_sim)
colnames(proportions_sim) <- colnames(purified_datasets) 




## 2. low epithelial (0.1 - 0.2)
set.seed(4)
proportions_sim <- rdirichlet(600, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(600, 0.1, 0.2)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]

mixture_sim_mat <- matrix(NA, nrow(purified_datasets), 600)
for (i in 1:nrow(purified_datasets)){
  mixture_sim_mat[i,] <- purified_datasets[i,] %*% t(proportions_sim)
}

rownames(mixture_sim_mat) <- rownames(purified_datasets)
#mixture_sim_mat <- purified_datasets %*% t(proportions_sim)
colnames(proportions_sim) <- colnames(purified_datasets) 




## 3. 
set.seed(2)
proportions_sim <- rdirichlet(600, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(600, 0.2, 0.5)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]

mixture_sim_mat <- matrix(NA, nrow(purified_datasets), 600)
for (i in 1:nrow(purified_datasets)){
  mixture_sim_mat[i,] <- purified_datasets[i,] %*% t(proportions_sim)
}

rownames(mixture_sim_mat) <- rownames(purified_datasets)
#mixture_sim_mat <- purified_datasets %*% t(proportions_sim)
colnames(proportions_sim) <- colnames(purified_datasets) 



# 4
set.seed(5)
proportions_sim <- rdirichlet(600, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(600, 0.5, 0.8)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]

mixture_sim_mat <- matrix(NA, nrow(purified_datasets), 600)
for (i in 1:nrow(purified_datasets)){
  mixture_sim_mat[i,] <- purified_datasets[i,] %*% t(proportions_sim)
}

rownames(mixture_sim_mat) <- rownames(purified_datasets)
#mixture_sim_mat <- purified_datasets %*% t(proportions_sim)
colnames(proportions_sim) <- colnames(purified_datasets) 


# 5
set.seed(7)
proportions_sim <- rdirichlet(600, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(600, 0.8, 0.9)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]

mixture_sim_mat <- matrix(NA, nrow(purified_datasets), 600)
for (i in 1:nrow(purified_datasets)){
  mixture_sim_mat[i,] <- purified_datasets[i,] %*% t(proportions_sim)
}

rownames(mixture_sim_mat) <- rownames(purified_datasets)
#mixture_sim_mat <- purified_datasets %*% t(proportions_sim)
colnames(proportions_sim) <- colnames(purified_datasets) 





##### analysis
# 1. using EPIC reference
load("FlowEPICProbesdefault.RData")
load("FlowEPICProbePreselect_multiclassGlmnet.RData")
source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)

probes <- probes_oneVsAllttest
benchmark_betamatrix <- mixture_sim_mat


probes_select <- intersect(probes, rownames(benchmark_betamatrix))  ## 600
library(EpiDISH)
source("projectCellType.R")

Houseman_res <- projectCellType(benchmark_betamatrix[probes_select,],as.matrix(compTable[probes_select,3:8]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:8]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:8]), method = "CBS")$estF

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(proportions_sim[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], CBS_res[i,], method = "spearman")
}


corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(proportions_sim[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i],Houseman_res[,i], method = "spearman")
}









#2 use EPIC + Epithelial reference
load("FlowEPIC_Epithelial_nocfDNAProbesdefault.RData")
load("FlowEPIC_Epithelial_nocfDNAProbePreselect_multiclassGlmnet.RData")
source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)

probes <- probes_oneVsAllttest
#probes <- ProbePreselect_multiclassGlmnet[[1]][-1]
benchmark_betamatrix <- mixture_sim_mat

probes_select <- intersect(probes, rownames(benchmark_betamatrix))  ## 600
library(EpiDISH)
source("projectCellType.R")
Houseman_res <- projectCellType(benchmark_betamatrix[probes_select,],as.matrix(compTable[probes_select,3:9]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:9]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:9]), method = "CBS")$estF
# 
corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(proportions_sim[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], CBS_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], method = "spearman")
}


corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(proportions_sim[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], Houseman_res[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], method = "spearman")
}




corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(proportions_sim[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK","Epithelial")], Houseman_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "Epithelial")], method = "spearman")
}

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(proportions_sim[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK","Epithelial")], RPC_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "Epithelial")], method = "spearman")
}


corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(proportions_sim[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK","Epithelial")], CBS_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "Epithelial")], method = "spearman")
}


corr <- rep(NA, 7)
for (i in 1:7){
  corr[i] <- cor(proportions_sim[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "Epithelial")][,i], Houseman_res[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "Epithelial")][,i], method = "spearman")
}



