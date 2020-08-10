## purified: immune + cfDNA
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

# cfDNA         cfDNA In vitro mix 
# 58                          5 
# **Colon epithelial cells           Cortical neurons 
# 3                          2 
# Hepatocytes               In vitro mix 
# 2                          9 
# Leukocytes      **Lung epithelial cells 
# 1                          3 
# **Pancreatic acinar cells      Pancreatic beta cells 
# 2                          1 
# ** Pancreatic duct cells Vascular endothelial cells 
# 2                          2 
set.seed(3)
betaMat_122126_sub <- betaMat_122126[,sample(which(phenotype_122126 == "cfDNA"),10,replace = FALSE)]
phenotype_122126_sub <- rep("cfDNA", 10)

ref_betamatrix <- cbind(ref_betamatrix, betaMat_122126_sub)
ref_phenotype <- c(ref_phenotype, phenotype_122126_sub)





#https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

estBetaParams_datRow <- function(datRow){
  mu = mean(datRow, na.rm = TRUE)
  var = var(datRow, na.rm = TRUE)
  return(estBetaParams(mu, var))
}

betaSim <- function(ref_betamatrix, ref_phenotype){
  phenotype <- levels(as.factor(ref_phenotype))
  dat_sim <- matrix(NA, nrow(ref_betamatrix), length(phenotype))
  for (i in 1:length(phenotype)){
    dat <- ref_betamatrix[,ref_phenotype == phenotype[i]]
    dat_sim[,i] <- apply(dat, 1, function(x){
      para <- estBetaParams_datRow(x)
      return(rbeta(1, para[[1]], para[[2]]))
    })
  }
  rownames(dat_sim) <- rownames(ref_betamatrix)
  colnames(dat_sim) <- phenotype
  return(dat_sim)
}



set.seed(2)
purified_datasets_sim <- list()
for (i in 1:20){
  purified_datasets_sim[[i]] <-  betaSim(ref_betamatrix, ref_phenotype)
}

save(purified_datasets_sim, file = "purified_datasets_sim_2.RData")

#save(purified_datasets_sim, file = "/sonas-hs/krasnitz/hpc/data/pfproj/purified_datasets_sim_2_new.RData")




## probes selection
set.seed(2)
source("refCompTableProbeSelection.R")
probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both")

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

save("probes_oneVsAllttest", "ProbePreselect_multiclassGlmnet", file = "/sonas-hs/krasnitz/hpc/data/pfproj/EPICcfDNAProbe.RData")





#### data simulation by multiplying purified dataset (simulated) with the proportions
save(purified_datasets_sim, file = "/sonas-hs/krasnitz/hpc/data/pfproj/purified_datasets_sim_2_new.RData")


## impute NAs with zeros
for (i in 1:20){
  purified_datasets_sim[[i]][which(is.na(purified_datasets_sim[[i]])==TRUE)] <- 0  
}

## Within in each of three groups (no epithelial, low epithelial and high epithelial)
## generate 600 mixture samples:

## combined the 20 ramdomly generated purified datasets (beta distribution) 
## with 30 randomly generated probablities to generate 600 mixture samples.


library(MCMCpack)

## 1. no cfDNA
set.seed(3)
proportions_sim <- rdirichlet(30, c(1,1,1,1,1,1))
proportions_sim <- cbind(proportions_sim, 0)
proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]

true_proportions_sim <- proportions_sim
for (i in 1:19){
  true_proportions_sim <- rbind(true_proportions_sim, proportions_sim)
}

mixture_sim_mat <- NULL
for (i in 1:20){
  #res = purified_datasets_sim[[i]] %*% t(proportions_sim)
  
  res <- matrix(NA, nrow(purified_datasets_sim[[i]]), 30)
  for (m in 1:nrow(res)){
    res[m,] <- purified_datasets_sim[[i]][m,] %*% t(proportions_sim)
  }
  
  mixture_sim_mat <- cbind(mixture_sim_mat, res)
}
rownames(mixture_sim_mat) <- rownames(purified_datasets_sim[[1]])

true_proportions_sim_1 <- true_proportions_sim
colnames(true_proportions_sim_1) <- colnames(purified_datasets_sim[[1]]) 
mixture_sim_mat_1 <- mixture_sim_mat



## 2. 
set.seed(4)
proportions_sim <- rdirichlet(30, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(30, 0.1, 0.2)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]

true_proportions_sim <- proportions_sim
for (i in 1:19){
  true_proportions_sim <- rbind(true_proportions_sim, proportions_sim)
}

mixture_sim_mat <- NULL
for (i in 1:20){
  #res = purified_datasets_sim[[i]] %*% t(proportions_sim)
  
  res <- matrix(NA, nrow(purified_datasets_sim[[i]]), 30)
  for (m in 1:nrow(res)){
    res[m,] <- purified_datasets_sim[[i]][m,] %*% t(proportions_sim)
  }
  
  mixture_sim_mat <- cbind(mixture_sim_mat, res)
}
rownames(mixture_sim_mat) <- rownames(purified_datasets_sim[[1]])

true_proportions_sim_2 <- true_proportions_sim
colnames(true_proportions_sim_2) <- colnames(purified_datasets_sim[[1]]) 
mixture_sim_mat_2 <- mixture_sim_mat



## 3
set.seed(11)
proportions_sim <- rdirichlet(30, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(30, 0.2, 0.5)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]


true_proportions_sim <- proportions_sim
for (i in 1:19){
  true_proportions_sim <- rbind(true_proportions_sim, proportions_sim)
}

mixture_sim_mat <- NULL
for (i in 1:20){
  #res = purified_datasets_sim[[i]] %*% t(proportions_sim)
  
  res <- matrix(NA, nrow(purified_datasets_sim[[i]]), 30)
  for (m in 1:nrow(res)){
    res[m,] <- purified_datasets_sim[[i]][m,] %*% t(proportions_sim)
  }
  
  mixture_sim_mat <- cbind(mixture_sim_mat, res)
}
rownames(mixture_sim_mat) <- rownames(purified_datasets_sim[[1]])

true_proportions_sim_3 <- true_proportions_sim
colnames(true_proportions_sim_3) <- colnames(purified_datasets_sim[[1]]) 
mixture_sim_mat_3 <- mixture_sim_mat


## 4
set.seed(7)
proportions_sim <- rdirichlet(30, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(30, 0.5, 0.8)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]


true_proportions_sim <- proportions_sim
for (i in 1:19){
  true_proportions_sim <- rbind(true_proportions_sim, proportions_sim)
}

mixture_sim_mat <- NULL
for (i in 1:20){
  #res = purified_datasets_sim[[i]] %*% t(proportions_sim)
  
  res <- matrix(NA, nrow(purified_datasets_sim[[i]]), 30)
  for (m in 1:nrow(res)){
    res[m,] <- purified_datasets_sim[[i]][m,] %*% t(proportions_sim)
  }
  
  mixture_sim_mat <- cbind(mixture_sim_mat, res)
}
rownames(mixture_sim_mat) <- rownames(purified_datasets_sim[[1]])

true_proportions_sim_4 <- true_proportions_sim
colnames(true_proportions_sim_4) <- colnames(purified_datasets_sim[[1]]) 
mixture_sim_mat_4 <- mixture_sim_mat



## 5. high cfDNA

set.seed(5)
proportions_sim <- rdirichlet(30, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(30, 0.8, 0.9)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]


true_proportions_sim <- proportions_sim
for (i in 1:19){
  true_proportions_sim <- rbind(true_proportions_sim, proportions_sim)
}

mixture_sim_mat <- NULL
for (i in 1:20){
  #res = purified_datasets_sim[[i]] %*% t(proportions_sim)
  
  res <- matrix(NA, nrow(purified_datasets_sim[[i]]), 30)
  for (m in 1:nrow(res)){
    res[m,] <- purified_datasets_sim[[i]][m,] %*% t(proportions_sim)
  }
  
  mixture_sim_mat <- cbind(mixture_sim_mat, res)
}
rownames(mixture_sim_mat) <- rownames(purified_datasets_sim[[1]])

true_proportions_sim_5 <- true_proportions_sim
colnames(true_proportions_sim_5) <- colnames(purified_datasets_sim[[1]]) 
mixture_sim_mat_5 <- mixture_sim_mat







##### analysis
# 1. using EPIC reference
load("FlowEPICProbesdefault.RData")
load("FlowEPICProbePreselect_multiclassGlmnet.RData")
source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)

probes <- probes_oneVsAllttest
benchmark_betamatrix <- mixture_sim_mat_5
true_proportions <- true_proportions_sim_5

probes_select <- intersect(probes, rownames(benchmark_betamatrix))  ## 600
library(EpiDISH)
source("projectCellType.R")
Houseman_res <- projectCellType(benchmark_betamatrix[probes_select,],as.matrix(compTable[probes_select,3:8]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:8]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:8]), method = "CBS")$estF

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], Houseman_res[i,], method = "spearman")
}
mean(corr)

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], RPC_res[i,], method = "spearman")
}
mean(corr)

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], CBS_res[i,], method = "spearman")
}
mean(corr)

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(true_proportions[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], Houseman_res[,i], method = "spearman")
}
corr

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(true_proportions[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], RPC_res[,i], method = "spearman")
}
corr


corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(true_proportions[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], CBS_res[,i], method = "spearman")
}
corr

#2 use EPIC + cfDNA reference
load("/sonas-hs/krasnitz/hpc/data/pfproj/EPICcfDNAProbe.RData")
source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)

#probes <- ProbePreselect_multiclassGlmnet[[1]][-1]
probes <- probes_oneVsAllttest

benchmark_betamatrix <- mixture_sim_mat_4
true_proportions <- true_proportions_sim_4


probes_select <- intersect(probes, rownames(benchmark_betamatrix))  ## 600
library(EpiDISH)
source("projectCellType.R")
Houseman_res <- projectCellType(benchmark_betamatrix[probes_select,],as.matrix(compTable[probes_select,3:9]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:9]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:9]), method = "CBS")$estF

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], 
                 Houseman_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], method = "spearman")
}
mean(corr)

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], 
                 RPC_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], method = "spearman")
}
mean(corr)

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], 
                 CBS_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")], method = "spearman")
}
mean(corr)


corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(true_proportions[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], 
                 Houseman_res[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], method = "spearman")
}
corr

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(true_proportions[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], 
                 RPC_res[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], method = "spearman")
}
corr



corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(true_proportions[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], 
                 CBS_res[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i], method = "spearman")
}
corr





corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK","cfDNA")], 
                 Houseman_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "cfDNA")], method = "spearman")
}
mean(corr)

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK","cfDNA")], 
                 RPC_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "cfDNA")], method = "spearman")
}
mean(corr)

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK","cfDNA")], 
                 CBS_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "cfDNA")], method = "spearman")
}
mean(corr)



corr <- rep(NA, 7)
for (i in 1:7){
  corr[i] <- cor(true_proportions[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "cfDNA")][,i], 
                 Houseman_res[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "cfDNA")][,i], method = "spearman")
}

corr


corr <- rep(NA, 7)
for (i in 1:7){
  corr[i] <- cor(true_proportions[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "cfDNA")][,i], 
                 RPC_res[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "cfDNA")][,i], method = "spearman")
}

corr



corr <- rep(NA, 7)
for (i in 1:7){
  corr[i] <- cor(true_proportions[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "cfDNA")][,i], 
                 CBS_res[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK", "cfDNA")][,i], method = "spearman")
}

corr




