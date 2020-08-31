# BetaMixture purified dataset (450k + Epi + Fib)

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


sim1 <- betaSim(ref_betamatrix, ref_phenotype)

set.seed(2)
purified_datasets_sim <- list()
for (i in 1:20){
  purified_datasets_sim[[i]] <-  betaSim(ref_betamatrix, ref_phenotype)
}

save(purified_datasets_sim, file = "purified_datasets_sim_1.RData")



for (i in 1:20){
  purified_datasets_sim[[i]][which(is.na(purified_datasets_sim[[i]])==TRUE)] <- 0  
}


library(MCMCpack)

##  around 0.2-0.3 epithelial; around 0.3-0.4 Fibroblasts
## 2. 
set.seed(4)
proportions_sim <- rdirichlet(30, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(30, 0.2, 0.3)
proportions_sim_fibroblast <- runif(30, 0.3, 0.4)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial - proportions_sim_fibroblast)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial, proportions_sim_fibroblast)
proportions_sim <- proportions_sim[,c(1,2,3,7,8,4,5,6)]

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

### 1. Flow450k + Epi + Fib
source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)

load("Flow450k_EpiFib_Probes.RData")

probes <- probes_oneVsAllttest
benchmark_betamatrix <- mixture_sim_mat_2
true_proportions <- true_proportions_sim_2

probes_select <- intersect(probes, rownames(benchmark_betamatrix))  ## 600
library(EpiDISH)
source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/projectCellType.R")
Houseman_res <- projectCellType(benchmark_betamatrix[probes_select,],as.matrix(compTable[probes_select,3:10]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:10]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:10]), method = "CBS")$estF

### 2. EpiDISH
data("centEpiFibIC.m")
use_probes <- intersect(rownames(benchmark_betamatrix),rownames(centEpiFibIC.m))
res <- projectCellType(benchmark_betamatrix[use_probes,], 
                       centEpiFibIC.m[use_probes,])
data("centDHSbloodDMC.m")
Houseman_res <- hepidish(benchmark_betamatrix, centEpiFibIC.m, centDHSbloodDMC.m, h.CT.idx = 3, method = "CP")
RPC_res <- hepidish(benchmark_betamatrix, centEpiFibIC.m, centDHSbloodDMC.m, h.CT.idx = 3, method = "RPC")
CBS_res <- hepidish(benchmark_betamatrix, centEpiFibIC.m, centDHSbloodDMC.m, h.CT.idx = 3, method = "CBS")
Houseman_res <- Houseman_res[, c(3,5,6,1,2,8,7,4)]
RPC_res <- RPC_res[, c(3,5,6,1,2,8,7,4)]
CBS_res <- CBS_res[, c(3,5,6,1,2,8,7,4)]


### 3. Flow450k
source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)


load("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/Flow450kProbesdefault.RData")
probes <- probes_oneVsAllttest
benchmark_betamatrix <- mixture_sim_mat_2
true_proportions <- true_proportions_sim_2

probes_select <- intersect(probes, rownames(benchmark_betamatrix))  ## 600
library(EpiDISH)
source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/projectCellType.R")
Houseman_res <- projectCellType(benchmark_betamatrix[probes_select,],as.matrix(compTable[probes_select,3:8]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:8]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:8]), method = "CBS")$estF


### 4 Flow450k + Epi
library(FlowSorted.Blood.450k)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]

load("ref_GSE40699_EpiFib.RData")
ref_betamatrix <- cbind(ref_betamatrix, ref_betamatrix_EpiFib)
ref_phenotype <- c(ref_phenotype, ref_phenotype_EpiFib)
ref_betamatrix <- ref_betamatrix[,1:47]
ref_phenotype <- ref_phenotype[1:47]

set.seed(2)
source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/refCompTableProbeSelection.R")
probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both")

compTable <- ref_compTable(ref_betamatrix, ref_phenotype)

probes <- probes_oneVsAllttest
benchmark_betamatrix <- mixture_sim_mat_2
true_proportions <- true_proportions_sim_2

probes_select <- intersect(probes, rownames(benchmark_betamatrix))  ## 600
library(EpiDISH)
source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/projectCellType.R")
Houseman_res <- projectCellType(benchmark_betamatrix[probes_select,],as.matrix(compTable[probes_select,3:9]))
RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:9]), method = "RPC")$estF
CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:9]), method = "CBS")$estF


## immune estimation
corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c(1,2,3,6,7,8)], Houseman_res[i,c(1,2,3,5,6,7)], method = "spearman")
}
mean(corr)

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c(1,2,3,6,7,8)], RPC_res[i,c(1,2,3,5,6,7)], method = "spearman")
}
mean(corr)

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c(1,2,3,6,7,8)], CBS_res[i,c(1,2,3,5,6,7)], method = "spearman")
}
mean(corr)

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(true_proportions[,c(1,2,3,6,7,8)][,i], Houseman_res[,c(1,2,3,5,6,7)][,i], method = "spearman")
}
corr

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(true_proportions[,c(1,2,3,6,7,8)][,i], RPC_res[,c(1,2,3,5,6,7)][,i], method = "spearman")
}
corr

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <- cor(true_proportions[,c(1,2,3,6,7,8)][,i], CBS_res[,c(1,2,3,5,6,7)][,i], method = "spearman")
}
corr





## immune + epithelial estimation
corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c(1,2,3,6,7,8,4)], Houseman_res[i,c(1,2,3,5,6,7,4)], method = "spearman")
}
mean(corr)

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c(1,2,3,6,7,8,4)], RPC_res[i,c(1,2,3,5,6,7,4)], method = "spearman")
}
mean(corr)

corr <- rep(NA, 600)
for (i in 1:600){
  corr[i] <- cor(true_proportions[i,c(1,2,3,6,7,8,4)], CBS_res[i,c(1,2,3,5,6,7,4)], method = "spearman")
}
mean(corr)

corr <- rep(NA, 7)
for (i in 1:7){
  corr[i] <- cor(true_proportions[,c(1,2,3,6,7,8,4)][,i], Houseman_res[,c(1,2,3,5,6,7,4)][,i], method = "spearman")
}
corr

corr <- rep(NA, 7)
for (i in 1:7){
  corr[i] <- cor(true_proportions[,c(1,2,3,6,7,8,4)][,i], RPC_res[,c(1,2,3,5,6,7,4)][,i], method = "spearman")
}
corr

corr <- rep(NA, 7)
for (i in 1:7){
  corr[i] <- cor(true_proportions[,c(1,2,3,6,7,8,4)][,i], CBS_res[,c(1,2,3,5,6,7,4)][,i], method = "spearman")
}
corr






### Normal mixture dataset  #############   (Have already saved on the server)!
## step 1 
## logit transform of the purified datasets
# ref_betamatrix min: 0.003; max: 0.997
ref_betamatrix_logit <- log(ref_betamatrix/(1-ref_betamatrix))
rownames(ref_betamatrix_logit) <- rownames(ref_betamatrix)
colnames(ref_betamatrix_logit) <- colnames(ref_betamatrix)

## step 2 
## purified dataset simulation, normal distribution

estNormalParams_datRow <- function(datRow){
  mu = mean(datRow, na.rm = TRUE)
  var = var(datRow, na.rm = TRUE)
  sd = sqrt(var)
  return(list(mu = mu, sd = sd))
}

normalSim <- function(ref_betamatrix, ref_phenotype){
  phenotype <- levels(as.factor(ref_phenotype))
  dat_sim <- matrix(NA, nrow(ref_betamatrix), length(phenotype))
  for (i in 1:length(phenotype)){
    dat <- ref_betamatrix[,ref_phenotype == phenotype[i]]
    dat_sim[,i] <- apply(dat, 1, function(x){
      para <- estNormalParams_datRow(x)
      return(rnorm(1, para[[1]], para[[2]]))
    })
  }
  rownames(dat_sim) <- rownames(ref_betamatrix)
  colnames(dat_sim) <- phenotype
  return(dat_sim)
}



set.seed(2)
purified_datasets_sim <- list()
for (i in 1:20){
  purified_datasets_sim[[i]] <-  normalSim(ref_betamatrix_logit, ref_phenotype)
}

save(purified_datasets_sim, file = "purified_datasets_sim_3.RData")

