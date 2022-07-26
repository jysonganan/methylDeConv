#### feature selection comparison on Benchmark dataset 1
library(methylDeConv)
####################################
#### Benchmark dataset 1
####################################
library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "FlowSorted.Blood.EPIC")
FlowSorted.Blood.EPIC <- hub[["EH1136"]]
annot <- as.data.frame(colData(FlowSorted.Blood.EPIC))
benchmark <- which(annot$CellType == "MIX")
tmp <- getBeta(preprocessNoob(FlowSorted.Blood.EPIC, dyeMethod = "single"))
benchmark_betamatrix <- tmp[,rownames(annot)[benchmark]]
benchmark_trueprop <- annot[benchmark, c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]

####################################
#### build reference
####################################
reference_EPIC <- build_reference_EPIC(extend = FALSE)
compTable_EPIC <- ref_compTable(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype)
compTable_EPIC <- compTable_EPIC[,3:8]

####################################
#### feature selection
####################################
set.seed(2)
probe_1 <- ref_probe_selection_oneVsAllttest(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype)
probe_2 <- ref_probe_selection_oneVsAllLimma(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype)
probe_3 <- ref_probe_selection_pairwiseLimma(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype)
probe_4 <- ref_probe_selection_pairwiseGlmnet_cv(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype)
model_multiEN <- ref_probe_selection_multiclassGlmnet_cv(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype)
probe_5 <- model_multiEN[[1]]

set.seed(2)
probe_6 <- ref_probe_selection_twoStage(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype,
                                        preselect = 300, ml_model = "elastic net")


set.seed(2)
tmp <- ref_probe_selection_twoStage(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype,
                                    preselect = 100, ml_model = "RF")
tmp <- ref_probe_selection_twoStage(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype,
                                    preselect = 200, ml_model = "RF")
probe_7 <- ref_probe_selection_twoStage(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype,
                                        preselect = 300, ml_model = "RF")

set.seed(2)
probe_8 <- ref_probe_selection_twoStage(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype,
                                        preselect = 300, ml_model = "rfe")

set.seed(2)
model_RF<- ref_probe_selection_multiclassRF_cv(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype,
                                               tune_grid = c(300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000))
importance <- varImp(model_RF, scale=FALSE)
probe_9 <- rownames(importance$importance)[1:600]



####################################
#### deconvolution
####################################
Houseman_res_1 <- Houseman_project(benchmark_betamatrix, compTable_EPIC, probe_1)
RPC_res_1 <- RPC(benchmark_betamatrix, compTable_EPIC, probe_1)
CBS_res_1 <- CBS(benchmark_betamatrix, compTable_EPIC, probe_1)
MethylResolver_res_1 <- MethylResolver(benchmark_betamatrix, compTable_EPIC, probe_1)

Houseman_res_2 <- Houseman_project(benchmark_betamatrix, compTable_EPIC, probe_2)
RPC_res_2 <- RPC(benchmark_betamatrix, compTable_EPIC, probe_2)
CBS_res_2 <- CBS(benchmark_betamatrix, compTable_EPIC, probe_2)
MethylResolver_res_2 <- MethylResolver(benchmark_betamatrix, compTable_EPIC, probe_2)

Houseman_res_3 <- Houseman_project(benchmark_betamatrix, compTable_EPIC, probe_3)
RPC_res_3 <- RPC(benchmark_betamatrix, compTable_EPIC, probe_3)
CBS_res_3 <- CBS(benchmark_betamatrix, compTable_EPIC, probe_3)
MethylResolver_res_3 <- MethylResolver(benchmark_betamatrix, compTable_EPIC, probe_3)

Houseman_res_4 <- Houseman_project(benchmark_betamatrix, compTable_EPIC, probe_4)
RPC_res_4 <- RPC(benchmark_betamatrix, compTable_EPIC, probe_4)
CBS_res_4 <- CBS(benchmark_betamatrix, compTable_EPIC, probe_4)
MethylResolver_res_4 <- MethylResolver(benchmark_betamatrix, compTable_EPIC, probe_4)

Houseman_res_5 <- Houseman_project(benchmark_betamatrix, compTable_EPIC, probe_5)
RPC_res_5 <- RPC(benchmark_betamatrix, compTable_EPIC, probe_5)
CBS_res_5 <- CBS(benchmark_betamatrix, compTable_EPIC, probe_5)
MethylResolver_res_5 <- MethylResolver(benchmark_betamatrix, compTable_EPIC, probe_5)

Houseman_res_6 <- Houseman_project(benchmark_betamatrix, compTable_EPIC, probe_6)
RPC_res_6 <- RPC(benchmark_betamatrix, compTable_EPIC, probe_6)
CBS_res_6 <- CBS(benchmark_betamatrix, compTable_EPIC, probe_6)
MethylResolver_res_6 <- MethylResolver(benchmark_betamatrix, compTable_EPIC, probe_6)

Houseman_res_7 <- Houseman_project(benchmark_betamatrix, compTable_EPIC, probe_7)
RPC_res_7 <- RPC(benchmark_betamatrix, compTable_EPIC, probe_7)
CBS_res_7 <- CBS(benchmark_betamatrix, compTable_EPIC, probe_7)
MethylResolver_res_7 <- MethylResolver(benchmark_betamatrix, compTable_EPIC, probe_7)

Houseman_res_8 <- Houseman_project(benchmark_betamatrix, compTable_EPIC, probe_8)
RPC_res_8 <- RPC(benchmark_betamatrix, compTable_EPIC, probe_8)
CBS_res_8 <- CBS(benchmark_betamatrix, compTable_EPIC, probe_8)
MethylResolver_res_8 <- MethylResolver(benchmark_betamatrix, compTable_EPIC, probe_8)

Houseman_res_9 <- Houseman_project(benchmark_betamatrix, compTable_EPIC, probe_9)
RPC_res_9 <- RPC(benchmark_betamatrix, compTable_EPIC, probe_9)
CBS_res_9 <- CBS(benchmark_betamatrix, compTable_EPIC, probe_9)
MethylResolver_res_9 <- MethylResolver(benchmark_betamatrix, compTable_EPIC, probe_9)
####################################
#### evaluation
####################################

## within sample correalations with true proportions; average over samples
within_sample_corr <- function(true_proportions, deconv_res){
  corr <- rep(NA, 12)
  for (i in 1:12){
    corr[i] <- cor(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]),
                   as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), method = "spearman")
  }
  print(mean(corr))
  return(mean(corr))
}


## Add RMSE, MAPE metrics
library(Metrics)
within_sample_RMSE <- function(true_proportions, deconv_res){
  RMSE <- rep(NA, 12)
  for (i in 1:12){
    RMSE[i] <- rmse(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), 
                    as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]))
  }
  print(mean(RMSE))
  return(mean(RMSE))
}



within_sample_MAPE <- function(true_proportions, deconv_res){
  MAPE <- rep(NA, 12)
  for (i in 1:12){
    MAPE[i] <- mape(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), 
                    as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]))
  }
  print(mean(MAPE))
  return(mean(MAPE))
}




print("oneVsAllttest")
cor_Houseman_1 <- within_sample_corr(benchmark_trueprop, Houseman_res_1)
cor_RPC_1 <-within_sample_corr(benchmark_trueprop, RPC_res_1)
cor_CBS_1 <- within_sample_corr(benchmark_trueprop, CBS_res_1)
cor_MethylResolver_1 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_1)

print("oneVsAllLimma")
cor_Houseman_2 <- within_sample_corr(benchmark_trueprop, Houseman_res_2)
cor_RPC_2 <- within_sample_corr(benchmark_trueprop, RPC_res_2)
cor_CBS_2 <- within_sample_corr(benchmark_trueprop, CBS_res_2)
cor_MethylResolver_2 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_2)


print("pairwiseLimma")
cor_Houseman_3 <- within_sample_corr(benchmark_trueprop, Houseman_res_3)
cor_RPC_3 <- within_sample_corr(benchmark_trueprop, RPC_res_3)
cor_CBS_3 <- within_sample_corr(benchmark_trueprop, CBS_res_3)
cor_MethylResolver_3 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_3)

print("pairwiseGlmnet")
cor_Houseman_4 <- within_sample_corr(benchmark_trueprop, Houseman_res_4)
cor_RPC_4 <- within_sample_corr(benchmark_trueprop, RPC_res_4)
cor_CBS_4 <- within_sample_corr(benchmark_trueprop, CBS_res_4)
cor_MethylResolver_4 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_4)

print("multiGlmnet on all probes")
cor_Houseman_5 <- within_sample_corr(benchmark_trueprop, Houseman_res_5)
cor_RPC_5 <- within_sample_corr(benchmark_trueprop, RPC_res_5)
cor_CBS_5 <- within_sample_corr(benchmark_trueprop, CBS_res_5)
cor_MethylResolver_5 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_5)

print("glmnetpreselect")
cor_Houseman_6 <- within_sample_corr(benchmark_trueprop, Houseman_res_6)
cor_RPC_6 <- within_sample_corr(benchmark_trueprop, RPC_res_6)
cor_CBS_6 <- within_sample_corr(benchmark_trueprop, CBS_res_6)
cor_MethylResolver_6 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_6)

print("RFpreselect")
cor_Houseman_7 <- within_sample_corr(benchmark_trueprop, Houseman_res_7)
cor_RPC_7 <- within_sample_corr(benchmark_trueprop, RPC_res_7)
cor_CBS_7 <- within_sample_corr(benchmark_trueprop, CBS_res_7)
cor_MethylResolver_7 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_7)

print("RFEpreselect")
cor_Houseman_8 <- within_sample_corr(benchmark_trueprop, Houseman_res_8)
cor_RPC_8 <- within_sample_corr(benchmark_trueprop, RPC_res_8)
cor_CBS_8 <- within_sample_corr(benchmark_trueprop, CBS_res_8)
cor_MethylResolver_8 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_8)

print("multiRF on all probes")
cor_Houseman_9 <- within_sample_corr(benchmark_trueprop, Houseman_res_9)
cor_RPC_9 <- within_sample_corr(benchmark_trueprop, RPC_res_9)
cor_CBS_9 <- within_sample_corr(benchmark_trueprop, CBS_res_9)
cor_MethylResolver_9 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_9)


