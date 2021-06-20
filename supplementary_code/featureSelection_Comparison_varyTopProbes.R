#### feature selection comparison on Benchmark dataset 1
#### OneVsAllttest varying number of probes (top 50, 100, 150, 200)
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
probe_1 <- ref_probe_selection_oneVsAllttest(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype, MaxDMRs = 50)
probe_2 <- ref_probe_selection_oneVsAllttest(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype, MaxDMRs = 100)
probe_3 <- ref_probe_selection_oneVsAllttest(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype, MaxDMRs = 150)
probe_4 <- ref_probe_selection_oneVsAllttest(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype, MaxDMRs = 200)


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


print("oneVsAllttest top 50")
cor_Houseman_1 <- within_sample_corr(benchmark_trueprop, Houseman_res_1)
cor_RPC_1 <-within_sample_corr(benchmark_trueprop, RPC_res_1)
cor_CBS_1 <- within_sample_corr(benchmark_trueprop, CBS_res_1)
cor_MethylResolver_1 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_1)

print("oneVsAllttest top 100")
cor_Houseman_2 <- within_sample_corr(benchmark_trueprop, Houseman_res_2)
cor_RPC_2 <- within_sample_corr(benchmark_trueprop, RPC_res_2)
cor_CBS_2 <- within_sample_corr(benchmark_trueprop, CBS_res_2)
cor_MethylResolver_2 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_2)


print("oneVsAllttest top 150")
cor_Houseman_3 <- within_sample_corr(benchmark_trueprop, Houseman_res_3)
cor_RPC_3 <- within_sample_corr(benchmark_trueprop, RPC_res_3)
cor_CBS_3 <- within_sample_corr(benchmark_trueprop, CBS_res_3)
cor_MethylResolver_3 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_3)

print("oneVsAllttest top 200")
cor_Houseman_4 <- within_sample_corr(benchmark_trueprop, Houseman_res_4)
cor_RPC_4 <- within_sample_corr(benchmark_trueprop, RPC_res_4)
cor_CBS_4 <- within_sample_corr(benchmark_trueprop, CBS_res_4)
cor_MethylResolver_4 <- within_sample_corr(benchmark_trueprop, MethylResolver_res_4)

