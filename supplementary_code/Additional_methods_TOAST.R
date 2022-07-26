##############################
### Reference-free algorithms
##############################
load("./Downloads/Benchmark1_featureSelection.rda")

library(TOAST)

############# 1. most variable probes,  RefFreeEWAS
refinx <- findRefinx(benchmark_betamatrix, nmarker = 1000)
Y <- benchmark_betamatrix[refinx,]
K <- 6
outT1 <- myRefFreeCellMix(Y, mu0=myRefFreeCellMixInitialize(Y, K = K))
estProp_RF_mostVariable <- outT1$Omega

############# 2. use oneVsAllttest probes (probe_1) from reference dataset, RefFreeEWAS
Y <- benchmark_betamatrix[probe_1,]
K <- 6
outT2 <- myRefFreeCellMix(Y, mu0=myRefFreeCellMixInitialize(Y, K = K))
estProp_RF_semi <- outT2$Omega

############# 3. most variable probes, improved reference-free deconvolution with cross-cell type differential analysis
refinx <- findRefinx(benchmark_betamatrix, nmarker = 1000)
InitNames <- rownames(benchmark_betamatrix)[refinx]
outRF1 <- csDeconv(benchmark_betamatrix, K = K, nMarker = 1000, 
                   InitMarker = InitNames, TotalIter = 30,
                   bound_negative = TRUE)
estProp_RF_mostVariable_Improved <- outRF1$estProp

############# 4. use oneVsAllttest probes from reference dataset, improved reference-free deconvolution with cross-cell type differential analysis

outRF2 <- csDeconv(benchmark_betamatrix, K = K, nMarker = 1000, 
                   InitMarker = probe_1, TotalIter = 30,
                   bound_negative = TRUE)
estProp_RF_semi_Improved <- outRF2$estProp




############# 5. Tsisal 
out1 = Tsisal(benchmark_betamatrix, K = K, knowRef = NULL)
estProp_Tsisal <- out1$estProp

############# 6. semi-supervised Tsisal
out2 = Tsisal(benchmark_betamatrix, K = K, knowRef = reference_EPIC$ref_betamatrix)
estProp_Tsisal_semi <- out2$estProp




##########################################
### post-processing TOAST deconv results
##########################################
postprocess <- function(estProp){
  estProp[estProp < 0] <- 0
  res <- apply(estProp, 1, function(x){return(x/sum(x))})
  res <- t(res)
  return(res)
}


estProp_RF_mostVariable <- postprocess(estProp_RF_mostVariable)
estProp_RF_semi <- postprocess(estProp_RF_semi)
estProp_RF_mostVariable_Improved <- postprocess(estProp_RF_mostVariable_Improved)
estProp_RF_semi_Improved <- postprocess(estProp_RF_semi_Improved)
estProp_Tsisal <- postprocess(estProp_Tsisal)
estProp_Tsisal_semi <- postprocess(estProp_Tsisal_semi)

####################################################################################
# we align the cell types from RF and RB estimations using pearson's correlation
####################################################################################
estProp_RF_mostVariable <- assignCellType(input = estProp_RF_mostVariable, reference = RPC_res_1) 
colnames(estProp_RF_mostVariable) <- colnames(RPC_res_1)
estProp_RF_semi <- assignCellType(input = estProp_RF_semi, reference = RPC_res_1) 
colnames(estProp_RF_semi) <- colnames(RPC_res_1)
estProp_RF_mostVariable_Improved <- assignCellType(input = estProp_RF_mostVariable_Improved, reference = RPC_res_1) 
colnames(estProp_RF_mostVariable_Improved) <- colnames(RPC_res_1)
estProp_RF_semi_Improved <- assignCellType(input = estProp_RF_semi_Improved, reference = RPC_res_1) 
colnames(estProp_RF_semi_Improved) <- colnames(RPC_res_1)
estProp_Tsisal <- assignCellType(input = estProp_Tsisal, reference = RPC_res_1) 
colnames(estProp_Tsisal) <- colnames(RPC_res_1)
estProp_Tsisal_semi <- assignCellType(input = estProp_Tsisal_semi, reference = RPC_res_1) 
colnames(estProp_Tsisal_semi) <- colnames(RPC_res_1)







####################################################################################
###### performance on benchmark dataset
####################################################################################
within_sample_corr(benchmark_trueprop, estProp_RF_mostVariable)   #0.1197837
within_sample_corr(benchmark_trueprop, estProp_RF_semi)     #0.6523315
within_sample_corr(benchmark_trueprop, estProp_RF_mostVariable_Improved)    #0.2326111
within_sample_corr(benchmark_trueprop, estProp_RF_semi_Improved)     #0.7123025
within_sample_corr(benchmark_trueprop, estProp_Tsisal)       #0.3080416
within_sample_corr(benchmark_trueprop, estProp_Tsisal_semi)     #0.3272978


benchmark_trueprop_percent <- t(apply(benchmark_trueprop, 1, function(x){return(x/sum(x))}))
within_sample_RMSE(benchmark_trueprop_percent, estProp_RF_mostVariable) #0.2940593
within_sample_RMSE(benchmark_trueprop_percent, estProp_RF_semi) #0.1732393
within_sample_RMSE(benchmark_trueprop_percent, estProp_RF_mostVariable_Improved) #0.3120699
within_sample_RMSE(benchmark_trueprop_percent, estProp_RF_semi_Improved) #0.1888424
within_sample_RMSE(benchmark_trueprop_percent, estProp_Tsisal) #0.1041043
within_sample_RMSE(benchmark_trueprop_percent, estProp_Tsisal_semi)  #0.1079503


