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

save("compTable_EPIC", "benchmark_betamatrix", "benchmark_trueprop", file = "Benchmark_EPICref.RData")





#################
## prepare data for ARIC algorithm
#################

#### get ARIC deconvolution results for simulated datasets with EPIC reference library
load("Benchmark1_featureSelection.rda")
probes <- list()
probes[[1]] <- probe_1
probes[[2]] <- probe_2
probes[[3]] <- probe_3
probes[[4]] <- probe_4[[-1]]
probes[[5]] <- probe_5[[-1]]
probes[[6]] <- probe_6[[-1]]
probes[[7]] <- probe_7
probes[[8]] <- probe_8
probes[[9]] <- probe_9

## 10 simulated datasets 5 for beta mixture, 5 for gaussian mixture
file_names <- c("beta_sim_data_nonimmune_level1", "beta_sim_data_nonimmune_level2",
                "beta_sim_data_nonimmune_level3", "beta_sim_data_nonimmune_level4",
                "beta_sim_data_nonimmune_level5",
                "gaussian_sim_data_nonimmune_level1", "gaussian_sim_data_nonimmune_level2",
                "gaussian_sim_data_nonimmune_level3", "gaussian_sim_data_nonimmune_level4",
                "gaussian_sim_data_nonimmune_level5")

for (i in 1:10){
  load(paste0(file_names[i],".RData"))
  for (i in 1:9){
    write.csv(sim_data[probes[[i]],],
              file = paste0(paste0("sim_data_probe_",i), ".csv"), row.names = TRUE)
  }
}


for (i in 1:9){
  write.csv(compTable_EPIC[probes[[i]],],
            file = paste0(paste0("compTable_EPIC_probe_",i), ".csv"), row.names = TRUE)
}

