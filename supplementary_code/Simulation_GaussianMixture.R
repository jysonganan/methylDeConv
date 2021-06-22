library(methylDeConv)

## Gaussian mixture simulation
### changing nonimmune_level = 1, 2, 3, 4, 5 corresponds to the 5 intervals in simulation.
### changing include = "Epithelial", "cfDNA".
nonimmune_level = 1
include = "Epithelial"

## generate purfied profiles from Extended reference library
reference_Extended <- build_reference_EPIC(extend = TRUE, include = include)
purified_datasets_sim <- gaussianSim_purifiedProfiles(reference_Extended$ref_betamatrix,
                                                  reference_Extended$ref_phenotype)

sim_data <- GaussianSim_mixtureProfiles(nonimmune_level = nonimmune_level, purified_datasets_sim)


## Deconvolution with EPIC reference library
reference_EPIC <- build_reference_EPIC(extend = FALSE)
compTable_EPIC <- ref_compTable(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype)
compTable_EPIC <- compTable_EPIC[,3:8]
probes_EPIC <- ref_probe_selection_oneVsAllttest(reference_EPIC$ref_betamatrix,
                                                 reference_EPIC$ref_phenotype)
Houseman_res_EPIC <- Houseman_project(sim_data$mixture_sim_mat, compTable_EPIC, probes_EPIC)
RPC_res_EPIC <- RPC(sim_data$mixture_sim_mat, compTable_EPIC, probes_EPIC)
CBS_res_EPIC <- CBS(sim_data$mixture_sim_mat, compTable_EPIC, probes_EPIC)
MethylResolver_res_EPIC <- MethylResolver(sim_data$mixture_sim_mat, compTable_EPIC, probes_EPIC)


## Deconvolution with Extended reference library
compTable_Extended <- ref_compTable(reference_Extended$ref_betamatrix, reference_Extended$ref_phenotype)
compTable_Extended <- compTable_Extended[,3:9]
probes_Extended <- ref_probe_selection_oneVsAllttest(reference_Extended$ref_betamatrix,
                                                     reference_Extended$ref_phenotype)
Houseman_res_Extended <- Houseman_project(sim_data$mixture_sim_mat, compTable_Extended, probes_Extended)
RPC_res_Extended <- RPC(sim_data$mixture_sim_mat, compTable_Extended, probes_Extended)
CBS_res_Extended <- CBS(sim_data$mixture_sim_mat, compTable_Extended, probes_Extended)
MethylResolver_res_Extended <- MethylResolver(sim_data$mixture_sim_mat, compTable_Extended, probes_Extended)



## within sample correalations with true proportions; average over samples
within_sample_corr <- function(true_proportions, deconv_res){
  corr <- rep(NA, 600)
  for (i in 1:600){
    corr[i] <- cor(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]),
                   as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), method = "spearman")
  }
  print(mean(corr))
}



within_sample_corr(sim_data$true_proportions_sim, Houseman_res_EPIC)
within_sample_corr(sim_data$true_proportions_sim, RPC_res_EPIC)
within_sample_corr(sim_data$true_proportions_sim, CBS_res_EPIC)
within_sample_corr(sim_data$true_proportions_sim, MethylResolver_res_EPIC)
within_sample_corr(sim_data$true_proportions_sim, Houseman_res_Extended)
within_sample_corr(sim_data$true_proportions_sim, RPC_res_Extended)
within_sample_corr(sim_data$true_proportions_sim, CBS_res_Extended)
within_sample_corr(sim_data$true_proportions_sim, MethylResolver_res_Extended)




## within cell type correalations with true proportions
within_celltype_corr <- function(true_proportions, deconv_res){
  corr <- rep(NA, 6)
  for (i in 1:6){
    corr[i] <- cor(as.numeric(true_proportions[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i]),
                   as.numeric(deconv_res[,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")][,i]), method = "spearman")
  }
  print(corr)
}

within_celltype_corr(sim_data$true_proportions_sim, Houseman_res_EPIC)
within_celltype_corr(sim_data$true_proportions_sim, RPC_res_EPIC)
within_celltype_corr(sim_data$true_proportions_sim, CBS_res_EPIC)
within_celltype_corr(sim_data$true_proportions_sim, MethylResolver_res_EPIC)
within_celltype_corr(sim_data$true_proportions_sim, Houseman_res_Extended)
within_celltype_corr(sim_data$true_proportions_sim, RPC_res_Extended)
within_celltype_corr(sim_data$true_proportions_sim, CBS_res_Extended)
within_celltype_corr(sim_data$true_proportions_sim, MethylResolver_res_Extended)
