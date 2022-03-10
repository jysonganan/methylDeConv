### stability analysis of different ref sample size


#################################### 
###### EPIC Epithelial (no cfDNA)
####################################
#save("betaMat_122126","phenotype_122126", file = "ref_122126_EPICEpithelial.RData")
data_type = "FlowEPIC_Epithelial_nocfDNA"
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
betaMat_122126_sub <- betaMat_122126[,phenotype_122126%in%c("Colon epithelial cells","Lung epithelial cells",
                                                            "Pancreatic acinar cells","Pancreatic duct cells")]
phenotype_122126_sub <- rep("Epithelial",10)

ref_betamatrix <- cbind(ref_betamatrix, betaMat_122126_sub)
ref_phenotype <- c(ref_phenotype, phenotype_122126_sub)





stabilityAnalysis <- function(ref_phenotype, ref_betamatrix, ref_sample_size = 5, benchmark_betamatrix = NULL, benchmark_trueprop = NULL, returnCompTable = FALSE){
  newID <- c(sample(which(ref_phenotype == "Bcell"),ref_sample_size,replace = FALSE),
             sample(which(ref_phenotype == "CD4T"),ref_sample_size,replace = FALSE),
             sample(which(ref_phenotype == "CD8T"),ref_sample_size,replace = FALSE),
             sample(which(ref_phenotype == "Epithelial"),ref_sample_size,replace = FALSE),
             sample(which(ref_phenotype == "Mono"),ref_sample_size,replace = FALSE),
             sample(which(ref_phenotype == "Neu"),ref_sample_size,replace = FALSE),
             sample(which(ref_phenotype == "NK"),ref_sample_size,replace = FALSE))
  
  ref_betamatrix_sub <- ref_betamatrix[,newID]
  ref_phenotype_sub <- ref_phenotype[newID]
  source("refCompTableProbeSelection.R")
  compTable <- ref_compTable(ref_betamatrix_sub, ref_phenotype_sub)
  
  nCores = 4
  
  probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix_sub, ref_phenotype_sub, probeSelect = "both")
  
  probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix_sub, ref_phenotype_sub, probeSelect = "both", MaxDMRs = 300)
  ProbePreselect_multiclassGlmnet <- ref_probe_selection_multiclassGlmnet_cv(ref_betamatrix_sub[probes,], ref_phenotype_sub)
  probes_select <- ProbePreselect_multiclassGlmnet[[1]][-1]
  
  if (!(is.null(benchmark_betamatrix))){
    library(EpiDISH)
    source("projectCellType.R")
    ## onevsall
    Houseman_res <- projectCellType(benchmark_betamatrix[probes_oneVsAllttest,],as.matrix(compTable[probes_oneVsAllttest,3:9]))
    RPC_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_oneVsAllttest,3:9]), method = "RPC")$estF
    CBS_res <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_oneVsAllttest,3:9]), method = "CBS")$estF
    #1
    corr <- rep(NA, ncol(benchmark_betamatrix))
    for (i in 1:ncol(benchmark_betamatrix)){
      corr[i] <-cor(as.numeric(Houseman_res[i,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]),
                    as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
    }
    print(mean(corr))
    #2
    corr <- rep(NA, ncol(benchmark_betamatrix))
    for (i in 1:ncol(benchmark_betamatrix)){
      corr[i] <-cor(as.numeric(RPC_res[i,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]),
                    as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
    }
    print(mean(corr))
    #3
    corr <- rep(NA, ncol(benchmark_betamatrix))
    for (i in 1:ncol(benchmark_betamatrix)){
      corr[i] <-cor(as.numeric(CBS_res[i,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]),
                    as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
    }
    print(mean(corr))
    
    #4
    corr <- rep(NA, ncol(benchmark_trueprop))
    for (i in 1:ncol(benchmark_trueprop)){
      corr[i] <-cor(as.numeric(Houseman_res[,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")][,i]),
                    as.numeric(as.character(benchmark_trueprop[,i])),method = "spearman")
    }
    print(corr)
    
    ## glmnet preselect
    Houseman_res_glmnet <- projectCellType(benchmark_betamatrix[probes_select,],as.matrix(compTable[probes_select,3:9]))
    RPC_res_glmnet <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:9]), method = "RPC")$estF
    CBS_res_glmnet <- epidish(benchmark_betamatrix, as.matrix(compTable[probes_select,3:9]), method = "CBS")$estF
    
    #1
    corr <- rep(NA, ncol(benchmark_betamatrix))
    for (i in 1:ncol(benchmark_betamatrix)){
      corr[i] <-cor(as.numeric(Houseman_res_glmnet[i,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]),
                    as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
    }
    print(mean(corr))
    #2
    corr <- rep(NA, ncol(benchmark_betamatrix))
    for (i in 1:ncol(benchmark_betamatrix)){
      corr[i] <-cor(as.numeric(RPC_res_glmnet[i,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]),
                    as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
    }
    print(mean(corr))
    #3
    corr <- rep(NA, ncol(benchmark_betamatrix))
    for (i in 1:ncol(benchmark_betamatrix)){
      corr[i] <-cor(as.numeric(CBS_res_glmnet[i,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]),
                    as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
    }
    print(mean(corr))
    
    #4
    corr <- rep(NA, ncol(benchmark_trueprop))
    for (i in 1:ncol(benchmark_trueprop)){
      corr[i] <-cor(as.numeric(Houseman_res_glmnet[,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")][,i]),
                    as.numeric(as.character(benchmark_trueprop[,i])),method = "spearman")
    }
    print(corr)
    
    ## predictive modeling of preselect + glmnet
    probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix_sub, ref_phenotype_sub,probeSelect = "both", MaxDMRs = 300)
    library(dplyr)
    multiGlmnet_predProb <- predict(ProbePreselect_multiclassGlmnet[[2]], newdata = t(benchmark_betamatrix[probes,]), type = "prob") %>% 
      mutate('class'=names(.)[apply(., 1, which.max)])
    rownames(multiGlmnet_predProb) <- colnames(benchmark_betamatrix)
    
    corr <- rep(NA, ncol(benchmark_betamatrix))
    for (i in 1:ncol(benchmark_betamatrix)){
      corr[i] <-cor(as.numeric(multiGlmnet_predProb[i,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]),as.numeric(as.character(benchmark_trueprop[i,])),method = "spearman")
    }
    print(mean(corr))
    corr <- rep(NA, ncol(benchmark_trueprop))
    for (i in 1:ncol(benchmark_trueprop)){
      corr[i] <-cor(as.numeric(multiGlmnet_predProb[,c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")][,i]),as.numeric(as.character(benchmark_trueprop[,i])),method = "spearman")
    }
    print(corr)
  
    
  }
    
  if (returnCompTable == TRUE){
    return(list(probes_oneVsAllttest, probes_select, compTable))
  }
  else{
    return(list(probes_oneVsAllttest, probes_select))
  }
  
}




set.seed(2)
probes = list()
for (i in 1:20){
  probes[[i]] = list()
  res <- stabilityAnalysis(ref_phenotype, ref_betamatrix, ref_sample_size = 5, benchmark_betamatrix, benchmark_trueprop)
  probes[[i]][[1]] <- res[[1]]
  probes[[i]][[2]] <- res[[2]]

}
### check consensus probes selected
a = probes[[1]][[1]]
for (i in 1:19){
  a <- intersect(a, probes[[i+1]][[1]])
}

b = probes[[1]][[2]]
for (i in 1:19){
  b <- intersect(b, probes[[i+1]][[2]])
}









### stability on EWAS data
set.seed(2)
library(EpiDISH)
source("projectCellType.R")
library(ggplot2)
library(tidyr)
#### input betaMat of EWAS

p <- list()
for (i in 1:10){
  p[[i]] = list()
  probes <- stabilityAnalysis(ref_phenotype, ref_betamatrix, ref_sample_size = 5, benchmark_betamatrix = NULL, benchmark_trueprop = NULL, returnCompTable = TRUE)
  compTable <- probes[[3]]
  probes_select <- probes[[1]]
  Houseman_res_epicEpithelial <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:9]))
  probes_select <- probes[[2]]
  Houseman_res_glmnetpreselect_epicEpithelial <- projectCellType(betaMat[probes_select,],as.matrix(compTable[probes_select,3:9]))
  
  
  group <- pD[,"condition:ch1"]
  gender_dat_blood <- matrix(NA, 30, 8)
  samples <- rownames(Houseman_res_epicEpithelial)
  rownames(gender_dat_blood) <- samples
  gender_dat_blood[,1:7] <- Houseman_res_epicEpithelial
  gender_dat_blood[,8] <- group
  gender_dat_blood <- as.data.frame(gender_dat_blood)
  gender_dat_blood[1:7] <- apply(gender_dat_blood[1:7], 2, as.numeric)
  gender_dat_blood[,8] <- as.character(gender_dat_blood[,8])
  colnames(gender_dat_blood) <- c(colnames(Houseman_res_epicEpithelial),"group")


  df <- gather(gender_dat_blood, series,value,-group)
  p[[i]][[1]]<- ggplot(df) + geom_boxplot(aes(series ,value,color=group)) +
    xlab('cell types')+
    ylab('proportions') +
    ggtitle("GSE112308-Melanoma_Houseman-onevsAllttest (EPICEpithelial)")

  
  
  group <- pD[,"condition:ch1"]
  gender_dat_blood <- matrix(NA,30,8)
  samples <- rownames(Houseman_res_glmnetpreselect_epicEpithelial)
  rownames(gender_dat_blood) <- samples
  gender_dat_blood[,1:7] <- Houseman_res_glmnetpreselect_epicEpithelial
  gender_dat_blood[,8] <- group
  gender_dat_blood <- as.data.frame(gender_dat_blood)
  gender_dat_blood[1:7] <- apply(gender_dat_blood[1:7], 2, as.numeric)
  gender_dat_blood[,8] <- as.character(gender_dat_blood[,8])
  colnames(gender_dat_blood) <- c(colnames(Houseman_res_glmnetpreselect_epicEpithelial),"group")
  

  df <- gather(gender_dat_blood, series,value,-group)
  p[[i]][[2]]<- ggplot(df) + geom_boxplot(aes(series ,value,color=group)) +
    xlab('cell types')+
    ylab('proportions') +
    ggtitle("GSE112308-Melanoma_Houseman-glmnetPreselect (EPICEpithelial)")
}


library(gridExtra)

pdf("plots_qsub.pdf", onefile = TRUE, width = 10, height = 8)
for (i in seq(length(p))) {
  do.call("grid.arrange", p[[i]])  
}
dev.off()
