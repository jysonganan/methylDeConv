## nonlinear simulation

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


