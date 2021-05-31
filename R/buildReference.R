#'Build the reference library for 450k arrays
#'
#'Build the reference library for 450k to include "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran",
#'And extend the reference library for 450k to include "Epithelial", "Fibroblast".
#'@param extend If TRUE, generate the extended reference library; otherwise, generate the reference library of
#'only 6 immune cell types. Default value is TRUE.
#'@return A list of beta value reference matrix (ref_betamatrix) and reference cell types (ref_phenotype).
#'@export


build_reference_450k <- function(extend = TRUE){

  library(FlowSorted.Blood.450k)
  CellLines.matrix = NULL
  cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
  ## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
  ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
  ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
  keep <- which(ref_phenotype %in% cellTypes)
  ref_betamatrix <- ref_betamatrix[,keep]
  ref_phenotype <- ref_phenotype[keep]

  if (extend == FALSE){
    return(list(ref_betamatrix = ref_betamatrix, ref_phenotype = ref_phenotype))
  }
  else{
    library(GEOquery)
    library(minfi)
    rgSet <- read.metharray.exp("GSE40699/idat")
    ref_betamatrix_1 <- getBeta(preprocessNoob(rgSet, dyeMethod = "single"))
    geoMat <- getGEO("GSE40699")
    pD.all <- pData(geoMat[[1]])

    ref_betamatrix_EpiFib <- cbind(ref_betamatrix_1[,match(rownames(pD.all)[c(2,12,27,21,28,35,44,46,50,51,56)],substr(colnames(rgSet),1,9))],
                                   ref_betamatrix_1[,match(rownames(pD.all)[c(6,8,10,11,14,16,60)],substr(colnames(rgSet),1,9))])
    ref_phenotype_EpiFib <- c(rep("Epithelial", 11), rep("Fibroblast", 7))

    ref_betamatrix <- cbind(ref_betamatrix, ref_betamatrix_EpiFib)
    ref_phenotype <- c(ref_phenotype, ref_phenotype_EpiFib)

    return(list(ref_betamatrix = ref_betamatrix, ref_phenotype = ref_phenotype))
  }
 }







#'Build the reference library for EPIC arrays
#'
#'Build the reference library for EPIC arrays to include "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu",
#'And extend the reference library for 450k to include "Epithelial".
#'@param extend If TRUE, generate the extended reference library; otherwise, generate the reference library of
#'only 6 immune cell types. Default value is TRUE.
#'@return A list of beta value reference matrix (ref_betamatrix) and reference cell types (ref_phenotype).
#'@export


build_reference_EPIC <- function(extend = TRUE){
  library(ExperimentHub)
  hub <- ExperimentHub()
  query(hub, "FlowSorted.Blood.EPIC")
  FlowSorted.Blood.EPIC <- hub[["EH1136"]]

  library(GEOquery)
  library(minfi)
  CellLines.matrix = NULL
  cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
  ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.EPIC, dyeMethod = "single"))
  ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.EPIC))$CellType
  keep <- which(ref_phenotype %in% cellTypes)
  ref_betamatrix <- ref_betamatrix[,keep]
  ref_phenotype <- ref_phenotype[keep]

  if (extend == FALSE){
    return(list(ref_betamatrix = ref_betamatrix, ref_phenotype = ref_phenotype))
  }
  else{
    getGEOSuppFiles("GSE122126")
    untar("GSE122126/GSE122126_RAW.tar", exdir = "GSE122126/idat")

    idatFiles <- list.files("GSE122126/idat", pattern = "idat.gz$", full = TRUE)
    sapply(idatFiles, gunzip, overwrite = TRUE)
    rgSet <- read.metharray.exp("GSE122126/idat", force = TRUE)
    grSet <- preprocessNoob(rgSet, dyeMethod = "single")
    betaMat <- getBeta(grSet) #866091     90

    geoMat <- getGEO("GSE122126")
    pD.all <- pData(geoMat[[2]])

    colnames(betaMat) <- substr(colnames(betaMat),1,10)
    betaMat_122126 <- betaMat[,rownames(pD.all)]
    phenotype_122126 <- pD.all[,"sample type:ch1"]

    betaMat_122126_sub <- betaMat_122126[,phenotype_122126%in%c("Colon epithelial cells","Lung epithelial cells",
                                                                "Pancreatic acinar cells","Pancreatic duct cells")]
    phenotype_122126_sub <- rep("Epithelial",10)

    ref_betamatrix <- cbind(ref_betamatrix, betaMat_122126_sub)
    ref_phenotype <- c(ref_phenotype, phenotype_122126_sub)
    return(list(ref_betamatrix = ref_betamatrix, ref_phenotype = ref_phenotype))
  }

}











#'Build the reference library for 450k arrays of brain tissue
#'
#'Build the reference library for 450k arrays of brain tissues (DLPFC) to include "NeuN_neg", "NeuN_pos".
#'@return A list of beta value reference matrix (ref_betamatrix) and reference cell types (ref_phenotype).
#'@export

build_reference_450k_neuron <- function(){
  library(FlowSorted.DLPFC.450k)
  library(minfi)
  GRset_frontCortex <-  preprocessNoob(FlowSorted.DLPFC.450k, dyeMethod = "single")
  CellLines.matrix = NULL
  cellTypes = c("NeuN_neg","NeuN_pos")
  ref_betamatrix <- getBeta(GRset_frontCortex)
  ref_phenotype <- as.data.frame(colData(FlowSorted.DLPFC.450k))$CellType
  keep <- which(ref_phenotype %in% cellTypes)
  ref_betamatrix <- ref_betamatrix[,keep]
  ref_phenotype <- ref_phenotype[keep]
  return(list(ref_betamatrix = ref_betamatrix, ref_phenotype = ref_phenotype))
}


#'Build the reference library for EPIC arrays of brain tissue
#'
#'Build the reference library for EPIC arrays of brain tissues to include "NeuN_neg", "NeuN_pos".
#'@return A list of beta value reference matrix (ref_betamatrix) and reference cell types (ref_phenotype).
#'@export

build_reference_EPIC_neuron <- function(){
  library(GEOquery)
  library(minfi)
  geoMat <- getGEO("GSE111165")

  pD.all <- pData(geoMat[[2]])
  pD_ref <- pD.all[pD.all[,"tissue:ch1"]%in%c("brain_neg","brain_pos"),]

  rgSet<- read.metharray.exp("GSE111165/idat",force = TRUE)
  grSet <- preprocessNoob(rgSet, dyeMethod = "single")
  betaMat <- getBeta(grSet)
  colnames(betaMat) <- substr(colnames(betaMat),1,10)
  ref_betamatrix <- betaMat[,rownames(pD_ref)]
  ref_phenotype <- c(rep("NeuN_neg",12),rep("NeuN_pos",5))
  return(list(ref_betamatrix = ref_betamatrix, ref_phenotype = ref_phenotype))
}





