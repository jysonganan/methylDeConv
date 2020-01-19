.isRGOrStop <- function(object){
  if (!is(object, "RGChannelSet")){
    stop("object is of class '", class(object), "', but needs to be of ",
         "class 'RGChannelSet' or 'RGChannelSetExtended'")
  }
}

.is450k <- function(object){
  annotation(object)["array"] == "IlluminaHumanMethylation450k"
}

.isEPIC <- function(object){
  annotation(object)["array"] == "IlluminaHumanMethylationEPIC"
}

.is450kOrStop <- function(object){
  if (!.is450k(object)){
    stop("object is not IlluminaHumanMethylation450k array!")
  }
}

.isEPICOrStop <- function(object){
  if (!.isEPIC(object)){
    stop("object is not IlluminaHumanMethylationEPIC array!")
  }
}

MethylDeconv <- function(input_methyl, method = "Houseman", normalized = TRUE, tissue = "Blood", custom_probes = NULL){
  if(normalized){
    # beta value input
    MethylDeconv_normalized(input_methyl, method, tissue, custom_probes)
  }
  else{
    ## input is RgSet
    library(minfi)
    RGset <- input_methyl
    .isRGOrStop(RGset)
    .is450kOrStop(RGset)
    if (tissue == "Breast"|tissue == "genericEpithelial"){
      stop("For breast tissue/ generic epithelial tissue, 
           only support normalized beta matrix, doesn't support RGset!")
    }
    if (!tissue == "CordBlood"){
      GRset.normalized <- preprocessQuantile(RGset, fixOutliers = TRUE,
                                           removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                           quantileNormalize = TRUE, stratified = TRUE, 
                                           mergeManifest = FALSE, sex = NULL)
      
    }else{
      GRset.normalized <- preprocessNoob(RGset)
    }
                                        
    output <- list()
    output[[1]] <- estimateCellCounts(RGset, compositeCellType = tissue, meanPlot = FALSE, returnAll = FALSE)
    output[[2]] <- MethylDeconv_normalized(getBeta(GRset.normalized), method, tissue, custom_probes)
    return(output)
  }
}

## Generate reference profiles
MethylDeconv_normalized <- function(input_methyl, method = "Houseman", tissue = "Blood", custom_probes = NULL){
  if (method == "Houseman" && tissue == "Blood"){
    library(FlowSorted.Blood.450k)
    source("projectCellType.R")
    if (is.null(custom_probes)){
      use_probes <- intersect(rownames(input_methyl), rownames(FlowSorted.Blood.450k.JaffeModelPars))
      res <- projectCellType(input_methyl[use_probes,], 
                             as.matrix(FlowSorted.Blood.450k.JaffeModelPars[use_probes,]))
    }else{
      use_probes <- intersect(intersect(rownames(input_methyl), custom_probes), rownames(FlowSorted.Blood.450k.compTable))
      res <- projectCellType(input_methyl[use_probes,],as.matrix(FlowSorted.Blood.450k.compTable[use_probes,3:8]))
    }
    return(res)
  }
  
  if (method == "Houseman" && tissue == "CordBlood"){
    library(FlowSorted.CordBlood.450k)
    source("projectCellType.R")
    if (is.null(custom_probes)){
      use_probes <- interect(rownames(input_methyl), rownames(FlowSorted.CordBlood.450k.ModelPars))
      res <- projectCellType(input_methyl[use_probes,], 
                             as.matrix(FlowSorted.CordBlood.450k.ModelPars[use_probes,]))
    }else{
      use_probes <- intersect(intersect(rownames(input_methyl), custom_probes), rownames(FlowSorted.CordBlood.450k.compTable))
      res <- projectCellType(input_methyl[use_probes,],as.matrix(FlowSorted.CordBlood.450k.compTable[use_probes,3:9]))
    }
    return(res)
  }
  
  if (method == "Houseman" && tissue == "DLPFC"){
    library(FlowSorted.DLPFC.450k)
    source("projectCellType.R")
    if (is.null(custom_probes)){
      load("FlowSorted.DLPFC.450k.ModelPars.RData")
      use_probes <- interect(rownames(input_methyl), rownames(FlowSorted.DLPFC.450k.ModelPars))
      res <- projectCellType(input_methyl[use_probes,], 
                             as.matrix(FlowSorted.DLPFC.450k.ModelPars[use_probes,]))
    }else{
      load("FlowSorted.DLPFC.450k.compTable.RData")
      use_probes <- intersect(intersect(rownames(input_methyl), custom_probes), rownames(FlowSorted.DLPFC.450k.compTable))
      res <- projectCellType(input_methyl[use_probes,],as.matrix(FlowSorted.DLPFC.450k.compTable[use_probes,3:4]))
    }
    return(res)
  }
  
  if (method == "Houseman" && tissue == "Breast"){
    library(EpiDISH)
    source("projectCellType.R")
    if (is.null(custom_probes)){
      data("centEpiFibFatIC.m")
      use_probes <- intersect(rownames(input_methyl),rownames(centEpiFibFatIC.m))
      res <- projectCellType(input_methyl[use_probes,], 
                             centEpiFibFatIC.m[use_probes,])
    }else{
      stop("For breast tissue, only support EpiDISH selected probes (centEpiFibFatIC.m), 
           custom probes do not work!")
    }
    return(res)
  }
  
  if (method == "Houseman" && tissue == "genericEpithelial"){
    library(EpiDISH)
    source("projectCellType.R")
    if (is.null(custom_probes)){
      data("centEpiFibIC.m")
      use_probes <- intersect(rownames(input_methyl),rownames(centEpiFibIC.m))
      res <- projectCellType(input_methyl[use_probes,], 
                             centEpiFibIC.m[use_probes,])
      data("centDHSbloodDMC.m")
      res2 <- hepidish(input_methyl, centEpiFibIC.m, centDHSbloodDMC.m, h.CT.idx = 3, method = "CP")
    }else{
      stop("For generic epithelial tissue, only support EpiDISH selected probes (centEpiFibIC.m), 
           custom probes do not work!")
    }
    return(list(res,res2))
    }
  
  if (method == "RPC" && tissue == "Blood"){
    library(FlowSorted.Blood.450k)
    library(EpiDISH)
    if (is.null(custom_probes)){
      res <- epidish(input_methyl, as.matrix(FlowSorted.Blood.450k.JaffeModelPars), method = "RPC")
    }
    else{
      res <- epidish(input_methyl, as.matrix(FlowSorted.Blood.450k.compTable[cutom_probes,3:8]), method = "RPC")
    }
    res <- res$estF
    return(res)
  }
  
  if (method == "RPC" && tissue == "CordBlood"){
    library(FlowSorted.CordBlood.450k)
    library(EpiDISH)
    if (is.null(custom_probes)){
      res <- epidish(input_methyl, as.matrix(FlowSorted.CordBlood.450k.ModelPars), method = "RPC")
    }else{
      res <- epidish(input_methyl, as.matrix(FlowSorted.CordBlood.450k.compTable[cutom_probes,3:9]), method = "RPC")
    }
    res <- res$estF
    return(res)
  }
  
  if (method == "RPC" && tissue == "DLPFC"){
    library(FlowSorted.DLPFC.450k)
    library(EpiDISH)
    if (is.null(custom_probes)){
      load("FlowSorted.DLPFC.450k.ModelPars.RData")
      res <- epidish(input_methyl, as.matrix(FlowSorted.DLPFC.450k.ModelPars), method = "RPC")
    }else{
      load("FlowSorted.DLPFC.450k.compTable.RData")
      res <- epidish(input_methyl, as.matrix(FlowSorted.DLPFC.450k.compTable[cutom_probes,3:4]), method = "RPC")
    }
    res <- res$estF
    return(res)
  }
  
  if (method == "RPC" && tissue == "Breast"){
    library(EpiDISH)
    if (is.null(custom_probes)){
      data("centEpiFibFatIC.m")
      res <- epidish(input_methyl, centEpiFibFatIC.m, method = "RPC")
      data("centDHSbloodDMC.m")
      res2 <- hepidish(input_methyl, centEpiFibFatIC.m, centDHSbloodDMC.m, h.CT.idx = 4, method = "RPC")
    }else{
      stop("For breast tissue, only support EpiDISH selected probes (centEpiFibFatIC.m), 
           custom probes do not work!")
    }
    res <- res$estF
    return(list(res, res2))
  }
  
  if (method == "RPC" && tissue == "genericEpithelial"){
    library(EpiDISH)
    if (is.null(custom_probes)){
      data("centEpiFibIC.m")
      res <- epidish(input_methyl, centEpiFibIC.m, method = "RPC")
      data("centDHSbloodDMC.m")
      res2 <- hepidish(input_methyl, centEpiFibIC.m, centDHSbloodDMC.m, h.CT.idx = 3, method = "RPC")
    }else{
      stop("For generic epithelial tissue, only support EpiDISH selected probes (centEpiFibIC.m), 
           custom probes do not work!")
    }
    res <- res$estF
    return(list(res, res2))
    }
  
  if (method == "CBS" && tissue == "Blood"){
    library(FlowSorted.Blood.450k)
    library(EpiDISH)
    if (is.null(custom_probes)){
      res <- epidish(input_methyl, as.matrix(FlowSorted.Blood.450k.JaffeModelPars), method = "CBS")
    }
    else{
      res <- epidish(input_methyl, as.matrix(FlowSorted.Blood.450k.compTable[cutom_probes,3:8]), method = "CBS")
    }
    res <- res$estF
    return(res)
  }
  
  if (method == "CBS" && tissue == "CordBlood"){
    library(FlowSorted.CordBlood.450k)
    library(EpiDISH)
    if (is.null(custom_probes)){
      res <- epidish(input_methyl, as.matrix(FlowSorted.CordBlood.450k.ModelPars), method = "CBS")
    }else{
      res <- epidish(input_methyl, as.matrix(FlowSorted.CordBlood.450k.compTable[cutom_probes,3:9]), method = "CBS")
    }
    res <- res$estF
    return(res)
  }
  
  if (method == "CBS" && tissue == "DLPFC"){
    library(FlowSorted.DLPFC.450k)
    library(EpiDISH)
    if (is.null(custom_probes)){
      load("FlowSorted.DLPFC.450k.ModelPars.RData")
      res <- epidish(input_methyl, as.matrix(FlowSorted.DLPFC.450k.ModelPars), method = "CBS")
    }else{
      load("FlowSorted.DLPFC.450k.compTable.RData")
      res <- epidish(input_methyl, as.matrix(FlowSorted.DLPFC.450k.compTable[cutom_probes,3:4]), method = "CBS")
    }
    res <- res$estF
    return(res)
  }
  
  if (method == "CBS" && tissue == "Breast"){
    library(EpiDISH)
    if (is.null(custom_probes)){
      data("centEpiFibFatIC.m")
      res <- epidish(input_methyl, centEpiFibFatIC.m, method = "CBS")
      data("centDHSbloodDMC.m")
      res2 <- hepidish(input_methyl, centEpiFibFatIC.m, centDHSbloodDMC.m, h.CT.idx = 4, method = "CBS")
    }else{
      stop("For breast tissue, only support EpiDISH selected probes (centEpiFibFatIC.m), 
           custom probes do not work!")
    }
    res <- res$estF
    return(list(res,res2))
  }
  
  if (method == "CBS" && tissue == "genericEpithelial"){
    library(EpiDISH)
    if (is.null(custom_probes)){
      data("centEpiFibIC.m")
      res <- epidish(input_methyl, centEpiFibIC.m, method = "CBS")
      data("centDHSbloodDMC.m")
      res2 <- hepidish(input_methyl, centEpiFibIC.m, centDHSbloodDMC.m, h.CT.idx = 3, method = "CBS")
    }else{
      stop("For generic epithelial tissue, only support EpiDISH selected probes (centEpiFibIC.m), 
           custom probes do not work!")
    }
    res <- res$estF
    return(list(res,res2))
    }
}


MethylDeconv_BloodEPIC <- function(input_methyl, method = "Houseman", normalized = TRUE, custom_probes = NULL){
  if(normalized){
    # beta value input
    MethylDeconv_normalized_BloodEPIC(input_methyl, method, custom_probes)
  }
  else{
    ## input is RgSet
    library(minfi)
    library(FlowSorted.Blood.EPIC)
    data (IDOLOptimizedCpGs)
    
    RGset <- input_methyl
    .isRGOrStop(RGset)
    .isEPICOrStop(RGset)
    GRset.normalized <- preprocessNoob(RGset)
    
    output <- list()
    output[[1]] <- estimateCellCounts2(RGset, compositeCellType = "Blood", processMethod = "preprocessNoob",
                                       probeSelect = "IDOL", 
                                       cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                                       referencePlatform = "IlluminaHumanMethylationEPIC",
                                       referenceset = NULL,
                                       IDOLOptimizedCpGs =IDOLOptimizedCpGs,
                                       returnAll = FALSE)
    output[[2]] <- MethylDeconv_normalized_BloodEPIC(getBeta(GRset.normalized), method, custom_probes)
    return(output)
  }
}

MethylDeconv_normalized_BloodEPIC <- function(input_methyl, method = "Houseman", custom_probes = NULL){
  if (method == "Houseman"){
    source("projectCellType.R")
    if (is.null(custom_probes)){
      load("FlowSorted.Blood.EPIC.IDOLModelPars.RData")
      use_probes <- intersect(rownames(input_methyl), rownames(FlowSorted.Blood.EPIC.IDOLModelPars))
      res <- projectCellType(input_methyl[use_probes,], 
                             as.matrix(FlowSorted.Blood.EPIC.IDOLModelPars[use_probes,]))}
    else{
      load("Users/junesong/Desktop/causal inference/CellProportion/methylDeconv_EPICdata/FlowSorted.Blood.EPIC.compTable.RData")
      use_probes <- intersect(intersect(rownames(FlowSorted.Blood.EPIC.compTable), custom_probes),rownames(input_methyl))
      res <- projectCellType(input_methyl[use_probes,],as.matrix(FlowSorted.Blood.EPIC.compTable[use_probes,3:8]))
    }
    return(res)
  }
  
  if (method == "RPC"){
    library(EpiDISH)
    if (is.null(custom_probes)){
      load("FlowSorted.Blood.EPIC.IDOLModelPars.RData")
      res <- epidish(input_methyl, as.matrix(FlowSorted.Blood.EPIC.IDOLModelPar), method = "RPC")
    }
    else{
      load("Users/junesong/Desktop/causal inference/CellProportion/methylDeconv_EPICdata/FlowSorted.Blood.EPIC.compTable.RData")
      res <- epidish(input_methyl, as.matrix(FlowSorted.Blood.EPIC.compTable[cutom_probes,3:8]), method = "RPC")
    }
    res <- res$estF
    return(res)
  }
  
  if (method == "CBS"){
    library(EpiDISH)
    if (is.null(custom_probes)){
      load("FlowSorted.Blood.EPIC.IDOLModelPars.RData")
      res <- epidish(input_methyl, as.matrix(FlowSorted.Blood.EPIC.IDOLModelPar), method = "CBS")
    }
    else{
      load("Users/junesong/Desktop/causal inference/CellProportion/methylDeconv_EPICdata/FlowSorted.Blood.EPIC.compTable.RData")
      res <- epidish(input_methyl, as.matrix(FlowSorted.Blood.EPIC.compTable[cutom_probes,3:8]), method = "CBS")
    }
    res <- res$estF
    return(res)
  }
}















