MethylDeconv <- function(input_methyl, method = "Houseman", normalized = TRUE, tissue = "Blood", custom_probes = NULL){
  if(normalized){
    # beta value input
    MethylDeconv_normalized(input_methyl, method, tissue, custom_probes)
  }
  else{
    ## input is RgSet
    library(minfi)
    RGset <- input_methyl
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
      res <- projectCellType(input_methyl[rownames(FlowSorted.Blood.450k.JaffeModelPars),], 
                             as.matrix(FlowSorted.Blood.450k.JaffeModelPars))
    }else{
      res <- projectCellType(input_methyl[custom_probes,],as.matrix(FlowSorted.Blood.450k.compTable[cutom_probes,3:8]))
    }
    return(res)
  }
  
  if (method == "Houseman" && tissue == "CordBlood"){
    library(FlowSorted.CordBlood.450k)
    source("projectCellType.R")
    if (is.null(custom_probes)){
      res <- projectCellType(input_methyl[rownames(FlowSorted.CordBlood.450k.ModelPars),], 
                             as.matrix(FlowSorted.CordBlood.450k.ModelPars))
    }else{
      res <- projectCellType(input_methyl[custom_probes,],as.matrix(FlowSorted.CordBlood.450k.compTable[cutom_probes,3:9]))
    }
    return(res)
  }
}