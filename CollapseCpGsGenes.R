# collapse CpGs into genes

## correct version TCGA KICH
CollapseCpGsGenes <- function(BetaMatrix, method = "average"){
  manifest <- read.csv("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/HumanMethylation450_15017482_v1-2.csv", header = T, skip = 7) #486428     33
  #manifest <- read.csv("HumanMethylation450_15017482_v1-2.csv", header = T, skip = 7)
  annot <- manifest[match(rownames(BetaMatrix),manifest[,1]),]
  annot_df <- as.data.frame(annot)
  annot_df_oneGene <- annot_df
  rownames(annot_df_oneGene) <- annot_df[,1]
  annot_df_oneGene[,"UCSC_RefGene_Name"] <- gsub(";.*$","",annot_df[,"UCSC_RefGene_Name"])
  annot_df_oneGene <- annot_df_oneGene[,c("CHR", "Infinium_Design_Type", "Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Name")]
  annot_df_oneGene <- annot_df_oneGene[!(annot_df_oneGene[,"UCSC_RefGene_Name"]== ""),]
  BetaMatrix <- BetaMatrix[match(rownames(annot_df_oneGene),rownames(BetaMatrix)),]
  
  if (method == "maxvar"){
    var_per_CpG <- apply(BetaMatrix,1,var)
    var_per_CpG <- as.data.frame(var_per_CpG)
    rownames(var_per_CpG) <- rownames(BetaMatrix)
    
    gene_list <- unique(annot_df_oneGene[,"UCSC_RefGene_Name"])
    CpG_maxVar_list <- rep(NA, length(gene_list))
    
    for (i in 1:length(gene_list)){
      gene <- gene_list[i]
      CpGs <- rownames(annot_df_oneGene[annot_df_oneGene[,"UCSC_RefGene_Name"] == gene,])
      if (length(CpGs) == 1){
        CpG_maxVar_list[i] <- CpGs
      }
      else{
        var_per_CpG_sub <- var_per_CpG[CpGs,1]
        if(sum(is.na(var_per_CpG_sub))==length(CpGs)){
          CpG_maxVar_list[i] <- NA
        }
        else{
          CpG_maxVar_list[i] <- CpGs[which.max(var_per_CpG_sub)]
        }
      }
    }
    miss_id <- which(is.na(CpG_maxVar_list) == TRUE)
    CpG_maxVar_list <- CpG_maxVar_list[-miss_id]
    gene_list <- gene_list[-miss_id]
    
    genelevel_BetaMatrix <- BetaMatrix[CpG_maxVar_list,]
    rownames(genelevel_BetaMatrix) <- gene_list
    genelevel_BetaMatrix <- genelevel_BetaMatrix[complete.cases(genelevel_BetaMatrix),]
    return(genelevel_BetaMatrix)
  }
  
  if (method == "average"){
    BetaMatrix <- cbind(BetaMatrix, annot_df_oneGene[,"UCSC_RefGene_Name"])
    BetaMatrix <- as.data.frame(BetaMatrix)
    BetaMatrix[,-dim(BetaMatrix)[2]] <- apply(BetaMatrix[,-dim(BetaMatrix)[2]],2,as.numeric)
    colnames(BetaMatrix)[dim(BetaMatrix)[2]] <- "UCSC_RefGene_Name"
    
    aggdata <-aggregate(BetaMatrix[,-dim(BetaMatrix)[2]], by= list(BetaMatrix$UCSC_RefGene_Name), FUN=mean, na.rm=TRUE)
    rownames(aggdata) <- aggdata[,1]
    aggdata <- aggdata[,-1]
    aggdata <- aggdata[complete.cases(aggdata),]
    return(aggdata)
    }
  
  if (method == "PCA"){
    BetaMatrix <- cbind(BetaMatrix, annot_df_oneGene[,"UCSC_RefGene_Name"])
    BetaMatrix <- as.data.frame(BetaMatrix)
    BetaMatrix[,-dim(BetaMatrix)[2]] <- apply(BetaMatrix[,-dim(BetaMatrix)[2]],2,as.numeric)
    colnames(BetaMatrix)[dim(BetaMatrix)[2]] <- "UCSC_RefGene_Name"
    no_na <- apply(BetaMatrix,1,function(x){return(sum(is.na(x)))})
    # delete probes with all missing values
    BetaMatrix <- BetaMatrix[which(no_na < dim(BetaMatrix)[2]-1),]
    # impute the remaining missing values

    BetaMatrix[,1:66] <- t(apply(BetaMatrix[,1:66],1,function(x){x[is.na(x)] <- median(x, na.rm = TRUE);return(x)}))
    BetaMatrix[,67] <- as.character(BetaMatrix[,67])
    pcaCollapse <- function(x){
      if(dim(x)[1] == 1){
        return(x)
      }
      else{
        #x <- apply(x,2,function(x){return(as.numeric(x))})
        return(prcomp(t(x), scale. = TRUE)$x[,1])
      }
    }
    BetaMatrix <- as.data.frame(BetaMatrix)
    res <- sapply(split(BetaMatrix[,1:66], BetaMatrix$UCSC_RefGene_Name), pcaCollapse)
    res <- t(res)
    gene_names <- rownames(res)
    res <- apply(res,2,function(x){return(as.numeric(x))})
    rownames(res) <- gene_names
    return(res)
  }
}

  
  

  


















  
    

