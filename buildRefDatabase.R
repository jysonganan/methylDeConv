


### 450k blood immune cells   

### FlowSorted.Blood.450k,  (Reinius 2012), GSE35069
library(FlowSorted.Blood.450k)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]

#1.  GSE43976, 43 samples,Mono
library(GEOquery)
library(minfi)
getGEOSuppFiles("GSE43976")
untar("GSE43976/GSE43976_RAW.tar", exdir = "GSE43976/idat")
head(list.files("GSE43976/idat", pattern = "idat"))

idatFiles <- list.files("GSE43976/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

geoMat <- getGEO("GSE43976")
pD.all <- pData(geoMat[[1]])
pD.all <- pD.all[pD.all[,"source_name_ch1"]=="CD14+ monocytes",]


rgSet <- read.metharray.exp("GSE43976/idat", force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
colnames(betaMat) <- substr(colnames(betaMat),1,10)
betaMat_43976 <- betaMat[,rownames(pD.all)]

phenotype_43976 <- rep("Mono", dim(pD.all)[1])
save("betaMat_43976","phenotype_43976", file = "ref_43976_450kBlood.RData")




#2.  GSE59065  99 CD4T, 100 CD8T
getGEOSuppFiles("GSE59065")
untar("GSE59065/GSE59065_RAW.tar", exdir = "GSE59065/idat")
head(list.files("GSE59065/idat", pattern = "idat"))

idatFiles <- list.files("GSE59065/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

geoMat <- getGEO("GSE59065")
pD.all <- pData(geoMat[[1]])
pD.all <- pD.all[pD.all[,"cell/tissue type:ch1"]%in%c("CD4","CD8"),]


rgSet <- read.metharray.exp("GSE59065/idat", force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
colnames(betaMat) <- substr(colnames(betaMat),1,10)
betaMat_59065 <- betaMat[,rownames(pD.all)]

phenotype_59065 <- pD.all[,"cell/tissue type:ch1"]
phenotype_59065 <- paste0(phenotype_59065,'T')
save("betaMat_59065","phenotype_59065", file = "ref_59065_450kBlood.RData")



#3  GSE71955 67 CD4T  68 CD8T
getGEOSuppFiles("GSE71955")
untar("GSE71955/GSE71955_RAW.tar", exdir = "GSE71955/idat")
head(list.files("GSE71955/idat", pattern = "idat"))

idatFiles <- list.files("GSE71955/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

geoMat <- getGEO("GSE71955")
pD.all <- pData(geoMat[[1]])
pD.all <- pD.all[pD.all[,"cell type:ch1"]%in%c("CD4 T cells","CD8 T cells"),]
######

rgSet <- read.metharray.exp("GSE71955/idat", force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
colnames(betaMat) <- substr(colnames(betaMat),1,10)
betaMat_71955 <- betaMat[,rownames(pD.all)]

phenotype_71955 <- pD.all[,"cell type:ch1"]
for (i in 1:length(phenotype_71955)){
  if (phenotype_71955[i] == "CD8 T cells"){
    phenotype_71955[i] <- "CD8T"
  }
  
  if (phenotype_71955[i] == "CD4 T cells"){
    phenotype_71955[i] <- "CD4T"
  }
}
save("betaMat_71955","phenotype_71955", file = "ref_71955_450kBlood.RData")



# 4 GSE59250    308 samples            CD14+ Monocytes 55                 CD19+ B-Cells 104
#          CD4+ T-cells 149


## there are NAs and may not preprocessNoob
geoMat <- getGEO("GSE59250")
pD.all <- pData(geoMat[[1]])
pD.all[,"cell type:ch1"]

getGEOSuppFiles("GSE59250")
untar("GSE59250/GSE59250_RAW.tar", exdir = "GSE59250/idat")
#head(list.files("GSE59250/idat", pattern = "idat"))

betaMat <- read.table("GSE59250/GSE59250_average_betas.txt", header = T, sep = "\t")
rownames(betaMat) <- betaMat[,1]
betaMat <- betaMat[,-1]
betaMat <- betaMat[rownames(ref_betamatrix),]
betaMat <- betaMat[,seq(1,867,2)]
colnames(betaMat) <- rownames(pD.all)

pD.all <- pD.all[pD.all[,"cell type:ch1"]%in%c("CD14+ Monocytes","CD19+ B-Cells","CD4+ T-cells"),]
betaMat_59250 <- betaMat[,rownames(pD.all)]

phenotype_59250 <- pD.all[,"cell type:ch1"]
for (i in 1:length(phenotype_59250)){
  if (phenotype_59250[i] == "CD4+ T-cells"){
    phenotype_59250[i] <- "CD4T"
  }
  
  if (phenotype_59250[i] == "CD19+ B-Cells"){
    phenotype_59250[i] <- "Bcell"
  }
  
  if (phenotype_59250[i] == "CD14+ Monocytes"){
    phenotype_59250[i] <- "Mono"
  }
}
save("betaMat_59250","phenotype_59250", file = "ref_59250_450kBlood.RData")





# GSE71244 # 20 samples 4 Bcell, 6 cd4t, 5 cd8t, 5 mono     ## normalized, may not be noob preprocessed. there are NAs!
getGEOSuppFiles("GSE71244")
untar("GSE71244/GSE71244_RAW.tar", exdir = "GSE71244/idat")
#head(list.files("GSE71244/idat", pattern = "idat"))
betaMat <- read.table("GSE71244_series_matrix.txt", skip = 73, header = T, sep = "\t", nrows = 485577)
rownames(betaMat) <- betaMat[,1]
betaMat <- betaMat[,-1]
betaMat <- betaMat[rownames(ref_betamatrix),]

geoMat <- getGEO("GSE71244")
pD.all <- pData(geoMat[[1]])
pD.all <- pD.all[pD.all[,"cell subset:ch1"]%in%c("Monocyte","B cells","CD8 T cells","CD4 T cells"),]

betaMat_71244 <- betaMat[,rownames(pD.all)]

phenotype_71244<- pD.all[,"cell subset:ch1"]
for (i in 1:length(phenotype_71244)){
  if (phenotype_71244[i] == "Monocyte"){
    phenotype_71244[i] <- "Mono"
  }
  
  if (phenotype_71244[i] == "B cells"){
    phenotype_71244[i] <- "Bcell"
  }
  
  if (phenotype_71244[i] == "CD8 T cells"){
    phenotype_71244[i] <- "CD8T"
  }
  
  if (phenotype_71244[i] == "CD4 T cells"){
    phenotype_71244[i] <- "CD4T"
  }
}
save("betaMat_71244","phenotype_71244", file = "ref_71244_450kBlood.RData")




#GSE50222   #32 samples, CD4T  #### but may not be preprocessNoob!! And many null probes!!
# getGEOSuppFiles("GSE50222")
# untar("GSE50222/GSE50222_RAW.tar", exdir = "GSE50222/idat")
# head(list.files("GSE50222/idat", pattern = "idat"))
betaMat <- read.table("GSE50222_series_matrix.txt", skip = 62, header = T, sep = "\t", nrows = 485577)
rownames(betaMat) <- betaMat[,1]
betaMat <- betaMat[,-1]
betaMat <- betaMat[rownames(ref_betamatrix),]

geoMat <- getGEO("GSE50222")
pD.all <- pData(geoMat[[1]])
pD.all <- pD.all[pD.all[,"source_name_ch1"]=="CD4+ T-cells",]

betaMat_50222 <- betaMat[,rownames(pD.all)]

phenotype_50222<- rep("CD4T", dim(pD.all)[1])
save("betaMat_50222","phenotype_50222", file = "ref_50222_450kBlood.RData")






#GSE56046 Mono 1202 samples  ### no good files for betamat
getGEOSuppFiles("GSE56046")
untar("GSE56046/GSE56046_RAW.tar", exdir = "GSE56046/idat")
#head(list.files("GSE56046/idat", pattern = "idat"))

betaMat <- read.table("GSE56046/GSE56046_methylome_normalized.txt", header = T, sep = "\t")
rownames(betaMat) <- betaMat[,1]
betaMat <- betaMat[,-1]
betaMat <- betaMat[rownames(ref_betamatrix),]

geoMat <- getGEO("GSE56046")
pD.all <- pData(geoMat[[1]])
pD.all <- pD.all[pD.all[,"cell subset:ch1"]%in%c("Monocyte","B cells","CD8 T cells","CD4 T cells"),]

betaMat_71244 <- betaMat[,rownames(pD.all)]

#GSE56581 T cells
getGEOSuppFiles("GSE56581")
untar("GSE56581/GSE56581_RAW.tar", exdir = "GSE56581/idat")




#################################### 
###### 450k fib, epithelial cells
####################################

# # ??? GSE31848   4 epithelial cell lines, 10 fibroblasts cell lines 
# geoMat <- getGEO("GSE31848")
# pD.all <- pData(geoMat[[1]])
# pD.all <- pD.all[pD.all[,"cell/tissue type:ch1"]%in%c("CD4","CD8"),]

### GSE122126 450K
##  Adipocytes 3; cfDNA: 1; Cortical neurons: 1; Hepatocytes: 1; Pancreatic acinar cells: 1; Pancreatic beta cells: 3; Pancreatic duct cells: 1
geoMat <- getGEO("GSE122126")
pD.all <- pData(geoMat[[1]])

getGEOSuppFiles("GSE122126")
untar("GSE122126/GSE122126_RAW.tar", exdir = "GSE122126/idat")
head(list.files("GSE122126/idat", pattern = "idat"))

idatFiles <- list.files("GSE122126/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)


rgSet <- read.metharray.exp("GSE122126/idat/450k", force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)

colnames(betaMat) <- substr(colnames(betaMat),1,10)
betaMat_122126 <- betaMat[,rownames(pD.all)]

phenotype_122126 <- pD.all[,"sample type:ch1"]

save("betaMat_122126","phenotype_122126", file = "ref_122126_450kEpithelial.RData")













#################################### 
###### 450k Neuron
####################################
library(FlowSorted.DLPFC.450k)
#library(minfi)
GRset_frontCortex <-  preprocessNoob(FlowSorted.DLPFC.450k, dyeMethod = "single")
CellLines.matrix = NULL
cellTypes = c("NeuN_neg","NeuN_pos")
ref_betamatrix <- getBeta(GRset_frontCortex)
ref_phenotype <- as.data.frame(colData(FlowSorted.DLPFC.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]


#################################### 
###### EPIC Neuron
####################################
data_type = "BrainEPIC"
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



















#################################### 
###### EPIC blood
####################################
## make FlowSorted.Blood.EPIC.RData
# library(ExperimentHub)  
# hub <- ExperimentHub()  
# query(hub, "FlowSorted.Blood.EPIC")  
# FlowSorted.Blood.EPIC <- hub[["EH1136"]]  
# FlowSorted.Blood.EPIC  

#Bcell  CD4T  CD8T  Mono   Neu    NK 
#6     7     6     6     6     6 
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
load("FlowSorted.Blood.EPIC.RData")
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.EPIC, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.EPIC))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]
#450 CpGs of six immune cell subtypes: neutrophils, B cells, monocytes, NK cells, CD4+ T cells, and CD 8+ T cells
# Consists of 37 magnetic sorted blood cell references and 12 artificial mixture samples.



#### EPIC benchmark with true proportions
####################################

###!!! 12 mixture with known proportions can be used as benchmark
load("FlowSorted.Blood.EPIC.RData")
annot <- as.data.frame(colData(FlowSorted.Blood.EPIC))
benchmark <- which(annot$CellType == "MIX")
tmp <- getBeta(preprocessNoob(FlowSorted.Blood.EPIC, dyeMethod = "single"))
benchmark_betamatrix <- tmp[,rownames(annot)[benchmark]]
benchmark_trueprop <- annot[benchmark, c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]



### GSE112618 6 samples of known proportions  ## maybe problematic as the sum wasn't 1.
geoMat <- getGEO("GSE112618")
pD.all <- pData(geoMat[[1]])
pD <- cbind(as.numeric(pD.all[,"bcell proportion:ch1"]), as.numeric(pD.all[,"cd4t proportion:ch1"]), 
                          as.numeric(pD.all[,"cd8t proportion:ch1"]), as.numeric(pD.all[,"monocytes proportion:ch1"]),
                          as.numeric(pD.all[,"neutrophils proportion:ch1"]),as.numeric(pD.all[,"nk proportion:ch1"]))
pD <- as.data.frame(pD)
colnames(pD) <- c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")
rownames(pD) <- rownames(pD.all)

getGEOSuppFiles("GSE112618")
untar("GSE112618/GSE112618_RAW.tar", exdir = "GSE112618/idat")
head(list.files("GSE112618/idat", pattern = "idat"))

idatFiles <- list.files("GSE112618/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

rgSet_112618 <- read.metharray.exp("GSE112618/idat",force = TRUE)

benchmark_betamatrix <- getBeta(preprocessNoob(rgSet_112618, dyeMethod = "single"))
benchmark_trueprop <- pD




#################################### 
###### EPIC epithelial
####################################

### 2. GSE122126 EPIC

#Purified pancreatic acinar cells, pancreatic duct cells, 
# pancreatic beta cells, vascular endothelial cells and colon epithelial cells of 450k and EPIC data.
# 
# cfDNA         cfDNA In vitro mix 
# 58                          5 
# **Colon epithelial cells           Cortical neurons 
# 3                          2 
# Hepatocytes               In vitro mix 
# 2                          9 
# Leukocytes      **Lung epithelial cells 
# 1                          3 
# **Pancreatic acinar cells      Pancreatic beta cells 
# 2                          1 
# ** Pancreatic duct cells Vascular endothelial cells 
# 2                          2 

library(GEOquery)
library(minfi)
getGEOSuppFiles("GSE122126")
untar("GSE122126/GSE122126_RAW.tar", exdir = "GSE122126/idat")
head(list.files("GSE122126/idat", pattern = "idat"))

idatFiles <- list.files("GSE122126/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)
rgSet <- read.metharray.exp("GSE122126/idat", force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
#866091     90

geoMat <- getGEO("GSE122126")
pD.all <- pData(geoMat[[2]])


colnames(betaMat) <- substr(colnames(betaMat),1,10)
betaMat_122126 <- betaMat[,rownames(pD.all)]

phenotype_122126 <- pD.all[,"sample type:ch1"]

save("betaMat_122126","phenotype_122126", file = "ref_122126_EPICEpithelial.RData")
















