## 40 samples from individuals hospitalized with acute mania as well as unaffected controls.

library(GEOquery)


#################################################################################
# supplemental files to GEO series (GSE), GEO platforms (GPL) and GEO samples (GSM)
getGEOSuppFiles("GSE68777")
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
head(list.files("GSE68777/idat", pattern = "idat"))

idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

library(minfi)
#targets <- read.metharray.sheet("GSE68777/idat") #however, no sample sheet!
rgSet <- read.metharray.exp("GSE68777/idat")
# 450k, (622399,40)

### phenotype data stored in GEO
geoMat <- getGEO("GSE68777")
pD.all <- pData(geoMat[[1]])
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)

names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis: ", "", pD$group)
pD$sex <- sub("^Sex: ", "", pD$sex)

## add phenotype data
sampleNames(rgSet) <- sub(".*_5", "5", sampleNames(rgSet))

rownames(pD) <- pD$title
pD <- pD[sampleNames(rgSet),]
pD <- as(pD, "DataFrame")
pData(rgSet) <- pD
rgSet

##
grSet <- preprocessQuantile(rgSet)
getBeta(grSet)[1:3,1:3]
head(getIslandStatus(grSet))
betaMat <- getBeta(grSet)

########################################################################################
geoMat <- getGEO("GSE68777")
class(exprs(geoMat[[1]]))
dim(exprs(geoMat[[1]]))
exprsMat <- exprs(geoMat[[1]])
colnames(exprsMat) <- colnames(betaMat)
exprsMat <- exprsMat[rownames(betaMat),]






### example
# (1) phenotype information
gse <- getGEO("GSE781",GSEMatrix=FALSE)
head(Meta(gse))
names(GSMList(gse))
GSMList(gse)[[1]]
names(GPLList(gse))

# (2) getting GSE series matrix as an ExpressionSet
## phenotype information
gse2553 <- getGEO('GSE2553', GSEMatrix = TRUE)
# gse_2553 <- getGEO('GSE2553') the same!!
class(gse2553[[1]])
dim(gse2553[[1]])
show(gse2553)
show(pData(phenoData(gse2553[[1]]))[1:5, c(1,6,8)])
class(exprs(gse2553[[1]]))
dim(exprs(gse2553[[1]])) #[1] 12600   181


### experiment data
## 1. FACS of neuronal nuclei and 450k in post mortem frontal cortex
##    of 29 major depression subjects and 29 matched controls.
geoMat_15014 <- getGEO("GSE15014")

pD.all <- pData(geoMat_15014[[1]])
pD <- pD.all[, c("marker:ch1","marker:ch2")]
head(pD)

## add phenotype data
getGEOSuppFiles("GSE15014")
untar("GSE15014/GSE15014_RAW.tar", exdir = "GSE15014/idat")
head(list.files("GSE15014/idat", pattern = "idat"))



## 2 whole blood for over 650 samples in RA cases and control (Liu data)
## geoMat_42861 <- getGEO("GSE42861")
## limit reached


##3 Test data of seven cell types (Gran, Mono, Bcell, CD4T, CD8T, NK and nRBCs) with matched FACS counts.
getGEOSuppFiles("GSE127824")
untar("GSE127824/GSE127824_RAW.tar", exdir = "GSE127824/idat")
head(list.files("GSE127824/idat", pattern = "idat"))

idatFiles <- list.files("GSE127824/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

library(minfi)
rgSet_127824 <- read.metharray.exp("GSE127824/idat")

geoMat_127824 <- getGEO("GSE127824")
pD.all <- pData(geoMat_127824[[1]])
pD <- pD.all[, c("title", "geo_accession", "b cells:ch1", "cd4t cells:ch1", "cd8t cells:ch1", "granulocytes:ch1", 
  "monocytes:ch1", "nk cells:ch1", "nrbcs:ch1", "Sex:ch1", "subject status:ch1", "tissue:ch1")]




## add phenotype data
sampleNames(rgSet_127824) <- substr(sampleNames(rgSet_127824), 1, 10)
pD <- pD[sampleNames(rgSet_127824),]
pD <- as(pD, "DataFrame")
pData(rgSet_127824) <- pD
rgSet_127824

##
grSet_127824 <- preprocessNoob(rgSet_127824)
getBeta(grSet_127824)[1:3,1:3]
head(getIslandStatus(grSet_127824))
betaMat_127824 <- getBeta(grSet_127824)
## 485512     24

res1 <- MethylDeconv(betaMat_127824, method = "Houseman", normalized = TRUE, tissue = "CordBlood")
### check correlation
pD_127824 <- pD.all[, c("title", "geo_accession", "b cells:ch1", "cd4t cells:ch1", "cd8t cells:ch1", "granulocytes:ch1", 
                 "monocytes:ch1", "nk cells:ch1", "nrbcs:ch1", "Sex:ch1", "subject status:ch1", "tissue:ch1")]
pD_127824 <- pD_127824[sampleNames(rgSet_127824),]
facs_127824 <- pD_127824[,3:9]
facs_127824_prop <- apply(facs_127824,1,function(x){return(as.numeric(x)/sum(as.numeric(x)))})
facs_127824_prop <- t(facs_127824_prop)
