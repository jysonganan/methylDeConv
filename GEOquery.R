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
# res1_1 <- MethylDeconv(rgSet_127824, method = "Houseman", normalized = FALSE, tissue = "CordBlood")
### check correlation
pD_127824 <- pD.all[, c("title", "geo_accession", "b cells:ch1", "cd4t cells:ch1", "cd8t cells:ch1", "granulocytes:ch1", 
                 "monocytes:ch1", "nk cells:ch1", "nrbcs:ch1", "Sex:ch1", "subject status:ch1", "tissue:ch1")]
pD_127824 <- pD_127824[sampleNames(rgSet_127824),]
facs_127824 <- pD_127824[,3:9]
facs_127824_prop <- apply(facs_127824,1,function(x){return(as.numeric(x)/sum(as.numeric(x)))})
facs_127824_prop <- t(facs_127824_prop)

rsquared <- rep(NA, 7)
for (i in 1:7){
  rsquared[i] <- summary(lm(res1[,i]~as.numeric(facs_127824[,i])))$r.squared
}
df <- data.frame(cellType = c("Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC"),
                 Rsquared = rsquared)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=Rsquared)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(rsquared,2)), vjust=-0.3, size=3.5)
p

corr <- rep(NA, 7)
for (i in 1:7){
  corr[i] <- cor(res1[,i], as.numeric(facs_127824_prop[,i]), method = "spearman")
}
df <- data.frame(cellType = c("Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC"),
                 SpearmanCorr = corr)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=SpearmanCorr)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(corr,2)), vjust=-0.3, size=3.5)
p


#
res2 <- MethylDeconv(betaMat_127824, method = "RPC", normalized = TRUE, tissue = "CordBlood")
rsquared <- rep(NA, 7)
for (i in 1:7){
  rsquared[i] <- summary(lm(res2[,i]~as.numeric(facs_127824[,i])))$r.squared
}
df <- data.frame(cellType = c("Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC"),
                 Rsquared = rsquared)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=Rsquared)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(rsquared,2)), vjust=-0.3, size=3.5)
p

corr <- rep(NA, 7)
for (i in 1:7){
  corr[i] <- cor(res2[,i], as.numeric(facs_127824_prop[,i]), method = "spearman")
}
df <- data.frame(cellType = c("Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC"),
                 SpearmanCorr = corr)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=SpearmanCorr)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(corr,2)), vjust=-0.3, size=3.5)
p

res3 <- MethylDeconv(betaMat_127824, method = "CBS", normalized = TRUE, tissue = "CordBlood")
rsquared <- rep(NA, 7)
for (i in 1:7){
  rsquared[i] <- summary(lm(res3[,i]~as.numeric(facs_127824[,i])))$r.squared
}
df <- data.frame(cellType = c("Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC"),
                 Rsquared = rsquared)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=Rsquared)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(rsquared,2)), vjust=-0.3, size=3.5)
p

corr <- rep(NA, 7)
for (i in 1:7){
  corr[i] <- cor(res3[,i], as.numeric(facs_127824_prop[,i]), method = "spearman")
}
df <- data.frame(cellType = c("Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC"),
                 SpearmanCorr = corr)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=SpearmanCorr)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(corr,2)), vjust=-0.3, size=3.5)
p



## 4. Methylation measured on EPIC platform whole blood, compared with FACS measured cell proportions.
getGEOSuppFiles("GSE112618")
untar("GSE112618/GSE112618_RAW.tar", exdir = "GSE112618/idat")
head(list.files("GSE112618/idat", pattern = "idat"))

idatFiles <- list.files("GSE112618/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)
library(minfi)
rgSet_112618 <- read.metharray.exp("GSE112618/idat")

geoMat_112618 <- getGEO("GSE112618")
pD.all <- pData(geoMat_112618[[1]])


pD <- pD.all[, c("title", "geo_accession", "bcell proportion:ch1", "cd4t proportion:ch1", "cd8t proportion:ch1",
                 "granulocytes proportion:ch1", "monocytes proportion:ch1", "neutrophils proportion:ch1", 
                 "nk proportion:ch1","Sex:ch1")]




## add phenotype data
sampleNames(rgSet_112618) <- substr(sampleNames(rgSet_112618), 1, 10)
pD <- pD[sampleNames(rgSet_112618),]
pD <- as(pD, "DataFrame")
pData(rgSet_112618) <- pD
rgSet_112618

# res1 <- MethylDeconv_BloodEPIC(rgSet_112618, method = "Houseman", normalized = FALSE)
# res2 <- MethylDeconv_BloodEPIC(rgSet_112618, method = "RPC", normalized = FALSE)
# res3 <- MethylDeconv_BloodEPIC(rgSet_112618, method = "CBS", normalized = FALSE)

##
grSet_112618 <- preprocessNoob(rgSet_112618)
getBeta(grSet_112618)[1:3,1:3]
head(getIslandStatus(grSet_112618))
betaMat_112618 <- getBeta(grSet_112618)
## 866091      6
res1 <- MethylDeconv_normalized_BloodEPIC(betaMat_112618, method = "Houseman")
res2 <- MethylDeconv_normalized_BloodEPIC(betaMat_112618, method = "RPC")
res3 <- MethylDeconv_normalized_BloodEPIC(betaMat_112618, method = "CBS")
###



pD_112618 <-pD.all[, c("title", "geo_accession", "bcell proportion:ch1", "cd4t proportion:ch1", "cd8t proportion:ch1",
                       "granulocytes proportion:ch1", "monocytes proportion:ch1", "neutrophils proportion:ch1", 
                       "nk proportion:ch1","Sex:ch1")]
pD_112618 <- pD_112618[sampleNames(rgSet_112618),]
facs_112618 <- pD_112618[,3:9]
colnames(facs_112618) <- c("Bcell", "CD4T", "CD8T","Gran","Mono","Neutro","NK")
facs_112618_prop <- cbind(facs_112618[,"CD8T"], facs_112618[,"CD4T"], facs_112618[,"NK"],
                          facs_112618[,"Bcell"],facs_112618[,"Mono"],
                          as.numeric(facs_112618[,"Gran"])+as.numeric(facs_112618[,"Neutro"]))
facs_112618_prop <- as.data.frame(facs_112618_prop)
rownames(facs_112618_prop) = rownames(res1)
colnames(facs_112618_prop) = colnames(res1)
facs_112618_prop

rsquared <- rep(NA, 6)
for (i in 1:6){
  rsquared[i] <- summary(lm(res1[,i]~as.numeric(facs_112618_prop[,i])))$r.squared
}
df <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"),
                 Rsquared = rsquared)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=Rsquared)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(rsquared,2)), vjust=-0.3, size=3.5)
p

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(res1[,i],as.numeric(facs_112618_prop[,i]),method = "spearman")
}
df <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"),
                 SpearmanCorr = corr)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=SpearmanCorr)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(corr,2)), vjust=-0.3, size=3.5)
p


rsquared <- rep(NA, 6)
for (i in 1:6){
  rsquared[i] <- summary(lm(res2[,i]~as.numeric(facs_112618_prop[,i])))$r.squared
}
df <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"),
                 Rsquared = rsquared)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=Rsquared)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(rsquared,2)), vjust=-0.3, size=3.5)
p

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(res2[,i],as.numeric(facs_112618_prop[,i]),method = "spearman")
}
df <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"),
                 SpearmanCorr = corr)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=SpearmanCorr)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(corr,2)), vjust=-0.3, size=3.5)
p



rsquared <- rep(NA, 6)
for (i in 1:6){
  rsquared[i] <- summary(lm(res3[,i]~as.numeric(facs_112618_prop[,i])))$r.squared
}
df <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"),
                 Rsquared = rsquared)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=Rsquared)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(rsquared,2)), vjust=-0.3, size=3.5)
p

corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(res3[,i],as.numeric(facs_112618_prop[,i]),method = "spearman")
}
df <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"),
                 SpearmanCorr = corr)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=SpearmanCorr)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(corr,2)), vjust=-0.3, size=3.5)
p



## 5 GSE110530
# Methylation measured on EPIC platform, compared with FACS measured cell proportions
getGEOSuppFiles("GSE110530")
untar("GSE110530/GSE110530_RAW.tar", exdir = "GSE110530/idat")
head(list.files("GSE110530/idat", pattern = "idat"))

idatFiles <- list.files("GSE110530/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)
library(minfi)
rgSet_110530 <- read.metharray.exp("GSE110530/idat")

geoMat_110530 <- getGEO("GSE110530")
pD.all <- pData(geoMat_110530[[1]])

## missing FACS counts for 5 samples of the total 12 samples



## 6 GSE69914
# 450k, whole blood, purified breast epithelial cells.
### limit reached
geoMat_69914<- getGEO("GSE69914")
pD.all <- pData(geoMat_69914[[1]])


##7 GSE67919
geoMat_67919<- getGEO("GSE67919")

## 8 