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
geoMat_42861 <- getGEO("GSE42861")
pD.all <- pData(geoMat_42861[[1]])

getGEOSuppFiles("GSE42861")
untar("GSE42861/GSE42861_RAW.tar", exdir = "GSE42861/idat")
head(list.files("GSE42861/idat", pattern = "idat"))

idatFiles <- list.files("GSE42861/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)
## too large !




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
for (i in 1:7){
  facs_127824[,i] <- as.numeric(facs_127824[,i])
}
facs_127824_prop <- apply(facs_127824,1,function(x){return(as.numeric(x)/sum(as.numeric(x)))})
facs_127824_prop <- t(facs_127824_prop)

rsquared <- rep(NA, 7)
for (i in 1:7){
  rsquared[i] <- summary(lm(res1[,i]~facs_127824_prop[,i]))$r.squared
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
  corr[i] <- cor(res1[,i], facs_127824_prop[,i], method = "spearman")
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
  rsquared[i] <- summary(lm(res2[,i]~facs_127824[,i]))$r.squared
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
  rsquared[i] <- summary(lm(res1[,i]~as.numeric(as.character(facs_112618_prop[,i]))))$r.squared
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
  corr[i] <-cor(res1[,i],as.numeric(as.character(facs_112618_prop[,i])),method = "spearman")
}
df <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"),
                 SpearmanCorr = corr)
library(ggplot2)
p<-ggplot(data=df, aes(x=cellType, y=SpearmanCorr)) +
  geom_bar(stat="identity",width=0.5)+
  geom_text(aes(label= round(corr,2)), vjust=-0.3, size=3.5)
p



library(tidyr)
library(ggplot2)
df1 <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"), Rsquared = rsquared,
                 SpearmanCorr = corr)
ggplot(data = df1 %>% gather(Variable, values, -cellType), 
       aes(x = cellType, y = values, fill = Variable)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(values,2)), position = position_dodge(0.9))



rsquared <- rep(NA, 6)
for (i in 1:6){
  rsquared[i] <- summary(lm(res2[,i]~as.numeric(as.character(facs_112618_prop[,i]))))$r.squared
}
corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(res2[,i],as.numeric(as.character(facs_112618_prop[,i])),method = "spearman")
}
df1 <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"), Rsquared = rsquared,
                  SpearmanCorr = corr)
ggplot(data = df1 %>% gather(Variable, values, -cellType), 
       aes(x = cellType, y = values, fill = Variable)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(values,2)), position = position_dodge(0.9))




rsquared <- rep(NA, 6)
for (i in 1:6){
  rsquared[i] <- summary(lm(res3[,i]~as.numeric(as.character(facs_112618_prop[,i]))))$r.squared
}
corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(res3[,i],as.numeric(as.character(facs_112618_prop[,i])),method = "spearman")
}
df1 <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"), Rsquared = rsquared,
                  SpearmanCorr = corr)
ggplot(data = df1 %>% gather(Variable, values, -cellType), 
       aes(x = cellType, y = values, fill = Variable)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(values,2)), position = position_dodge(0.9))








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
geoMat_69914<- getGEO("GSE69914")
pD.all <- pData(geoMat_69914[[1]])
# 407 47
pD <- pD.all[, c("title", "geo_accession", 
                 "status(0=normal 1=normal-adjacent 2=breast cancer 3=normal-brca1 4=cancer-brca1):ch1")]

#getGEOSuppFiles("GSE69914")
#untar("GSE69914/GSE69914_RAW.tar", exdir = "GSE69914/idat")
## head(list.files("GSE69914/idat", pattern = "idat"))
## not idat.gz files
betaMat_69914 <- read.table("GSE69914_series_matrix.txt", skip = 73, header = T, sep = "\t", nrows = 485577)
rownames(betaMat_69914) <- betaMat_69914[,1]
betaMat_69914 <- betaMat_69914[,-1]

###





## 7 GSE67919 #96 samples
# 450k, breast epithelial cells with race, sex, age, smoking information.
geoMat_67919<- getGEO("GSE67919")
pD.all <- pData(geoMat_67919[[1]])
pD <- pD.all[, c("title", "geo_accession", "age at surgery:ch1",
                 "alcohol use:ch1","gender:ch1","menopausal status:ch1","race:ch1","tissue type:ch1")]
# 96 49

getGEOSuppFiles("GSE67919")
untar("GSE67919/GSE67919_RAW.tar", exdir = "GSE67919/idat")
head(list.files("GSE67919/idat", pattern = "idat"))
## not idat.gz file, only GPL file
betaMat_67919 <- read.table("GSE67919_series_matrix.txt", skip = 70, header = T, sep = "\t", nrows = 485577)
rownames(betaMat_67919) <- betaMat_67919[,1]
betaMat_67919 <- betaMat_67919[,-1]

res1 <- MethylDeconv(betaMat_67919, tissue = "Breast", method = "Houseman")
res2 <- MethylDeconv(betaMat_67919, tissue = "Breast", method = "RPC")
#res3 <- MethylDeconv(betaMat_67919, tissue = "Breast", method = "CBS")
tissue_dat_blood <- matrix(NA, 96, 5)
samples <- rownames(res1)
rownames(tissue_dat_blood) <- samples
tissue_dat_blood[,1:4] <- res1
tissue_dat_blood[,5] <- as.character(pD[,"tissue type:ch1"])
tissue_dat_blood <- as.data.frame(tissue_dat_blood)
tissue_dat_blood[1:4] <- apply(tissue_dat_blood[1:4], 2, as.numeric)
tissue_dat_blood[,5] <- as.character(tissue_dat_blood[,5])

colnames(tissue_dat_blood) <- c(colnames(res1),"tissue")

library(ggplot2)
library(tidyr)
df <- gather(tissue_dat_blood, series,value,-tissue)
ggplot(df) + geom_boxplot(aes(series ,value,color= tissue)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE67919-Breast-Houseman")



#res1 <- MethylDeconv(betaMat_67919, tissue = "genericEpithelial", method = "Houseman")
res2 <- MethylDeconv(betaMat_67919, tissue = "genericEpithelial", method = "RPC")
res2 <-res2[[2]]
#res3 <- MethylDeconv(betaMat_67919, tissue = "Breast", method = "CBS")
tissue_dat_blood <- matrix(NA, 96, 10)
samples <- rownames(res2)
rownames(tissue_dat_blood) <- samples
tissue_dat_blood[,1:9] <- res2
tissue_dat_blood[,10] <- as.character(pD[,"tissue type:ch1"])
tissue_dat_blood <- as.data.frame(tissue_dat_blood)
tissue_dat_blood[1:9] <- apply(tissue_dat_blood[1:9], 2, as.numeric)
tissue_dat_blood[,10] <- as.character(tissue_dat_blood[,10])

colnames(tissue_dat_blood) <- c(colnames(res2),"tissue")

library(ggplot2)
library(tidyr)
df <- gather(tissue_dat_blood, series,value,-tissue)
ggplot(df) + geom_boxplot(aes(series ,value,color= tissue)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE67919-genericEpithelial-RPC")







## 8 GSE77797 whole blood 18 samples
# 12 artificial reconstructed mixtures of purified leukocyte subtypes measured on 450k platform 
# with mixture proportions from a six-component Dirichlet distribution.
geoMat_77797<- getGEO("GSE77797")
pD.all <- pData(geoMat_77797[[1]])
# 18 47
getGEOSuppFiles("GSE77797")
untar("GSE77797/GSE77797_RAW.tar", exdir = "GSE77797/idat")
head(list.files("GSE77797/idat", pattern = "idat"))

idatFiles <- list.files("GSE77797/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

library(minfi)
rgSet_77797 <- read.metharray.exp("GSE77797/idat")

pD <- pD.all[, c("title", "geo_accession", "b cell (%):ch1", "cd4+ t cell (%):ch1", "cd8+ t cell (%):ch1",
                 "granulocyte (%):ch1", "monocyte (%):ch1", "natural killer cell (%):ch1")]


## add phenotype data
sampleNames(rgSet_77797) <- substr(sampleNames(rgSet_77797), 1, 10)
pD <- pD[sampleNames(rgSet_77797),]
pD <- as(pD, "DataFrame")
pData(rgSet_77797) <- pD
rgSet_77797

grSet_77797 <- preprocessQuantile(rgSet_77797)
getBeta(grSet_77797)[1:3,1:3]
head(getIslandStatus(grSet_77797))
betaMat_77797 <- getBeta(grSet_77797)
## 485512      18
res1 <- MethylDeconv(betaMat_77797, method = "Houseman")
res2 <- MethylDeconv(betaMat_77797, method = "RPC")
res3 <- MethylDeconv(betaMat_77797, method = "CBS")
###
res1_1 <- MethylDeconv(rgSet_77797, method = "Houseman", normalized = FALSE)
#res2_1 <- MethylDeconv(rgSet_77797, method = "RPC", normalized = FALSE)
#res3_1 <- MethylDeconv(rgSet_77797, method = "RPC", normalized = FALSE)

pD_77797 <- pD.all[, c("title", "geo_accession", "b cell (%):ch1", "cd4+ t cell (%):ch1", "cd8+ t cell (%):ch1",
                             "granulocyte (%):ch1", "monocyte (%):ch1", "natural killer cell (%):ch1")]
pD_77797 <- pD_77797[sampleNames(rgSet_77797),]
facs_77797 <- pD_77797[,3:8]
colnames(facs_77797) <- c("Bcell", "CD4T", "CD8T","Gran","Mono","NK")
facs_77797_prop <- cbind(facs_77797[,"CD8T"], facs_77797[,"CD4T"], facs_77797[,"NK"],
                          facs_77797[,"Bcell"],facs_77797[,"Mono"],facs_77797[,"Gran"])
facs_77797_prop <- as.data.frame(facs_77797_prop)
rownames(facs_77797_prop) = rownames(res1)
colnames(facs_77797_prop) = colnames(res1)
for (i in 1:6){
  facs_77797_prop[,i] <- as.numeric(as.character(facs_77797_prop[,i]))
}
facs_77797_prop <- facs_77797_prop/100

rsquared = corr = rmse = rep(NA, 6)
library(hydroGOF)
for (i in 1:6){
  rsquared[i] <- summary(lm(res1[,i]~facs_77797_prop[,i]))$r.squared
  corr[i] <- cor(res1[,i], facs_77797_prop[,i], method = "spearman")
}
rmse <- mse(res1, facs_77797_prop)

library(tidyr)
library(ggplot2)
df <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"), Rsquared = rsquared,
                 SpearmanCorr = corr)
ggplot(data = df %>% gather(Variable, values, -cellType), 
       aes(x = cellType, y = values, fill = Variable)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(values,2)), position = position_dodge(0.9))



rsquared = corr = rmse = rep(NA, 6)
library(hydroGOF)
for (i in 1:6){
  rsquared[i] <- summary(lm(res2[,i]~facs_77797_prop[,i]))$r.squared
  corr[i] <- cor(res2[,i], facs_77797_prop[,i], method = "spearman")
}
rmse <- mse(res2, facs_77797_prop)

library(tidyr)
library(ggplot2)
df <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"), Rsquared = rsquared,
                 SpearmanCorr = corr)
ggplot(data = df %>% gather(Variable, values, -cellType), 
       aes(x = cellType, y = values, fill = Variable)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(values,2)), position = position_dodge(0.9))



rsquared = corr = rmse = rep(NA, 6)
library(hydroGOF)
for (i in 1:6){
  rsquared[i] <- summary(lm(res3[,i]~facs_77797_prop[,i]))$r.squared
  corr[i] <- cor(res3[,i], facs_77797_prop[,i], method = "spearman")
}
rmse <- mse(res3, facs_77797_prop)

library(tidyr)
library(ggplot2)
df <- data.frame(cellType = c("CD8T","CD4T","NK","Bcell","Mono","Neu"), Rsquared = rsquared,
                 SpearmanCorr = corr)
ggplot(data = df %>% gather(Variable, values, -cellType), 
       aes(x = cellType, y = values, fill = Variable)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(values,2)), position = position_dodge(0.9))











################################################################################################################
################################################################################################################
#######EWAS 
# GSE30870 age 40 samples 450k
library(GEOquery)
geoMat_30870<- getGEO("GSE30870")
pD.all <- pData(geoMat_30870[[1]])
# 40 37

pD <- pD.all[, c("title", "geo_accession", "source_name_ch1",  "age:ch1", "cell type:ch1", "disease status:ch1",
                 "tissue:ch1")]

getGEOSuppFiles("GSE30870")
untar("GSE30870/GSE30870_RAW.tar", exdir = "GSE30870/idat")
head(list.files("GSE30870/idat", pattern = "idat"))
## not idat files, gpl files

betaMat_30870 <- read.table("GSE30870_series_matrix.txt", skip = 63, header = T, sep = "\t", nrows = 485577)

rownames(betaMat_30870) <- betaMat_30870[,1]
betaMat_30870 <- betaMat_30870[,-1]
# dim 485577     40

res1 <- MethylDeconv(betaMat_30870, method = "Houseman")
res2 <- MethylDeconv(betaMat_30870, method = "RPC")
#res3 <- MethylDeconv(betaMat_30870, method = "CBS")

age_dat_blood <- matrix(NA, 40, 7)
samples <- rownames(res1)
rownames(age_dat_blood) <- samples
age_dat_blood[,1:6] <- res1
age_dat_blood[,7] <- as.character(pD[,3])
age_dat_blood <- as.data.frame(age_dat_blood)
age_dat_blood[1:6] <- apply(age_dat_blood[1:6], 2, as.numeric)
age_dat_blood[,7] <- as.character(age_dat_blood[,7])
colnames(age_dat_blood) <- c(colnames(res1),"age")

library(ggplot2)
library(tidyr)
df <- gather(age_dat_blood, series,value,-age)
ggplot(df) + geom_boxplot(aes(series ,value,color=age)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE30870-Houseman")

# betaMat_30870_cordBlood <- betaMat_30870[,c(2,22:40)]
# betaMat_30870_wholeBlood <- betaMat_30870[,c(1,3:21)]
# res1_cordBlood <- MethylDeconv(betaMat_30870_cordBlood, tissue = "CordBlood", method = "Houseman")
# res1[c(2,22:40),]












# GSE50660 smoking 464 samples, blood
library(GEOquery)
geoMat_50660<- getGEO("GSE50660")
pD.all <- pData(geoMat_50660[[1]])
pD <- pD.all[, c("title", "geo_accession",  "age:ch1", "gender:ch1", 
                 "smoking (0, 1 and 2, which represent never, former and current smokers):ch1",
                 "tissue:ch1")]
getGEOSuppFiles("GSE50660")
untar("GSE50660/GSE50660_RAW.tar", exdir = "GSE50660/idat")
head(list.files("GSE50660/idat", pattern = "idat"))
# no idat files ----> processed.txt
betaMat_50660 <- read.table("GSE50660/GSE50660_matrix_processed.txt", header = T, sep = "\t")
rownames(betaMat_50660) <- betaMat_50660[,1]
betaMat_50660 <- betaMat_50660[,-1]
colnames(betaMat_50660) <- rownames(pD)
res1 <- MethylDeconv(betaMat_50660, method = "Houseman")
res2 <- MethylDeconv(betaMat_50660, method = "RPC")
res3 <- MethylDeconv(betaMat_50660, method = "CBS")

gender_dat_blood <- matrix(NA, 464, 7)
samples <- rownames(res1)
rownames(gender_dat_blood) <- samples
gender_dat_blood[,1:6] <- res1
gender_dat_blood[,7] <- pD[,4]
gender_dat_blood <- as.data.frame(gender_dat_blood)
gender_dat_blood[1:6] <- apply(gender_dat_blood[1:6], 2, as.numeric)
gender_dat_blood[,7] <- as.character(gender_dat_blood[,7])
colnames(gender_dat_blood) <- c(colnames(res1),"gender")

library(ggplot2)
library(tidyr)
df <- gather(gender_dat_blood, series,value,-gender)
ggplot(df) + geom_boxplot(aes(series ,value,color=gender)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE50660-Houseman")



smoking_dat_blood <- matrix(NA, 464, 7)
samples <- rownames(res1)
rownames(smoking_dat_blood) <- samples
smoking_dat_blood[,1:6] <- res1
smoking_dat_blood[,7] <- pD[,5]
smoking_dat_blood <- as.data.frame(smoking_dat_blood)
smoking_dat_blood[1:6] <- apply(smoking_dat_blood[1:6], 2, as.numeric)
smoking_dat_blood[,7] <- as.character(smoking_dat_blood[,7])
for (i in 1:464){
  if (smoking_dat_blood[i,7] == "0"){
    smoking_dat_blood[i,7] <- "never"
  }
  if (smoking_dat_blood[i,7] == "1"){
    smoking_dat_blood[i,7] <- "former"
  }
  if (smoking_dat_blood[i,7] == "2"){
    smoking_dat_blood[i,7] <- "current"
  }
}
colnames(smoking_dat_blood) <- c(colnames(res1),"smoking")

library(ggplot2)
library(tidyr)
df <- gather(smoking_dat_blood, series,value,-smoking)
ggplot(df) + geom_boxplot(aes(series ,value,color= smoking)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE50660-Houseman")













# GSE42861 smoking ## 689 samples, have idat files!
library(GEOquery)
geoMat_42861<- getGEO("GSE42861")
pD.all <- pData(geoMat_42861[[1]])












# whole blood of 30 SLE patients and normal 25 controls
# GSE82221
library(GEOquery)
geoMat_82221<- getGEO("GSE82221")
pD.all <- pData(geoMat_82221[[1]])
# 55 samples
pD <- pD.all[, c("title", "geo_accession",  "age:ch1", "Sex:ch1", "tissue:ch1")]

betaMat_82221 <- read.table("GSE82221-GPL13534_series_matrix.txt", skip = 60, header = T, sep = "\t", nrows = 485577)

rownames(betaMat_82221) <- betaMat_82221[,1]
betaMat_82221 <- betaMat_82221[,-1]
# dim 485577     55

### sample names are not the same!!
res1 <- MethylDeconv(betaMat_82221, method = "Houseman")
res2 <- MethylDeconv(betaMat_82221, method = "RPC")
#res3 <- MethylDeconv(betaMat_82221, method = "CBS")

# patient_dat_blood <- matrix(NA, 55, 7)
# samples <- rownames(res1)
# rownames(patient_dat_blood) <- samples
# patient_dat_blood[,1:6] <- res1
# patient_dat_blood[,7] <- pD[,1]
# 
# patient_dat_blood <- as.data.frame(patient_dat_blood)
# patient_dat_blood[1:6] <- apply(patient_dat_blood[1:6], 2, as.numeric)
# patient_dat_blood[,7] <- as.character(patient_dat_blood[,7])
# 
# colnames(patient_dat_blood) <- c(colnames(res1),"group")
# 
# library(ggplot2)
# library(tidyr)
# df <- gather(patient_dat_blood, series,value,-group)
# ggplot(df) + geom_boxplot(aes(series ,value,color= group)) +
#   xlab('cell types')+
#   ylab('proportions') +
#   ggtitle("GSE82221-Houseman")
# 






















## EPIC GBM, GSE116298  47 samples
library(GEOquery)
geoMat_116298<- getGEO("GSE116298")
pD.all <- pData(geoMat_116298[[1]])
pD <- pD.all[, c("title", "geo_accession", "gender:ch1", "tissue:ch1")]

getGEOSuppFiles("GSE116298")
head(list.files("/Users/junesong/Desktop/causal inference/CellProportion/methylDeConv/GSE116298/GSE116298_RAW"))

idatFiles <- list.files("/Users/junesong/Desktop/causal inference/CellProportion/methylDeConv/GSE116298/GSE116298_RAW", 
                        pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)


library(minfi)
rgSet_116298 <- read.metharray.exp("/Users/junesong/Desktop/causal inference/CellProportion/methylDeConv/GSE116298/GSE116298_RAW",
                                   force=TRUE)
#dim: 1051539 47 

sampleNames(rgSet_116298) <- substr(sampleNames(rgSet_116298), 1, 10)
pD <- pD[sampleNames(rgSet_116298),]
pD <- as(pD, "DataFrame")
pData(rgSet_116298) <- pD
rgSet_116298

grSet_116298 <- preprocessNoob(rgSet_116298)
getBeta(grSet_116298)[1:3,1:3]
head(getIslandStatus(grSet_116298))
betaMat_116298 <- getBeta(grSet_116298)
## 865859     47
res1 <- MethylDeconv_BloodEPIC(betaMat_116298, method = "Houseman")
res2 <- MethylDeconv_BloodEPIC(betaMat_116298, method = "RPC")
res3 <- MethylDeconv_BloodEPIC(betaMat_116298, method = "CBS")

pD.all <- pData(geoMat_116298[[1]])
pD <- pD.all[, c("title", "geo_accession", "gender:ch1", "tissue:ch1")]
tissue_dat_blood <- matrix(NA, 47, 7)
samples <- rownames(res2)
rownames(tissue_dat_blood) <- samples
tissue_dat_blood[,1:6] <- res2
tissue_dat_blood[,7] <- pD[,4]
tissue_dat_blood <- as.data.frame(tissue_dat_blood)
tissue_dat_blood[1:6] <- apply(tissue_dat_blood[1:6], 2, as.numeric)
tissue_dat_blood[,7] <- as.character(tissue_dat_blood[,7])

colnames(tissue_dat_blood) <- c(colnames(res2),"tissue")

library(ggplot2)
library(tidyr)
df <- gather(tissue_dat_blood, series,value,-tissue)
ggplot(df) + geom_boxplot(aes(series ,value,color= tissue)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE116298-EPIC-RPC")


pD <- pD.all[, c("title", "geo_accession", "gender:ch1", "tissue:ch1")]
gender_dat_blood <- matrix(NA, 47, 7)
samples <- rownames(res2)
rownames(gender_dat_blood) <- samples
gender_dat_blood[,1:6] <- res2
gender_dat_blood[,7] <- pD[,3]
gender_dat_blood <- as.data.frame(gender_dat_blood)
gender_dat_blood[1:6] <- apply(gender_dat_blood[1:6], 2, as.numeric)
gender_dat_blood[,7] <- as.character(gender_dat_blood[,7])

colnames(gender_dat_blood) <- c(colnames(res2),"gender")

library(ggplot2)
library(tidyr)
df <- gather(gender_dat_blood, series,value,-gender)
ggplot(df) + geom_boxplot(aes(series ,value,color= gender)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE116298-EPIC-RPC")
