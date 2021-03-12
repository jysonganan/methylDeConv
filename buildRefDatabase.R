## 450k
library(FlowSorted.Blood.450k)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]

###### add epithelial
# library(GEOquery)
# library(minfi)
# rgSet <- read.metharray.exp("GSE40699/idat")
# ref_betamatrix_1 <- getBeta(preprocessNoob(rgSet, dyeMethod = "single"))
# 
# geoMat <- getGEO("GSE40699")
# pD.all <- pData(geoMat[[1]])
# 
# ref_betamatrix_EpiFib <- cbind(ref_betamatrix_1[,match(rownames(pD.all)[c(2,12,27,21,28,35,44,46,50,51,56)],substr(colnames(rgSet),1,9))],
#                                ref_betamatrix_1[,match(rownames(pD.all)[c(6,8,10,11,14,16,60)],substr(colnames(rgSet),1,9))])
# ref_phenotype_EpiFib <- c(rep("Epithelial", 11), rep("Fibroblast", 7))
# save("ref_betamatrix_EpiFib", "ref_phenotype_EpiFib", file = "ref_GSE40699_EpiFib.RData")


load("ref_GSE40699_EpiFib.RData")
ref_betamatrix <- cbind(ref_betamatrix, ref_betamatrix_EpiFib)
ref_phenotype <- c(ref_phenotype, ref_phenotype_EpiFib)




# ##### 1. oneVsAllttest
# set.seed(2)
# source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/refCompTableProbeSelection.R")
# probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both")
# 
# ##### 2. glmnetpreselect
# library(dplyr)
# library(caret)
# library(doParallel)
# library(matrixStats)
# library(plyr)
# library(recipes)
# library(adabag)
# 
# set.seed(2)
# nCores = 4
# 
# source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/refCompTableProbeSelection.R")
# probes <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype,probeSelect = "both", MaxDMRs = 300)
# ProbePreselect_multiclassGlmnet <- ref_probe_selection_multiclassGlmnet_cv(ref_betamatrix[probes,], ref_phenotype)

save("probes_oneVsAllttest", "ProbePreselect_multiclassGlmnet", file = "Flow450k_EpiFib_Probes.RData")
########################################
## 800 and 774 probes (310 intersection)




#### TCGA SKCM analysis

## or

#https://xenabrowser.net/datapages
x <- read.table("TCGA-SKCM.methylation450.tsv", header = T, sep = "\t")
rownames(x) <- x[,1]
x <- x[,-1]
## 485577 475
BetaMatrix_noNA <- na.omit(x)





source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)


disease <- "SKCM"
methyl_path <- paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_methyl/", disease, "_methyl.RData")
load(methyl_path)
dim(BetaMatrix) ## 485577 475
BetaMatrix_noNA <- na.omit(BetaMatrix) ## 373814    475

probes_select <- probes_oneVsAllttest
library(EpiDISH)
source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/projectCellType.R")
Houseman_res <- projectCellType(BetaMatrix_noNA[intersect(rownames(BetaMatrix_noNA),probes_select),],
                                as.matrix(compTable[intersect(rownames(BetaMatrix_noNA),probes_select),3:10]))
probes_select <- ProbePreselect_multiclassGlmnet[[1]][-1]
Houseman_res_glmnetpreselect <- projectCellType(BetaMatrix_noNA[intersect(rownames(BetaMatrix_noNA),probes_select),],
                                                as.matrix(compTable[intersect(rownames(BetaMatrix_noNA),probes_select),3:10]))
library(ggplot2)
library(tidyr)
pdf(file = "Plot1.pdf",  width = 10, height = 8) 
boxplot(Houseman_res, main = "Flow450k_EpiFib, Houseman, oneVsAllttest")
dev.off()


clinicalTab <- read.table("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_clinical/gdac.broadinstitute.org_SKCM.Merge_Clinical.Level_1.2016012800.0.0/SKCM.clin.merged.txt",
                          header = F, sep = "\t", fill = T, quote = "")
rownames(clinicalTab) <- as.character(clinicalTab[,1])
clinicalTab <- clinicalTab[,-1]
colnames(clinicalTab) <- NULL
colnames(clinicalTab) <- as.character(as.matrix(clinicalTab[12,]))
colnames(clinicalTab) <- toupper(colnames(clinicalTab))

gender_dat_epithelial <- matrix(NA, 475, 9)
samples <- substr(rownames(Houseman_res),1,12)
rownames(gender_dat_epithelial) <- samples
gender_dat_epithelial[,1:8] <- Houseman_res
gender_dat_epithelial[,9] <- as.character(as.matrix(clinicalTab[931, match(samples,colnames(clinicalTab))]))
gender_dat_epithelial <- as.data.frame(gender_dat_epithelial)
gender_dat_epithelial[1:8] <- apply(gender_dat_epithelial[1:8], 2, as.numeric)
colnames(gender_dat_epithelial) <- c(colnames(Houseman_res),"tumor_tissue_site")

pdf(file = "Plot1.pdf",  width = 16, height = 8) 
df <- gather(gender_dat_epithelial, series,value,-tumor_tissue_site)
ggplot(df) + geom_boxplot(aes(series ,value,color=tumor_tissue_site)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("SKCM, Houseman, oneVsAllttest")+theme(legend.position="bottom")
dev.off()



gender_dat_epithelial <- matrix(NA, 475, 9)
samples <- substr(rownames(Houseman_res),1,12)
rownames(gender_dat_epithelial) <- samples
gender_dat_epithelial[,1:8] <- Houseman_res
gender_dat_epithelial[,9] <- as.character(as.matrix(clinicalTab[1780, match(samples,colnames(clinicalTab))]))
gender_dat_epithelial <- as.data.frame(gender_dat_epithelial)
gender_dat_epithelial[1:8] <- apply(gender_dat_epithelial[1:8], 2, as.numeric)
colnames(gender_dat_epithelial) <- c(colnames(Houseman_res),"malignancy_type")

pdf(file = "Plot1.pdf",  width = 12, height = 8) 
df <- gather(gender_dat_epithelial, series,value,-malignancy_type)
ggplot(df) + geom_boxplot(aes(series ,value,color=malignancy_type)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("SKCM, Houseman, oneVsAllttest")+theme(legend.position="bottom")
dev.off()



table(as.character(as.matrix(clinicalTab[1780,])))
#prior malignancy synchronous malignancy 
#20                      6 



table(as.character(as.matrix(clinicalTab[931,])))  ## tissue source site
#distant metastasis 
#68 
#primary tumor 
#103 
#regional cutaneous or subcutaneous tissue (includes satellite and in-transit metastasis) 
#74 
#regional lymph node 
#222



gender_dat_epithelial <- matrix(NA, 475, 9)
samples <- substr(rownames(Houseman_res),1,12)
rownames(gender_dat_epithelial) <- samples
gender_dat_epithelial[,1:8] <- Houseman_res
gender_dat_epithelial[,9] <- as.character(as.matrix(clinicalTab[906, match(samples,colnames(clinicalTab))]))
gender_dat_epithelial <- as.data.frame(gender_dat_epithelial)
gender_dat_epithelial[1:8] <- apply(gender_dat_epithelial[1:8], 2, as.numeric)
colnames(gender_dat_epithelial) <- c(colnames(Houseman_res),"patient.sites_of_primary_melanomas.site.tumor_tissue_site")

pdf(file = "Plot1.pdf",  width = 12, height = 8) 
df <- gather(gender_dat_epithelial, series,value,-patient.sites_of_primary_melanomas.site.tumor_tissue_site)
ggplot(df) + geom_boxplot(aes(series ,value,color=patient.sites_of_primary_melanomas.site.tumor_tissue_site)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("SKCM, Houseman, oneVsAllttest")+theme(legend.position="bottom")
dev.off()


table(as.character(as.matrix(clinicalTab[917,]))) #patient.stage_event.pathologic_stage

#i/ii nos    stage 0    stage i   stage ia   stage ib   stage ii  stage iia 
#14          7         30         18         29         30         18 
#stage iib  stage iic  stage iii stage iiia stage iiib stage iiic   stage iv 
#28         64         41         16         46         68         23 



table(as.character(as.matrix(clinicalTab[906,])))

#extremities  head and neck other  specify          trunk 
#196             37             13            172 



table(as.character(as.matrix(clinicalTab[772,]))) #radiation therapy
#no yes 
#438  23 



gender_dat_epithelial <- matrix(NA, 475, 9)
samples <- substr(rownames(Houseman_res),1,12)
rownames(gender_dat_epithelial) <- samples
gender_dat_epithelial[,1:8] <- Houseman_res
gender_dat_epithelial[,9] <- as.character(as.matrix(clinicalTab[767, match(samples,colnames(clinicalTab))]))
gender_dat_epithelial <- as.data.frame(gender_dat_epithelial)
gender_dat_epithelial[1:8] <- apply(gender_dat_epithelial[1:8], 2, as.numeric)
colnames(gender_dat_epithelial) <- c(colnames(Houseman_res),"patient.person_neoplasm_cancer_status")

pdf(file = "Plot1.pdf",  width = 12, height = 8) 
df <- gather(gender_dat_epithelial, series,value,-patient.person_neoplasm_cancer_status)
ggplot(df) + geom_boxplot(aes(series ,value,color=patient.person_neoplasm_cancer_status)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("SKCM, Houseman, oneVsAllttest")+theme(legend.position="bottom")
dev.off()



table(as.character(as.matrix(clinicalTab[771,])))   ### race
#asian black or african american                     white 
#12                         1                       447 


table(as.character(as.matrix(clinicalTab[767,]))) ### patient.person_neoplasm_cancer_status
#tumor free with tumor 
#222        220 

table(as.character(as.matrix(clinicalTab[756,]))) ## patient.new_tumor_events.new_tumor_event.new_neoplasm_event_type

#distant metastasis locoregional recurrence    new primary melanoma 
#91                      46                      19 
#regional lymph node 
#44 


gender_dat_epithelial <- matrix(NA, 475, 9)
samples <- substr(rownames(Houseman_res),1,12)
rownames(gender_dat_epithelial) <- samples
gender_dat_epithelial[,1:8] <- Houseman_res
gender_dat_epithelial[,9] <- as.character(as.matrix(clinicalTab[756, match(samples,colnames(clinicalTab))]))
gender_dat_epithelial <- as.data.frame(gender_dat_epithelial)
gender_dat_epithelial[1:8] <- apply(gender_dat_epithelial[1:8], 2, as.numeric)
colnames(gender_dat_epithelial) <- c(colnames(Houseman_res),"new_tumor_event.new_neoplasm_event_type")

pdf(file = "Plot1.pdf",  width = 12, height = 8) 
df <- gather(gender_dat_epithelial, series,value,-new_tumor_event.new_neoplasm_event_type)
ggplot(df) + geom_boxplot(aes(series ,value,color=new_tumor_event.new_neoplasm_event_type)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("SKCM, Houseman, oneVsAllttest")+theme(legend.position="bottom")
dev.off()




table(as.character(as.matrix(clinicalTab[436,]))) #patient.melanoma_ulceration_indicator

#no yes 
#146 167 


gender_dat_epithelial <- matrix(NA, 475, 9)
samples <- substr(rownames(Houseman_res),1,12)
rownames(gender_dat_epithelial) <- samples
gender_dat_epithelial[,1:8] <- Houseman_res
gender_dat_epithelial[,9] <- as.character(as.matrix(clinicalTab[434, match(samples,colnames(clinicalTab))]))
gender_dat_epithelial <- as.data.frame(gender_dat_epithelial)
gender_dat_epithelial[1:8] <- apply(gender_dat_epithelial[1:8], 2, as.numeric)
colnames(gender_dat_epithelial) <- c(colnames(Houseman_res),"melanoma_clark_level_value")

pdf(file = "Plot1.pdf",  width = 12, height = 8) 
df <- gather(gender_dat_epithelial, series,value,-melanoma_clark_level_value)
ggplot(df) + geom_boxplot(aes(series ,value,color=melanoma_clark_level_value)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("SKCM, Houseman, oneVsAllttest")+theme(legend.position="bottom")
dev.off()




table(as.character(as.matrix(clinicalTab[434,]))) #patient.melanoma_clark_level_value

#i  ii iii  iv   v 
#6  18  77 168  52 


table(as.character(as.matrix(clinicalTab[426,]))) #patient.history_of_neoadjuvant_treatment

#no yes 
#445  25 


gender_dat_epithelial <- matrix(NA, 475, 9)
samples <- substr(rownames(Houseman_res),1,12)
rownames(gender_dat_epithelial) <- samples
gender_dat_epithelial[,1:8] <- Houseman_res
gender_dat_epithelial[,9] <- as.character(as.matrix(clinicalTab[426, match(samples,colnames(clinicalTab))]))
gender_dat_epithelial <- as.data.frame(gender_dat_epithelial)
gender_dat_epithelial[1:8] <- apply(gender_dat_epithelial[1:8], 2, as.numeric)
colnames(gender_dat_epithelial) <- c(colnames(Houseman_res),"history_of_neoadjuvant_treatment")

pdf(file = "Plot1.pdf",  width = 12, height = 8) 
df <- gather(gender_dat_epithelial, series,value,-history_of_neoadjuvant_treatment)
ggplot(df) + geom_boxplot(aes(series ,value,color=history_of_neoadjuvant_treatment)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("SKCM, Houseman, oneVsAllttest")+theme(legend.position="bottom")
dev.off()


table(as.character(as.matrix(clinicalTab[418,]))) #patient.follow_ups.follow_up.person_neoplasm_cancer_status

#tumor free with tumor 
#182        105 


gender_dat_epithelial <- matrix(NA, 475, 9)
samples <- substr(rownames(Houseman_res),1,12)
rownames(gender_dat_epithelial) <- samples
gender_dat_epithelial[,1:8] <- Houseman_res
gender_dat_epithelial[,9] <- as.character(as.matrix(clinicalTab[424, match(samples,colnames(clinicalTab))]))
gender_dat_epithelial <- as.data.frame(gender_dat_epithelial)
gender_dat_epithelial[1:8] <- apply(gender_dat_epithelial[1:8], 2, as.numeric)
colnames(gender_dat_epithelial) <- c(colnames(Houseman_res),"gender")

pdf(file = "Plot1.pdf",  width = 12, height = 8) 
df <- gather(gender_dat_epithelial, series,value,-gender)
ggplot(df) + geom_boxplot(aes(series ,value,color=gender)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("SKCM, Houseman, oneVsAllttest")+theme(legend.position="bottom")
dev.off()



table(as.character(as.matrix(clinicalTab[424,])))

#female   male 
#180    290 

##### Survival plots
survTab <- clinicalTab[c(17,19, 932),]
samples <- substr(rownames(Houseman_res),1,12)
## patient.days_to_death; patient.days_to_last_followup; patient.follow_ups.follow_up.vital_status 
survTab <- t(survTab)
survTab <- as.data.frame(survTab)
survTab[,1] <- as.numeric(as.character(survTab[,1]))
survTab[,2] <- as.numeric(as.character(survTab[,2]))
survTab <- cbind(survTab,Houseman_res[match(rownames(survTab),samples),])
survTab$status <- as.numeric(survTab$patient.vital_status == "dead")




library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(survminer)

survTab$days <- rep(NA,dim(survTab)[1])
survTab[which(survTab$status == 1),13] <- survTab[which(survTab$status == 1),1]
survTab[which(survTab$status == 0),13] <- survTab[which(survTab$status == 0),2]


survTab$trt <- factor((survTab$Fibroblast < 0.347),labels = c("high_Fibroblast","low_Fibroblast"))
## 0 cencored 1 observed.
km_trt_fit <- survfit(Surv(days, status) ~ trt, data=survTab)
# summary(km_trt_fit)
# autoplot(km_trt_fit)
# autoplot(km_trt_fit,main = "KICH: Blood-Houseman")
## or ggplot

ggsurvplot(km_trt_fit, data = survTab, pval = TRUE, conf.int = TRUE,title = "SKCM, Houseman, high 50% vs low 50%")



survTab_sub <- survTab[survTab$CD8T < 0.114|survTab$CD8T > 0.261,] ## high 25% low 25%
survTab_sub$trt <- factor((survTab_sub$CD8T < 0.114),labels = c("high_CD8T","low_CD8T"))
km_trt_fit <- survfit(Surv(days, status) ~ trt, data=survTab_sub)
ggsurvplot(km_trt_fit, data = survTab_sub, pval = TRUE, conf.int = TRUE,title = "SKCM, Houseman, high 25% vs low 25%")


## 282 patient.follow_ups.follow_up.days_to_death
## 17 patient.days_to_death
## 11 patient.age_at_initial_pathologic_diagnosis



clinicalTab <- clinicalTab[c(11,244,265,276,323,325),]


  
  


######### add epithelial 
#save("betaMat_122126","phenotype_122126", file = "ref_122126_450kEpithelial.RData")

# epithelial cell lines HRPTEC, HREC and HRCEC
library(GEOquery)
library(minfi)
geoMat <- getGEO("GSE31848")
pD.all <- pData(geoMat[[1]])
rownames(pD.all)[which(pD.all[,"cell line:ch1"] %in% c("HRPTEpiC257","HRCEpiC255","HREpiC256"))]
#setwd("/sonas-hs/krasnitz/hpc/data/pfproj/")
# getGEOSuppFiles("GSE31848")
# untar("GSE31848/GSE31848_RAW.tar", exdir = "GSE31848/idat")
# head(list.files("GSE31848/idat", pattern = "idat"))

##  GSE31848_raw.txt but may not preprocessNoob... ranges: 0~5000, they are not beta values

betamat <- read.table("/sonas-hs/krasnitz/hpc/data/pfproj/GSE31848/GSE31848_raw.txt")









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



