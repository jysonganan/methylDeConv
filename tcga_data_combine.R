#setwd("/Users/junesong/Desktop/causal inference/TCGA_SKCM")
#dat <- read.table("all.tsv", header = T, sep = "\t")

#1.
setwd("/Users/junesong/Desktop/causal inference/TCGA_skcm_LinkedOmics")
dat <- read.table("Human__TCGA_SKCM__JHU_USC__Methylation__Meth450__01_28_2016__BI__CpG__Firehose_Methylation_Prepocessor.cct", 
                  header = T, sep = "\t")
# 335645 CpGs, 105 samples


# a single sample downloaded
dat1 <- read.table("/Users/junesong/Desktop/causal inference/TCGA_skcm_/jhu-usc.edu_SKCM.HumanMethylation450.13.lvl-3.TCGA-Z2-A8RT-06A-11D-A373-05.gdc_hg38.txt",
                   header = T, sep = "\t")

#RTCGA
library(RTCGA)
checkTCGA('Dates')
# Download Clinical data
downloadTCGA(cancerTypes = "SKCM", destDir = "/Users/junesong/Desktop/causal inference", date = "2016-01-28")
clinicalTab <- read.table("/Users/junesong/Desktop/causal inference/gdac.broadinstitute.org_SKCM.Merge_Clinical.Level_1.2016012800.0.0/SKCM.clin.merged.txt",
                          header = F, sep = "\t", fill = T, quote = "")

dim(clinicalTab) #1912 471 (470 samples)
#barcode
clinicalTab[12,]

##
datInfo <- checkTCGA("DataSets","SKCM", date = "2016-01-28")

# Download methylation data
downloadTCGA(cancerTypes = "SKCM", dataSet = "Methylation_Preprocess", 
             destDir = "/Users/junesong/Desktop/causal inference", date = "2016-01-28")
downloadTCGA(cancerTypes = "SKCM", dataSet = "Merge_methylation__humanmethylation450", 
             destDir = "/Users/junesong/Desktop/causal inference", date = "2016-01-28")

methylTab <- read.table("/Users/junesong/Desktop/causal inference/gdac.broadinstitute.org_SKCM.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/SKCM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                        header = F, sep = "\t", fill = T, quote = "")
# 485579
path <- list.files(path = "/Users/junesong/Desktop/causal inference/gdac.broadinstitute.org_SKCM.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/MethylDataset", full.names = TRUE, recursive = TRUE)
MethylTab <- readTCGA(path, dataType = "methylation")






## KICH
downloadTCGA(cancerTypes = "KICH", dataSet = "Merge_methylation__humanmethylation450", 
             destDir = "/Users/junesong/Desktop/causal inference", date = "2016-01-28")
path_KICH <- list.files(path = "/Users/junesong/Desktop/causal inference/gdac.broadinstitute.org_KICH.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/MethylDataset", full.names = TRUE, recursive = TRUE)
MethylTab_KICH <- readTCGA(path_KICH,dataType = "methylation")
#  dim(MethylTab_KICH)
#  66 485578
BetaMatrix_KICH <- matrix(NA,485577,66)
BetaMatrix_KICH <- t(MethylTab_KICH[,-1])
colnames(BetaMatrix_KICH) <- as.character(MethylTab_KICH[,1])
rownames(BetaMatrix_KICH) <- colnames(MethylTab_KICH)[-1]

## dim: 485577 66

BetaMatrix_KICH_noNA <- na.omit(BetaMatrix_KICH)
# dim 391298     66
res1 <- MethylDeconv_normalized(BetaMatrix_KICH_noNA, method = "Houseman", tissue = "Blood", custom_probes = NULL)
res2 <- MethylDeconv_normalized(BetaMatrix_KICH_noNA, method = "RPC", tissue = "Blood", custom_probes = NULL)
res3 <- MethylDeconv_normalized(BetaMatrix_KICH_noNA, method = "CBS", tissue = "Blood", custom_probes = NULL)
library(hydroGOF)
mse(res1,res2)
mse(res1,res3)
mse(res2,res3)

res4 <- MethylDeconv_normalized(BetaMatrix_KICH_noNA, method = "Houseman", tissue = "genericEpithelial", custom_probes = NULL)
res5 <- MethylDeconv_normalized(BetaMatrix_KICH_noNA, method = "RPC", tissue = "genericEpithelial", custom_probes = NULL)
res6 <- MethylDeconv_normalized(BetaMatrix_KICH_noNA, method = "CBS", tissue = "genericEpithelial", custom_probes = NULL)
mse(res4[[2]],res5[[2]])
mse(res4[[2]],res6[[2]])
mse(res5[[2]],res6[[2]])


boxplot(res1, main = "KICH: Blood-Houseman")
boxplot(res2, main = "KICH: Blood-RPC")
boxplot(res3, main = "KICH: Blood-CIBERSORT")


boxplot(res4[[2]][,c("CD8T","CD4T","NK","B","Mono","Neutro","Eosino","Epi","Fib")], main = "KICH: Epithelial-Houseman")
boxplot(res5[[2]][,c("CD8T","CD4T","NK","B","Mono","Neutro","Eosino","Epi","Fib")], main = "KICH: Epithelial-RPC")
boxplot(res6[[2]][,c("CD8T","CD4T","NK","B","Mono","Neutro","Eosino","Epi","Fib")], main = "KICH: Epithelial-CIBERSORT")

# Download Clinical data
downloadTCGA(cancerTypes = "KICH", destDir = "/Users/junesong/Desktop/causal inference", date = "2016-01-28")
clinicalTab <- read.table("/Users/junesong/Desktop/causal inference/gdac.broadinstitute.org_KICH.Merge_Clinical.Level_1.2016012800.0.0/KICH.clin.merged.txt",
                          header = F, sep = "\t", fill = T, quote = "")
rownames(clinicalTab) <- as.character(clinicalTab[,1])
clinicalTab <- clinicalTab[,-1]
colnames(clinicalTab) <- NULL
colnames(clinicalTab) <- as.character(as.matrix(clinicalTab[12,]))
clinicalTab <- clinicalTab[c(11,244,265,276,323,325),]
colnames(clinicalTab) <- toupper(colnames(clinicalTab))

age_dat_blood <- matrix(NA, 66, 7)
samples <- substr(rownames(res1),1,12)
rownames(age_dat_blood) <- samples
age_dat_blood[,1:6] <- res1
age_dat_blood[,7] <- as.character(clinicalTab[1, match(samples,colnames(clinicalTab))])
age_dat_blood <- apply(age_dat_blood, 2, as.numeric)
colnames(age_dat_blood) <- c(colnames(res1),"age")
age_dat_blood <- as.data.frame(age_dat_blood)

age_dat_blood_epithelial <- matrix(NA, 66, 10)
samples <- substr(rownames(res4[[2]]),1,12)
rownames(age_dat_blood_epithelial) <- samples
age_dat_blood_epithelial[,1:9] <- res4[[2]]
age_dat_blood_epithelial[,10] <- as.character(clinicalTab[1, match(samples,colnames(clinicalTab))])
age_dat_blood_epithelial <- apply(age_dat_blood_epithelial, 2, as.numeric)
colnames(age_dat_blood_epithelial) <- c(colnames(res4[[2]]),"age")
age_dat_blood_epithelial <- as.data.frame(age_dat_blood_epithelial)


library(ggplot2)
library(tidyr)
# p = ggplot() + 
#   geom_point(data = age_dat_blood, aes(x = age, y = CD8T), color = "blue") +
#   geom_point(data = age_dat_blood, aes(x = age, y = CD4T), color = "red") +
#   geom_point(data = age_dat_blood, aes(x = age, y = NK), color = "black") +
#   geom_point(data = age_dat_blood, aes(x = age, y = Bcell), color = "green") +
#   geom_point(data = age_dat_blood, aes(x = age, y = Mono), color = "yellow") +
#   geom_point(data = age_dat_blood, aes(x = age, y = Gran), color = "purple") +


df <- gather(age_dat_blood, series,value,-age)
ggplot(df) + geom_point(aes(age ,value,color=series)) +
  xlab('age') +
  ylab('proportions') +
  ggtitle("KICH: Blood-Houseman")

ggplot(df) + geom_point(aes(age ,value,color=series)) +
  geom_smooth(aes(age ,value,color=series)) +
  xlab('age') +
  ylab('proportions') +
  ggtitle("KICH: Blood-Houseman")

df <- gather(age_dat_blood_epithelial, series,value,-age)
ggplot(df) + geom_point(aes(age ,value,color=series)) +
  geom_smooth(aes(age ,value,color=series)) +
  xlab('age') +
  ylab('proportions') +
  ggtitle("KICH: Epithelial-Houseman")


gender_dat_blood <- matrix(NA, 66, 7)
samples <- substr(rownames(res1),1,12)
rownames(gender_dat_blood) <- samples
gender_dat_blood[,1:6] <- res1
gender_dat_blood[,7] <- as.character(as.matrix(clinicalTab[2, match(samples,colnames(clinicalTab))]))
gender_dat_blood <- as.data.frame(gender_dat_blood)
gender_dat_blood[1:6] <- apply(gender_dat_blood[1:6], 2, as.numeric)
colnames(gender_dat_blood) <- c(colnames(res1),"gender")

df <- gather(gender_dat_blood, series,value,-gender)
ggplot(df) + geom_boxplot(aes(series ,value,color=gender)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("KICH: Blood-Houseman")

gender_dat_epithelial <- matrix(NA, 66, 10)
samples <- substr(rownames(res1),1,12)
rownames(gender_dat_epithelial) <- samples
gender_dat_epithelial[,1:9] <- res4[[2]]
gender_dat_epithelial[,10] <- as.character(as.matrix(clinicalTab[2, match(samples,colnames(clinicalTab))]))
gender_dat_epithelial <- as.data.frame(gender_dat_epithelial)
gender_dat_epithelial[1:9] <- apply(gender_dat_epithelial[1:9], 2, as.numeric)
colnames(gender_dat_epithelial) <- c(colnames(res4[[2]]),"gender")

df <- gather(gender_dat_epithelial, series,value,-gender)
ggplot(df) + geom_boxplot(aes(series ,value,color=gender)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("KICH: Epithelial-Houseman")





gender_dat_epithelial <- matrix(NA, 66, 10)
samples <- substr(rownames(res4[[2]]),1,12)
rownames(age_dat_blood_epithelial) <- samples
age_dat_blood_epithelial[,1:9] <- res4[[2]]
age_dat_blood_epithelial[,10] <- as.character(clinicalTab[1, match(samples,colnames(clinicalTab))])
age_dat_blood_epithelial <- apply(age_dat_blood_epithelial, 2, as.numeric)
colnames(age_dat_blood_epithelial) <- c(colnames(res4[[2]]),"age")
age_dat_blood_epithelial <- as.data.frame(age_dat_blood_epithelial)


##### Survival plots
clinicalTab <- read.table("/Users/junesong/Desktop/causal inference/gdac.broadinstitute.org_KICH.Merge_Clinical.Level_1.2016012800.0.0/KICH.clin.merged.txt",
                          header = F, sep = "\t", fill = T, quote = "")
rownames(clinicalTab) <- as.character(clinicalTab[,1])
clinicalTab <- clinicalTab[,-1]
colnames(clinicalTab) <- NULL
colnames(clinicalTab) <- as.character(as.matrix(clinicalTab[12,]))
survTab <- clinicalTab[c(17,19, 325),]
colnames(survTab) <- toupper(colnames(survTab))
samples <- substr(rownames(res1),1,12)
survTab <- survTab[,match(samples, colnames(survTab))]
survTab <- t(survTab)
survTab <- as.data.frame(survTab)
survTab[,1] <- as.numeric(as.character(survTab[,1]))
survTab[,2] <- as.numeric(as.character(survTab[,2]))
survTab <- cbind(survTab,res1)
survTab$status <- as.numeric(survTab$patient.vital_status == "dead")

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

survTab$days <- rep(NA,66)
survTab[which(survTab$status == 1),11] <- survTab[which(survTab$status == 1),1]
survTab[which(survTab$status == 0),11] <- survTab[which(survTab$status == 0),2]

survTab$trt <- factor((survTab$CD4T < 0.367),labels = c("high_CD4T","low_CD4T"))
## 0 cencored 1 observed.
km_trt_fit <- survfit(Surv(days, status) ~ trt, data=survTab)
summary(km_trt_fit)
autoplot(km_trt_fit)
autoplot(km_trt_fit,main = "KICH: Blood-Houseman")









library(magrittr)
infoTCGA()
head()
knitr::kable()

installTCGA("RTCGA.methylation.20160128") 
data(BRCA.methylation.20160128)

if (!require(devtools)) {
  install.packages("devtools")
  require(devtools)
}
install_github("RTCGA/RTCGA", build_vignettes = TRUE)




install.packages("remotes")
remotes::install_github("RTCGA/RTCGA.clinical")

library(RTCGA)
installTCGA("RTCGA.clinical")

## bioconductor only 2015-11-01
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RTCGA.clinical")