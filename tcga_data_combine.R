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