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