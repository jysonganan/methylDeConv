# collapse CpGs into genes

## correct version TCGA KICH

manifest <- read.csv("HumanMethylation450_15017482_v1-2.csv", header = T, skip = 7) #486428     33
# column names: UCSC_RefGene_Name, UCSC_RefGene_Group

# For those CpGs with multiple mapping, you can use the first one (e.g use TSS200 for the sample gene below)

# UCSC_RefGene_Name: Gene1;Gene2;Gene1

# UCSC_RefGene_Group: TSS200;Body;1stExon


library(minfi)
##3 Test data of seven cell types (Gran, Mono, Bcell, CD4T, CD8T, NK and nRBCs) with matched FACS counts.
rgSet_127824 <- read.metharray.exp("GSE127824/idat")
class(rgSet_127824)
annot <- getAnnotation(rgSet_127824)
# DataFrame with 485512 rows and 33 columns
annot_df <- as.data.frame(annot)
head(annot_df[,"UCSC_RefGene_Name"],100)

grSet_127824 <- preprocessNoob(rgSet_127824)
getBeta(grSet_127824)[1:3,1:3]
head(getIslandStatus(grSet_127824))
betaMat_127824 <- getBeta(grSet_127824)

var_per_CpG <- apply(betaMat_127824,1,var)
var_per_CpG <- as.data.frame(var_per_CpG)
rownames(var_per_CpG) <- rownames(betaMat_127824)


annot_df_oneGene <- annot_df
annot_df_oneGene[,"UCSC_RefGene_Name"] <- gsub(";.*$","",annot_df[,"UCSC_RefGene_Name"])
annot_df_oneGene <- annot_df_oneGene[,c("chr", "Type", "Relation_to_Island", "UCSC_RefGene_Name")]
head(annot_df_oneGene)
## unique UCSC_RefGene_Name 20622
## For each gene, select the probe with highest variance
gene_list <- unique(annot_df_oneGene[,"UCSC_RefGene_Name"])
CpG_maxVar_list <- rep(NA, length(gene_list))

for (i in 1:length(gene_list)){
  gene <- gene_list[i]
  CpGs <- rownames(annot_df_oneGene[annot_df_oneGene[,"UCSC_RefGene_Name"] == gene,])
  var_per_CpG_sub <- var_per_CpG[CpGs,]
  CpG_maxVar_list[i] <- CpGs[which.max(var_per_CpG_sub)]
}

save("CpG_maxVar_list", file = "CpG_maxVar_list_GSE127824.RData")
genelevel_betaMat_127824 <- betaMat_127824[CpG_maxVar_list,]
rownames(genelevel_betaMat_127824) <- gene_list
genelevel_betaMat_127824 <- genelevel_betaMat_127824[-3,]
# 20621 genes


#genelevel_betaMat_127824_rank <- apply(genelevel_betaMat_127824,2,rank)
# larger, rank larger -- expression smaller
# reverse: smaller rank, smaller expression
genelevel_betaMat_127824_rank <- apply(-genelevel_betaMat_127824,2,rank)


#averaging all probes within 200KB of transcription start site (TSS) to represent a gene-level methylation. 
# Obviously 200KB is arbitrary and also the assumption of something close to TSS is very important 
#than other gene-regions made me to look for alternate ways.

# gene level methylation integrating expression and methylation data
# e.g. compute average of all CpGs within 1500 BPs of TSS.

# https://nsaunders.wordpress.com/2013/02/12/genes-x-samples-please-explain/
# several publications which tried to integrate measurements of methylation and gene expression
# 1. selecting one probe per gene using particular criteria (eg.highest variance)
# 2. complex clustering procedures based on chromosome coordinates


# TCGA does it to assign a methylation score for each is to 
#correlate all probe values in the proximity of a gene with the gene expression, 
# and pick the one that best negative correlation. 
#This is of course having some phenotype of interest in mind, 
#i.e. it works only if you want to see methylation probes that can help explain the gene expression levels across many patients




load("xCell.data.rda")
# 10808 genes and 489 signatures
length(intersect(xCell.data$genes,rownames(genelevel_betaMat_127824_rank)))
## 10694 shared
source("xCell_custom.R")
xCellScores<- xCellAnalysis(genelevel_betaMat_127824_rank, rnaseq = FALSE)
## 67 24
save(xCellScores, file = "xCellScores_127824.RData")


library(GEOquery)
geoMat_127824 <- getGEO("GSE127824")
pD.all <- pData(geoMat_127824[[1]])
pD_127824 <- pD.all[, c("title", "geo_accession", "b cells:ch1", "cd4t cells:ch1", "cd8t cells:ch1", "granulocytes:ch1", 
                        "monocytes:ch1", "nk cells:ch1", "nrbcs:ch1", "Sex:ch1", "subject status:ch1", "tissue:ch1")]
sampleNames(rgSet_127824) <- substr(sampleNames(rgSet_127824), 1, 10)
pD_127824 <- pD_127824[sampleNames(rgSet_127824),]
facs_127824 <- pD_127824[,3:9]
for (i in 1:7){
  facs_127824[,i] <- as.numeric(facs_127824[,i])
}
facs_127824_prop <- apply(facs_127824,1,function(x){return(as.numeric(x)/sum(as.numeric(x)))})
facs_127824_prop <- t(facs_127824_prop)
# 24 7

source("pipeline_cellProp.R")
res1 <- MethylDeconv(betaMat_127824, method = "Houseman", normalized = TRUE, tissue = "CordBlood")
res2 <- MethylDeconv(betaMat_127824, method = "RPC", normalized = TRUE, tissue = "CordBlood")
res3 <- MethylDeconv(betaMat_127824, method = "CBS", normalized = TRUE, tissue = "CordBlood")
# 24 7

colnames(facs_127824_prop) <- colnames(res1)

xCellScores <- t(xCellScores)
rownames(xCellScores) <- substr(rownames(xCellScores), 1, 10)

# plots
boxplot(xCellScores, las = 2,cex.axis = 0.5)

# df_plot <- cbind(cbind(res1[,1], xCellScores[,"B-cells"]),"Bcell")
# df_plot <- rbind(df_plot, cbind(cbind(res1[,2],xCellScores[,"CD4+ T-cells"]), "CD4T"))
# df_plot <- rbind(df_plot, cbind(cbind(res1[,3],xCellScores[,"CD8+ T-cells"]), "CD8T"))
# df_plot <- rbind(df_plot, cbind(cbind(res1[,5],xCellScores[,"Monocytes"]), "Mono"))
# #df_plot <- rbind(df_plot, cbind(cbind(res1[,6],xCellScores[,"NK cells"]), "NK"))
# df_plot <- as.data.frame(df_plot)
# df_plot[,1] <- as.numeric(as.character(df_plot[,1]))
# df_plot[,2] <- as.numeric(as.character(df_plot[,2]))
# df_plot[,3] <- as.character(df_plot[,3])
# colnames(df_plot) <- c("EstProp_pipeline", "xCellScore", "type")
# library(ggplot2)
# ggplot(df_plot) + geom_point(aes(EstProp_pipeline, xCellScore, color=type))+
#   geom_smooth(aes(EstProp_pipeline, xCellScore, color=type)) +
#   xlab('EstProp_pipeline') +
#   ylab('xCellScore') 

## check correlation within each sample
# df_xCellScore <- xCellScores[,c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "Monocytes")]
# df_xCellScore <- apply(df_xCellScore,1,function(x){return(x/sum(x))})
# df_xCellScore <- t(df_xCellScore)
# df_res1 <- res1[,c(1,2,3,5)]
# df_res1 <- apply(df_res1,1,function(x){return(x/sum(x))})
# df_res1 <- t(df_res1)
# df_res2 <- res2[,c(1,2,3,5)]
# df_res2 <- apply(df_res2,1,function(x){return(x/sum(x))})
# df_res2 <- t(df_res2)
# df_res3 <- res3[,c(1,2,3,5)]
# df_res3 <- apply(df_res3,1,function(x){return(x/sum(x))})
# df_res3 <- t(df_res3)
# corr <- rep(NA, 24)
# for (i in 1:24){
#   corr[i] <- cor(df_res1[i,], df_xCellScore[i,], method = "spearman")
# }
# 
# corr2 <- rep(NA, 24)
# for (i in 1:24){
#   corr2[i] <- cor(df_res2[i,], df_xCellScore[i,], method = "spearman")
# }
# 
# corr3 <- rep(NA, 24)
# for (i in 1:24){
#   corr3[i] <- cor(df_res3[i,], df_xCellScore[i,], method = "spearman")
# }
# 
# df_facs <- facs_127824_prop[,c(1,2,3,5)]
# df_facs <- apply(df_facs,1,function(x){return(x/sum(x))})
# df_facs <- t(df_facs)
# corr4 <- rep(NA, 24)
# for (i in 1:24){
#   corr4[i] <- cor(df_facs[i,], df_xCellScore[i,], method = "spearman")
# }
# 
# df <- data.frame("Houseman" = corr, "RPC" = corr2, "CBS" = corr3, "FACS" = corr4)
# stripchart(df, frame = TRUE, vertical = TRUE,
#            method = "jitter", pch = c(21, 18, 16, 14),
#            col = c("#999999", "#E69F00", "#56B4E9", "009E73"),
#            main = "Correlation with xCellScores within each sample", ylab = "Spearman correlation", ylim = c(0,1))




####################
# 67919 breast epithelial cells
#####################

betaMat_67919 <- read.table("GSE67919_series_matrix.txt", skip = 70, header = T, sep = "\t", nrows = 485577)
rownames(betaMat_67919) <- betaMat_67919[,1]
betaMat_67919 <- betaMat_67919[,-1]
# 485577 96

manifest <- read.csv("HumanMethylation450_15017482_v1-2.csv", header = T, skip = 7)
# 486428     33
annot <- manifest[match(rownames(betaMat_67919),manifest[,1]),]
# 485577     96
annot_df <- as.data.frame(annot)
head(annot_df[,"UCSC_RefGene_Name"],100)

var_per_CpG <- apply(betaMat_67919,1,var)
var_per_CpG <- as.data.frame(var_per_CpG)
rownames(var_per_CpG) <- rownames(betaMat_67919)


annot_df_oneGene <- annot_df
rownames(annot_df_oneGene) <- annot_df[,1]
annot_df_oneGene[,"UCSC_RefGene_Name"] <- gsub(";.*$","",annot_df[,"UCSC_RefGene_Name"])
annot_df_oneGene <- annot_df_oneGene[,c("CHR", "Infinium_Design_Type", "Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Name")]
head(annot_df_oneGene)

length(unique(annot_df_oneGene[,4]))
## unique UCSC_RefGene_Name 20622
## For each gene, select the probe with highest variance
gene_list <- unique(annot_df_oneGene[,"UCSC_RefGene_Name"])
CpG_maxVar_list <- rep(NA, length(gene_list))






for (i in 1:length(gene_list)){
  gene <- gene_list[i]
  CpGs <- rownames(annot_df_oneGene[annot_df_oneGene[,"UCSC_RefGene_Name"] == gene,])
  var_per_CpG_sub <- var_per_CpG[CpGs,]
  if(!is.na(var_per_CpG_sub)){
    CpG_maxVar_list[i] <- CpGs[which.max(var_per_CpG_sub)]
  }else{
    CpG_maxVar_list[i] <- NA
  }
}





save("CpG_maxVar_list", file = "CpG_maxVar_list_GSE67919.RData")

genelevel_betaMat_67919 <- betaMat_67919[CpG_maxVar_list,]
rownames(genelevel_betaMat_67919) <- gene_list
genelevel_betaMat_67919 <- genelevel_betaMat_67919[-4,]
genelevel_betaMat_67919 <- genelevel_betaMat_67919[complete.cases(genelevel_betaMat_67919),]
# 20621 genes
# remove missing
# 13495 96
genelevel_betaMat_67919_rank <- apply(-genelevel_betaMat_67919,2,rank)

load("xCell.data.rda")
# 10808 genes and 489 signatures
length(intersect(xCell.data$genes,rownames(genelevel_betaMat_67919_rank)))
## 7194 shared
source("xCell_custom.R")
xCellScores<- xCellAnalysis(genelevel_betaMat_67919_rank, rnaseq = FALSE)
## 67 24
save(xCellScores, file = "xCellScores_67919.RData")

xCellScores <- t(xCellScores)
rownames(xCellScores) <- substr(rownames(xCellScores), 1, 10)
# plots
boxplot(xCellScores, las = 2,cex.axis = 0.5)

xCellScores_breast <- xCellScores
load("xCellScores_127824.RData")
xCellScores <- t(xCellScores)
rownames(xCellScores) <- substr(rownames(xCellScores), 1, 10)
xCellScores_cordblood <- xCellScores


## plot
df_plot <- rbind(xCellScores_breast, xCellScores_cordblood)
df_plot <- cbind(df_plot, c(rep("Breast", 96), rep("CordBlood",24)))
df_plot <- as.data.frame(df_plot)
df_plot[,1:67] <- apply(df_plot[,1:67],2,function(x){return(as.numeric(as.character(x)))})
df_plot[,68] <- as.character(df_plot[,68])
colnames(df_plot)[68] <- "tissue"

df_plot_2 <- matrix(NA, 7680, 3)
df_plot_2[,1] <- as.vector(as.matrix(df_plot[,1:64]))
df_plot_2[,2] <- rep(df_plot[,68],64)
df_plot_2[,3] <- rep(colnames(df_plot)[1:64],each = 120)
colnames(df_plot_2) <- c("scores", "tissue","CellType")
df_plot_2 <- as.data.frame(df_plot_2)
df_plot_2[,1] <- as.numeric(as.character(df_plot_2[,1]))



library(ggplot2)
ggplot(df_plot_2) + geom_boxplot(aes(CellType ,scores,color= tissue)) +
  xlab('cell types')+
  ylab('scores') + theme(axis.text.x = element_text(angle = 90, hjust = 1))



df_plot_2 <- matrix(NA, 360, 3)
df_plot_2[,1] <- as.vector(as.matrix(df_plot[,65:67]))
df_plot_2[,2] <- rep(df_plot[,68],3)
df_plot_2[,3] <- rep(colnames(df_plot)[65:67],each = 120)
colnames(df_plot_2) <- c("scores", "tissue","types")
df_plot_2 <- as.data.frame(df_plot_2)
df_plot_2[,1] <- as.numeric(as.character(df_plot_2[,1]))


library(ggplot2)
ggplot(df_plot_2) + geom_boxplot(aes(types ,scores,color= tissue)) +
  xlab('types')+
  ylab('scores') + theme(axis.text.x = element_text(angle = 90, hjust = 1))












###################################
# TCGA KICH
###################################
library(RTCGA)
path_KICH <- list.files(path = "/Users/junesong/Desktop/causal inference/gdac.broadinstitute.org_KICH.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/MethylDataset", full.names = TRUE, recursive = TRUE)
MethylTab_KICH <- readTCGA(path_KICH,dataType = "methylation")
BetaMatrix_KICH <- matrix(NA,485577,66)
BetaMatrix_KICH <- t(MethylTab_KICH[,-1])
colnames(BetaMatrix_KICH) <- as.character(MethylTab_KICH[,1])
rownames(BetaMatrix_KICH) <- colnames(MethylTab_KICH)[-1]
## dim: 485577 66

manifest <- read.csv("HumanMethylation450_15017482_v1-2.csv", header = T, skip = 7)
# 486428     33
annot <- manifest[match(rownames(BetaMatrix_KICH),manifest[,1]),]
# 485577     33
annot_df <- as.data.frame(annot)
head(annot_df[,"UCSC_RefGene_Name"],100)

annot_df_oneGene <- annot_df
rownames(annot_df_oneGene) <- annot_df[,1]
annot_df_oneGene[,"UCSC_RefGene_Name"] <- gsub(";.*$","",annot_df[,"UCSC_RefGene_Name"])
annot_df_oneGene <- annot_df_oneGene[,c("CHR", "Infinium_Design_Type", "Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Name")]
head(annot_df_oneGene)
length(unique(annot_df_oneGene[,4]))
## unique UCSC_RefGene_Name 20622

annot_df_oneGene <- annot_df_oneGene[!(annot_df_oneGene[,"UCSC_RefGene_Name"]== ""),]
#365860      4
BetaMatrix_KICH <- BetaMatrix_KICH[match(rownames(annot_df_oneGene),rownames(BetaMatrix_KICH)),]
#365860      66

var_per_CpG <- apply(BetaMatrix_KICH,1,var)
var_per_CpG <- as.data.frame(var_per_CpG)
rownames(var_per_CpG) <- rownames(BetaMatrix_KICH)






## For each gene, select the probe with highest variance
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



save("CpG_maxVar_list", file = "CpG_maxVar_list_TCGAKICH.RData")

miss_id <- which(is.na(CpG_maxVar_list) == TRUE)
CpG_maxVar_list <- CpG_maxVar_list[-miss_id]
gene_list <- gene_list[-miss_id]

genelevel_BetaMatrix_KICH <- BetaMatrix_KICH[CpG_maxVar_list,]
rownames(genelevel_BetaMatrix_KICH) <- gene_list
#20473 66
genelevel_BetaMatrix_KICH <- genelevel_BetaMatrix_KICH[complete.cases(genelevel_BetaMatrix_KICH),]
# remove missing
# 20347 66
genelevel_BetaMatrix_KICH_rank <- apply(-genelevel_BetaMatrix_KICH,2,rank)

load("xCell.data.rda")
# 10808 genes and 489 signatures
length(intersect(xCell.data$genes,rownames(genelevel_BetaMatrix_KICH_rank)))
## 10667 shared
source("xCell_custom.R")
xCellScores<- xCellAnalysis(genelevel_BetaMatrix_KICH_rank, rnaseq = FALSE)
## 67 66
save(xCellScores, file = "xCellScores_KICH_methyl.RData")

xCellScores <- t(xCellScores)
rownames(xCellScores) <- substr(rownames(xCellScores), 1, 10)

# plots
boxplot(xCellScores, las = 2,cex.axis = 0.5)





### alternatively, take the average
library(RTCGA)
path_KICH <- list.files(path = "/Users/junesong/Desktop/causal inference/gdac.broadinstitute.org_KICH.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/MethylDataset", full.names = TRUE, recursive = TRUE)
MethylTab_KICH <- readTCGA(path_KICH,dataType = "methylation")
BetaMatrix_KICH <- matrix(NA,485577,66)
BetaMatrix_KICH <- t(MethylTab_KICH[,-1])
colnames(BetaMatrix_KICH) <- as.character(MethylTab_KICH[,1])
rownames(BetaMatrix_KICH) <- colnames(MethylTab_KICH)[-1]
## dim: 485577 66

manifest <- read.csv("HumanMethylation450_15017482_v1-2.csv", header = T, skip = 7)
# 486428     33
annot <- manifest[match(rownames(BetaMatrix_KICH),manifest[,1]),]
# 485577     33
annot_df <- as.data.frame(annot)
head(annot_df[,"UCSC_RefGene_Name"],100)

annot_df_oneGene <- annot_df
rownames(annot_df_oneGene) <- annot_df[,1]
annot_df_oneGene[,"UCSC_RefGene_Name"] <- gsub(";.*$","",annot_df[,"UCSC_RefGene_Name"])
annot_df_oneGene <- annot_df_oneGene[,c("CHR", "Infinium_Design_Type", "Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Name")]
head(annot_df_oneGene)

# aggregate
annot_df_oneGene <- annot_df_oneGene[!(annot_df_oneGene[,"UCSC_RefGene_Name"]== ""),]
#365860      4
BetaMatrix_KICH <- BetaMatrix_KICH[match(rownames(annot_df_oneGene),rownames(BetaMatrix_KICH)),]
#365860      66
BetaMatrix_KICH <- cbind(BetaMatrix_KICH,annot_df_oneGene[,"UCSC_RefGene_Name"])
BetaMatrix_KICH <- as.data.frame(BetaMatrix_KICH)
BetaMatrix_KICH[,1:66] <- apply(BetaMatrix_KICH[,1:66],2,as.numeric)
colnames(BetaMatrix_KICH)[67] <- "UCSC_RefGene_Name"


aggdata <-aggregate(BetaMatrix_KICH[,1:66], by= list(BetaMatrix_KICH$UCSC_RefGene_Name),
                    FUN=mean, na.rm=TRUE)
rownames(aggdata) <- aggdata[,1]
aggdata <- aggdata[,-1]
#20621    66
sum(complete.cases(aggdata))
#20348
aggdata <- aggdata[complete.cases(aggdata),]
aggdata <- apply(-aggdata,2,rank)
load("xCell.data.rda")
# 10808 genes and 489 signatures
length(intersect(xCell.data$genes,rownames(aggdata)))
#10694
source("xCell_custom.R")
xCellScores<- xCellAnalysis(aggdata, rnaseq = FALSE)
save(xCellScores, file = "xCellScores_KICH_methyl_mean.RData")


### alternatively, take the PCA

library(RTCGA)
path_KICH <- list.files(path = "/Users/junesong/Desktop/causal inference/gdac.broadinstitute.org_KICH.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/MethylDataset", full.names = TRUE, recursive = TRUE)
MethylTab_KICH <- readTCGA(path_KICH,dataType = "methylation")
BetaMatrix_KICH <- matrix(NA,485577,66)
BetaMatrix_KICH <- t(MethylTab_KICH[,-1])
colnames(BetaMatrix_KICH) <- as.character(MethylTab_KICH[,1])
rownames(BetaMatrix_KICH) <- colnames(MethylTab_KICH)[-1]
## dim: 485577 66

manifest <- read.csv("HumanMethylation450_15017482_v1-2.csv", header = T, skip = 7)
# 486428     33
annot <- manifest[match(rownames(BetaMatrix_KICH),manifest[,1]),]
# 485577     33
annot_df <- as.data.frame(annot)
head(annot_df[,"UCSC_RefGene_Name"],100)

annot_df_oneGene <- annot_df
rownames(annot_df_oneGene) <- annot_df[,1]
annot_df_oneGene[,"UCSC_RefGene_Name"] <- gsub(";.*$","",annot_df[,"UCSC_RefGene_Name"])
annot_df_oneGene <- annot_df_oneGene[,c("CHR", "Infinium_Design_Type", "Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Name")]
head(annot_df_oneGene)

annot_df_oneGene <- annot_df_oneGene[!(annot_df_oneGene[,"UCSC_RefGene_Name"]== ""),]
#365860      4
BetaMatrix_KICH <- BetaMatrix_KICH[match(rownames(annot_df_oneGene),rownames(BetaMatrix_KICH)),]
#365860      66
BetaMatrix_KICH <- cbind(BetaMatrix_KICH,annot_df_oneGene[,"UCSC_RefGene_Name"])
BetaMatrix_KICH <- as.data.frame(BetaMatrix_KICH)
BetaMatrix_KICH[,1:66] <- apply(BetaMatrix_KICH[,1:66],2,as.numeric)
colnames(BetaMatrix_KICH)[67] <- "UCSC_RefGene_Name"

# prcomp








############ KICH RNA expression
library(RTCGA)
checkTCGA('Dates')
datInfo <- checkTCGA("DataSets","KICH", date = "2016-01-28")
downloadTCGA(cancerTypes = "KICH", dataSet = "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level", 
             destDir = "/Users/junesong/Desktop/causal inference", date = "2016-01-28")
path_KICH <- list.files(path = "/Users/junesong/Desktop/causal inference/gdac.broadinstitute.org_KICH.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/RNADataset", full.names = TRUE, recursive = TRUE)
RNATab_KICH <- readTCGA(path_KICH,dataType = "rnaseq")
rownames(RNATab_KICH) <- RNATab_KICH[,1]
RNATab_KICH <- RNATab_KICH[,-1]
RNATab_KICH <- t(RNATab_KICH)

sample_shared <- intersect(substr(as.character(MethylTab_KICH[,1]),1,19),substr(colnames(RNATab_KICH),1,19))
colnames(RNATab_KICH) <- substr(colnames(RNATab_KICH),1,19)
RNATab_KICH <- RNATab_KICH[,sample_shared]
# 20531    66

genes_all <- rownames(RNATab_KICH)
genes_all <- gsub("\\|.*","", genes_all)
# 20502 unqiue
RNATab_KICH <- RNATab_KICH[-(1:29),]
# 20502    66
rownames(RNATab_KICH) <- genes_all[-(1:29)]


RNA_KICH_rank <- apply(RNATab_KICH,2,rank)

load("xCell.data.rda")
# 10808 genes and 489 signatures
length(intersect(xCell.data$genes,rownames(RNA_KICH_rank)))
## 10808 shared
source("xCell_custom.R")
xCellScores<- xCellAnalysis(RNA_KICH_rank, rnaseq = TRUE)
## 67 66
save(xCellScores, file = "xCellScores_KICH_RNA.RData")

xCellScores <- t(xCellScores)
rownames(xCellScores) <- substr(rownames(xCellScores), 1, 10)

# plots
boxplot(xCellScores, las = 2,cex.axis = 0.5)

load("xCellScores_KICH_methyl.RData")
#load("xCellScores_KICH_methyl_mean.RData")
xCellScores <- t(xCellScores)
scores_methyl <- xCellScores

load("xCellScores_KICH_RNA.RData")
xCellScores <- t(xCellScores)
scores_RNA <- xCellScores

df <- rbind(scores_methyl, scores_RNA)
df_plot <- matrix(NA,8448,3)
df_plot[,1] <- as.vector(as.matrix(df[,1:64]))
method <- c(rep("methyl",66), rep("RNA",66))
df_plot[,2] <- rep(method,64)
df_plot[,3] <- rep(colnames(df)[1:64],each = 132)
colnames(df_plot) <- c("scores", "method","CellType")
df_plot <- as.data.frame(df_plot)
df_plot[,1] <- as.numeric(as.character(df_plot[,1]))

library(ggplot2)
ggplot(df_plot) + geom_boxplot(aes(CellType ,scores,color= method)) +
  ylab('scores') + theme(axis.text.x = element_text(angle = 90, hjust = 1))



## check average correlation scores 
load("xCellScores_KICH_methyl.RData")
xCellScores <- t(xCellScores)
scores_methyl <- xCellScores

load("xCellScores_KICH_RNA.RData")
xCellScores <- t(xCellScores)
scores_RNA <- xCellScores

load("xCellScores_KICH_methyl_mean.RData")
xCellScores <- t(xCellScores)
scores_methyl_mean <- xCellScores

corr <- rep(NA, 66)
# for (i in 1:66){
#   corr[i] <- cor(scores_methyl[i,1:64],scores_RNA[i,1:64], method = "spearman")
# }

for (i in 1:66){
  corr[i] <- cor(scores_methyl[i,1:64],scores_RNA[i,1:64])
}

for (i in 1:66){
  corr[i] <- cor(scores_methyl_mean[i,1:64],scores_RNA[i,1:64])
}

for (i in 1:66){
  corr[i] <- cor(scores_methyl_mean[i,1:64],scores_methyl[i,1:64])
}
