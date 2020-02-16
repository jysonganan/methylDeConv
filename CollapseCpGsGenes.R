# collapse CpGs into genes
library(IlluminaHumanMethylation450kmanifest)
data(IlluminaHumanMethylation450kmanifest)
library(minfi)

df <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = c("I", "II", "Control",
                              "I-Green", "I-Red", "SnpI", "SnpII"))
class(df)
getManifestInfo(IlluminaHumanMethylation450kmanifest, type = c("nLoci", "locusNames"))

head(getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "I"))
head(IlluminaHumanMethylation450kmanifest@data$TypeI)
head(IlluminaHumanMethylation450kmanifest@data$TypeII)
head(IlluminaHumanMethylation450kmanifest@data$TypeControl)


getProbeInfo(IlluminaHumanMethylation450kmanifest) #DataFrame with 135476 rows and 8 columns
getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "II") #DataFrame with 350036 rows and 4 columns
# or
IlluminaHumanMethylation450kmanifest@data$TypeI #DataFrame with 135476 rows and 8 columns
IlluminaHumanMethylation450kmanifest@data$TypeII #DataFrame with 350036 rows and 4 columns

IlluminaHumanMethylation450kmanifest@data$TypeControl #DataFrame with 850 rows and 4 columns

getManifestInfo(IlluminaHumanMethylation450kmanifest) #485512
getManifestInfo(IlluminaHumanMethylation450kmanifest, type = "locusNames") #chracter with length 485512
getManifest(IlluminaHumanMethylation450kmanifest)

getControlAddress(IlluminaHumanMethylation450kmanifest, controlType = "NORM_A")
getControlAddress(IlluminaHumanMethylation450kmanifest, controlType = "NORM_C")
getControlAddress(IlluminaHumanMethylation450kmanifest, controlType = "NORM_G")
getControlAddress(IlluminaHumanMethylation450kmanifest, controlType = "NORM_T")





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




library(IlluminaHumanMethylation450k.db)
CpG_annotation <- as.list(IlluminaHumanMethylation450kSYMBOL[mappedkeys(IlluminaHumanMethylation450kSYMBOL)])

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










