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