### compare 450k, EPIC blood ref

library(FlowSorted.Blood.450k)
dat_450k <- pickCompProbes(preprocessQuantile(FlowSorted.Blood.450k), cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),
                      probeSelect = "both")

library(FlowSorted.Blood.EPIC)
load("/Users/junesong/Desktop/causal inference/CellProportion/methylDeconv_EPICdata/FlowSorted.Blood.EPIC.RData")
GRset_EPIC <- preprocessNoob(FlowSorted.Blood.EPIC)
dat_EPIC<- pickCompProbes(GRset_EPIC, cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"), probeSelect = "both")





########################## CIBERSORT Signatures
cibersort_sig_matrix <- read.table("LM22.txt",header=T,sep="\t",row.names=1,check.names=F)
# 547 genes
## Check back genes for 450k selected signature probes (600)
# manifest annotations
manifest <- read.csv("HumanMethylation450_15017482_v1-2.csv", header = T, skip = 7)
annot <- manifest[match(dat_450k[["trainingProbes"]],manifest[,1]),] 
annot_df <- as.data.frame(annot)

annot_df_oneGene <- annot_df
rownames(annot_df_oneGene) <- annot_df[,1]
annot_df_oneGene[,"UCSC_RefGene_Name"] <- gsub(";.*$","",annot_df[,"UCSC_RefGene_Name"])
annot_df_oneGene[,"UCSC_RefGene_Group"] <- gsub(";.*$","",annot_df[,"UCSC_RefGene_Group"])
annot_df_oneGene[1:5,c("UCSC_RefGene_Name","UCSC_RefGene_Group")]

intersect(rownames(cibersort_sig_matrix),annot_df_oneGene[,"UCSC_RefGene_Name"]) 
## length: 36; if 450k for each cell type 100
intersect(rownames(cibersort_sig_matrix),annot_df_oneGene[,"UCSC_RefGene_Name"])
## length: 53; if 450k for each cell type 200

### check for the details of genes
## enrichment of the 600 genes (corresponds to 600 signature CpGs)

## map cibersort signatures to CpGs
annot_df <- as.data.frame(manifest)

annot_df_oneGene <- annot_df
rownames(annot_df_oneGene) <- annot_df[,1]
annot_df_oneGene[,"UCSC_RefGene_Name"] <- gsub(";.*$","",annot_df[,"UCSC_RefGene_Name"])
annot_df_oneGene[,"UCSC_RefGene_Group"] <- gsub(";.*$","",annot_df[,"UCSC_RefGene_Group"])
annot_df_oneGene[1:5,c("UCSC_RefGene_Name","UCSC_RefGene_Group")]
cibersort_sig_cpgs <- rownames(annot_df_oneGene[annot_df_oneGene[,"UCSC_RefGene_Name"]%in%rownames(cibersort_sig_matrix),])
## length 8100

### use cibersort_sig_cpgs for methylation deconvolution
library(minfi)
library(GEOquery)
rgSet_77797 <- read.metharray.exp("GSE77797/idat")
geoMat_77797<- getGEO("GSE77797")
pD.all <- pData(geoMat_77797[[1]])
pD <- pD.all[, c("title", "geo_accession", "b cell (%):ch1", "cd4+ t cell (%):ch1", "cd8+ t cell (%):ch1",
                 "granulocyte (%):ch1", "monocyte (%):ch1", "natural killer cell (%):ch1")]
sampleNames(rgSet_77797) <- substr(sampleNames(rgSet_77797), 1, 10)
pD <- pD[sampleNames(rgSet_77797),]
pD <- as(pD, "DataFrame")
pData(rgSet_77797) <- pD
grSet_77797 <- preprocessQuantile(rgSet_77797)
betaMat_77797 <- getBeta(grSet_77797)
## 485512      18

res1 <- MethylDeconv(betaMat_77797, method = "Houseman")
res2 <- MethylDeconv(betaMat_77797, method = "RPC")
res3 <- MethylDeconv(betaMat_77797, method = "CBS")
### use top 100 * 2 *6
#dat_450k <- pickCompProbes(preprocessQuantile(FlowSorted.Blood.450k), cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"), probeSelect = "both",numProbes = 100)
#res1 <- MethylDeconv(betaMat_77797, method = "Houseman",custom_probes = dat_450k[[2]])
#res2 <- MethylDeconv(betaMat_77797, method = "RPC",custom_probes = dat_450k[[2]])
#res3 <- MethylDeconv(betaMat_77797, method = "CBS",custom_probes = dat_450k[[2]])

### Use cibersort cpgs instead
res1_cbs <- MethylDeconv(betaMat_77797, method = "Houseman",custom_probes = intersect(cibersort_sig_cpgs,dat_450k[[2]]))
res2_cbs <- MethylDeconv(betaMat_77797, method = "RPC", custom_probes = intersect(cibersort_sig_cpgs,dat_450k[[2]]))
res3_cbs <- MethylDeconv(betaMat_77797, method = "CBS", custom_probes = intersect(cibersort_sig_cpgs,dat_450k[[2]]))


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
facs_prop <- as.matrix(facs_77797_prop)

##### analyze Collaspe on Benchmark data (GEO with FACS)
corr <- rep(NA,18)
for (i in 1:18){
  corr[i] <- cor(res1[i,], as.numeric(facs_prop[i,]),method = "spearman")
}
mean(corr)

corr <- rep(NA,18)
for (i in 1:18){
  corr[i] <- cor(res2[i,], as.numeric(facs_prop[i,]),method = "spearman")
}
mean(corr)

corr <- rep(NA,18)
for (i in 1:18){
  corr[i] <- cor(res3[i,], as.numeric(facs_prop[i,]),method = "spearman")
}
mean(corr)


corr <- rep(NA,18)
for (i in 1:18){
  corr[i] <- cor(res1_cbs[i,], as.numeric(facs_prop[i,]),method = "spearman")
}
mean(corr)

corr <- rep(NA,18)
for (i in 1:18){
  corr[i] <- cor(res2_cbs[i,], as.numeric(facs_prop[i,]),method = "spearman")
}
mean(corr)

corr <- rep(NA,18)
for (i in 1:18){
  corr[i] <- cor(res3_cbs[i,], as.numeric(facs_prop[i,]),method = "spearman")
}
mean(corr)

corr <- rep(NA,6)
for (i in 1:6){
  corr[i] <- cor(res1[,i], as.numeric(facs_prop[,i]),method = "spearman")
}
corr

corr <- rep(NA,6)
for (i in 1:6){
  corr[i] <- cor(res2[,i], as.numeric(facs_prop[,i]),method = "spearman")
}
corr

corr <- rep(NA,6)
for (i in 1:6){
  corr[i] <- cor(res3[,i], as.numeric(facs_prop[,i]),method = "spearman")
}
corr


corr <- rep(NA,6)
for (i in 1:6){
  corr[i] <- cor(res1_cbs[,i], as.numeric(facs_prop[,i]),method = "spearman")
}
corr

corr <- rep(NA,6)
for (i in 1:6){
  corr[i] <- cor(res2_cbs[,i], as.numeric(facs_prop[,i]),method = "spearman")
}
corr

corr <- rep(NA,6)
for (i in 1:6){
  corr[i] <- cor(res3_cbs[,i], as.numeric(facs_prop[,i]),method = "spearman")
}
corr


########################## 450k vs EPIC Signatures
# B cells versus all else
a <- dat_450k[["tstatList"]][["Bcell"]][order(dat_450k[["tstatList"]][["Bcell"]]$statistic),]
b <- dat_EPIC[["tstatList"]][["Bcell"]][order(dat_EPIC[["tstatList"]][["Bcell"]]$statistic),]
length(intersect(rownames(head(a,50)),rownames(head(b,50))))
length(intersect(rownames(tail(a,50)),rownames(tail(b,50))))
length(intersect(rownames(head(a,100)),rownames(head(b,100))))
length(intersect(rownames(tail(a,100)),rownames(tail(b,100))))
# > dim(a)
# [1] 485512      3
# > dim(b)
# [1] 866091      3
# > length(intersect(rownames(a),rownames(b)))
# [1] 452567
a1 <- a[intersect(rownames(a),rownames(b)),]
b1 <- b[intersect(rownames(a),rownames(b)),]  
a1 <- a1[order(a1$statistic),]
b1 <- b1[order(b1$statistic),]
length(intersect(rownames(head(a1,50)),rownames(head(b1,50))))
length(intersect(rownames(tail(a1,50)),rownames(tail(b1,50))))
length(intersect(rownames(head(a1,100)),rownames(head(b1,100))))
length(intersect(rownames(tail(a1,100)),rownames(tail(b1,100))))

## check back the gene signatures
annot_df <- as.data.frame(manifest)

annot_df_oneGene <- annot_df
rownames(annot_df_oneGene) <- annot_df[,1]

annot_df_oneGene <- annot_df_oneGene[intersect(rownames(head(a,50)),rownames(head(b,50))),]
annot_df_oneGene[,"UCSC_RefGene_Name"] <- gsub(";.*$","",annot_df_oneGene[,"UCSC_RefGene_Name"])
annot_df_oneGene[,"UCSC_RefGene_Group"] <- gsub(";.*$","",annot_df_oneGene[,"UCSC_RefGene_Group"])
annot_df_oneGene[,c("UCSC_RefGene_Name","UCSC_RefGene_Group")]




### check whether the top CpGs in EPIC are those CpGs not among 450000 CpGs
length(intersect(rownames(head(b,50)),rownames(a)))
length(intersect(rownames(tail(b,50)),rownames(a)))
length(intersect(rownames(head(b,100)),rownames(a)))
length(intersect(rownames(tail(b,100)),rownames(a)))

length(intersect(rownames(head(a,50)),rownames(b)))
length(intersect(rownames(tail(a,50)),rownames(b)))
length(intersect(rownames(head(a,100)),rownames(b)))
length(intersect(rownames(tail(a,100)),rownames(b)))


a <- dat_450k[["tstatList"]][["CD8T"]][order(dat_450k[["tstatList"]][["CD8T"]]$statistic),]
b <- dat_EPIC[["tstatList"]][["CD8T"]][order(dat_EPIC[["tstatList"]][["CD8T"]]$statistic),]
length(intersect(rownames(head(a,50)),rownames(head(b,50))))
length(intersect(rownames(tail(a,50)),rownames(tail(b,50))))
length(intersect(rownames(head(a,100)),rownames(head(b,100))))
length(intersect(rownames(tail(a,100)),rownames(tail(b,100))))
# > dim(a)
# [1] 485512      3
# > dim(b)
# [1] 866091      3
# > length(intersect(rownames(a),rownames(b)))
# [1] 452567
a1 <- a[intersect(rownames(a),rownames(b)),]
b1 <- b[intersect(rownames(a),rownames(b)),]  
a1 <- a1[order(a1$statistic),]
b1 <- b1[order(b1$statistic),]
length(intersect(rownames(head(a1,50)),rownames(head(b1,50))))
length(intersect(rownames(tail(a1,50)),rownames(tail(b1,50))))
length(intersect(rownames(head(a1,100)),rownames(head(b1,100))))
length(intersect(rownames(tail(a1,100)),rownames(tail(b1,100))))


#### 









########################## 450k vs EPIC Signatures
# B cells versus all else
a <- dat_450k[["tstatList"]][["Bcell"]][order(dat_450k[["tstatList"]][["Bcell"]]$statistic),]
b <- dat_EPIC[["tstatList"]][["Bcell"]][order(dat_EPIC[["tstatList"]][["Bcell"]]$statistic),]

### check the rank (in EPIC) of the signatures selected based on 450k data
rank_EPIC <- rank(dat_EPIC[["tstatList"]][["Bcell"]][,1])
names(rank_EPIC) <- rownames(dat_EPIC[["tstatList"]][["Bcell"]])
#Check how the top (50/100 both up-/down-regulated) CpGs from 450k ref profiles ranked in EPIC profiles
rank_EPIC[c(rownames(head(a,50)),rownames(tail(a,50)))]
hist(rank_EPIC[rownames(head(a,50))])
mean(rank_EPIC[rownames(head(a,50))],na.rm = TRUE)
median(rank_EPIC[rownames(head(a,50))],na.rm = TRUE)
mean(rank_EPIC[rownames(tail(a,50))],na.rm = TRUE)
median(rank_EPIC[rownames(tail(a,50))],na.rm = TRUE)

rank_EPIC[c(rownames(head(a,100)),rownames(tail(a,100)))]
mean(rank_EPIC[rownames(head(a,100))],na.rm = TRUE)
median(rank_EPIC[rownames(head(a,100))],na.rm = TRUE)
mean(rank_EPIC[rownames(tail(a,100))],na.rm = TRUE)
median(rank_EPIC[rownames(tail(a,100))],na.rm = TRUE)

dat_EPIC[["tstatList"]][["Bcell"]][c(rownames(head(a,50)),rownames(tail(a,50))),]


### check the top 100 CpGs selected from 450k 
### check the top 100 CpGs selected from EPIC

## when they both applied as the selected CpGs on EPIC benchmark dataset with Flow.EPIC reference
library(minfi)
rgSet_112618 <- read.metharray.exp("GSE112618/idat")
library(GEOquery)
geoMat_112618 <- getGEO("GSE112618")
pD.all <- pData(geoMat_112618[[1]])


pD <- pD.all[, c("title", "geo_accession", "bcell proportion:ch1", "cd4t proportion:ch1", "cd8t proportion:ch1",
                 "granulocytes proportion:ch1", "monocytes proportion:ch1", "neutrophils proportion:ch1", 
                 "nk proportion:ch1","Sex:ch1")]

sampleNames(rgSet_112618) <- substr(sampleNames(rgSet_112618), 1, 10)
pD <- pD[sampleNames(rgSet_112618),]
pD <- as(pD, "DataFrame")
pData(rgSet_112618) <- pD
grSet_112618 <- preprocessNoob(rgSet_112618, dyeMethod = "single")
betaMat_112618 <- getBeta(grSet_112618)



#### preprocessNoob dyeMethod = "single"

library(FlowSorted.Blood.450k)
dat_450k <- pickCompProbes(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"), cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),
                           probeSelect = "both")

library(FlowSorted.Blood.EPIC)
load("/Users/junesong/Desktop/causal inference/CellProportion/methylDeconv_EPICdata/FlowSorted.Blood.EPIC.RData")
GRset_EPIC <- preprocessNoob(FlowSorted.Blood.EPIC, dyeMethod = "single")
dat_EPIC<- pickCompProbes(GRset_EPIC, cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"), probeSelect = "both")


# # B cells versus all else
# a <- dat_450k[["tstatList"]][["Bcell"]][order(dat_450k[["tstatList"]][["Bcell"]]$statistic),]
# b <- dat_EPIC[["tstatList"]][["Bcell"]][order(dat_EPIC[["tstatList"]][["Bcell"]]$statistic),]
# 
# # ### check the rank (in EPIC) of the signatures selected based on 450k data
# rank_EPIC <- rank(dat_EPIC[["tstatList"]][["Bcell"]][,1])
# names(rank_EPIC) <- rownames(dat_EPIC[["tstatList"]][["Bcell"]])
# #Check how the top (50/100 both up-/down-regulated) CpGs from 450k ref profiles ranked in EPIC profiles
# rank_EPIC[c(rownames(head(a,50)),rownames(tail(a,50)))]
# hist(rank_EPIC[rownames(head(a,50))])
# mean(rank_EPIC[rownames(head(a,50))],na.rm = TRUE)
# median(rank_EPIC[rownames(head(a,50))],na.rm = TRUE)
# mean(rank_EPIC[rownames(tail(a,50))],na.rm = TRUE)
# median(rank_EPIC[rownames(tail(a,50))],na.rm = TRUE)
# 
# rank_EPIC[c(rownames(head(a,100)),rownames(tail(a,100)))]
# mean(rank_EPIC[rownames(head(a,100))],na.rm = TRUE)
# median(rank_EPIC[rownames(head(a,100))],na.rm = TRUE)
# mean(rank_EPIC[rownames(tail(a,100))],na.rm = TRUE)
# median(rank_EPIC[rownames(tail(a,100))],na.rm = TRUE)
# 
# dat_EPIC[["tstatList"]][["Bcell"]][c(rownames(head(a,50)),rownames(tail(a,50))),]



probes_450k <- dat_450k$trainingProbes
probes_EPIC <- dat_EPIC$trainingProbes


res1_450k <- MethylDeconv_normalized_BloodEPIC(betaMat_112618, method = "Houseman", custom_probes = probes_450k)
res2_450k <- MethylDeconv_normalized_BloodEPIC(betaMat_112618, method = "RPC", custom_probes = probes_450k)
res3_450k <- MethylDeconv_normalized_BloodEPIC(betaMat_112618, method = "CBS", custom_probes = probes_450k)

res1_EPIC <- MethylDeconv_normalized_BloodEPIC(betaMat_112618, method = "Houseman",custom_probes = probes_EPIC)
res2_EPIC <- MethylDeconv_normalized_BloodEPIC(betaMat_112618, method = "RPC",custom_probes = probes_EPIC)
res3_EPIC <- MethylDeconv_normalized_BloodEPIC(betaMat_112618, method = "CBS", custom_probes = probes_EPIC)





pD_112618 <-pD.all[, c("title", "geo_accession", "bcell proportion:ch1", "cd4t proportion:ch1", "cd8t proportion:ch1",
                       "granulocytes proportion:ch1", "monocytes proportion:ch1", "neutrophils proportion:ch1", 
                       "nk proportion:ch1","Sex:ch1")]
pD_112618 <- pD_112618[sampleNames(rgSet_112618),]
facs_112618 <- pD_112618[,3:9]
colnames(facs_112618) <- c("Bcell", "CD4T", "CD8T","Gran","Mono","Neutro","NK")
facs_112618_prop <- cbind(as.numeric(facs_112618[,"CD8T"]), as.numeric(facs_112618[,"CD4T"]), 
                          as.numeric(facs_112618[,"NK"]),
                          as.numeric(facs_112618[,"Bcell"]),as.numeric(facs_112618[,"Mono"]),
                          as.numeric(facs_112618[,"Gran"])+as.numeric(facs_112618[,"Neutro"]))
facs_112618_prop <- as.data.frame(facs_112618_prop)
rownames(facs_112618_prop) = rownames(res1_450k)
colnames(facs_112618_prop) = colnames(res1_450k)
facs_112618_prop <- as.matrix(facs_112618_prop)
### check correlation within each sample
corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(res1_EPIC[i,],as.numeric(as.character(facs_112618_prop[i,])),method = "spearman")
}
mean(corr)

corr2 <- rep(NA,6)
for (i in 1:6){
  corr2[i] <-cor(res2_EPIC[i,],as.numeric(as.character(facs_112618_prop[i,])),method = "spearman")
}
mean(corr2)

corr3 <- rep(NA,6)
for (i in 1:6){
  corr3[i] <-cor(res3_EPIC[i,],as.numeric(as.character(facs_112618_prop[i,])),method = "spearman")
}
mean(corr3)
### within each cell type
corr <- rep(NA, 6)
for (i in 1:6){
  corr[i] <-cor(res3_EPIC[,i],as.numeric(as.character(facs_112618_prop[,i])),method = "spearman")
}
corr




##### another benchmark data
### GSE122126
library(GEOquery)
library(minfi)
getGEOSuppFiles("GSE122126")
untar("GSE122126/GSE122126_RAW.tar", exdir = "GSE122126/idat")
head(list.files("GSE122126/idat", pattern = "idat"))

idatFiles <- list.files("GSE122126/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)
rgSet <- read.metharray.exp("GSE122126/idat", force = TRUE)


########### Enrichment Analysis of overlap signatures
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(missMethyl)
gst <- gometh(sig.cpg=dat_450k[["trainingProbes"]], all.cpg=rownames(dat_450k[["compTable"]]), collection="GO")
gst <- gometh(sig.cpg=dat_450k[["trainingProbes"]], all.cpg=rownames(dat_450k[["compTable"]]), collection="KEGG")
topGSA(gst)
#change the array type, the array.type argument can be specified as either “450K” or “EPIC”. The default is “450K”.
# GO or KEGG
gst <- gometh(sig.cpg=dat_EPIC[["trainingProbes"]], all.cpg=rownames(dat_EPIC[["compTable"]]), collection="GO", array.type = "EPIC")
gst <- gometh(sig.cpg=dat_EPIC[["trainingProbes"]], all.cpg=rownames(dat_EPIC[["compTable"]]), collection="KEGG", array.type = "EPIC")
topGSA(gst)
