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
#dat_450k <- pickCompProbes(preprocessQuantile(FlowSorted.Blood.450k), cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),
                      probeSelect = "both",numProbes = 100)
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
