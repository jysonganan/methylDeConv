disease <- "KICH"
methyl_path <- paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_methyl/", disease, "_methyl.RData")
RNA_path <- paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_rna/", disease, "_rna.RData")
load(methyl_path)
load(RNA_path)
sample_shared <- intersect(substr(colnames(BetaMatrix),1,19),substr(colnames(RNAMatrix),1,19))
colnames(BetaMatrix) <- substr(colnames(BetaMatrix),1,19)
colnames(RNAMatrix) <- substr(colnames(RNAMatrix),1,19)
BetaMatrix <- BetaMatrix[,sample_shared]  
RNAMatrix <- RNAMatrix[,sample_shared] 

print(dim(BetaMatrix))
print(dim(RNAMatrix))

source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/CollapseCpGsGenes.R")
BetaMatrix_average <- CollapseCpGsGenes(BetaMatrix, include = "all")
BetaMatrix_average_nearest <- CollapseCpGsGenes(BetaMatrix, NearestToCpG = TRUE)
BetaMatrix_average_nearest_TSS <- CollapseCpGsGenes(BetaMatrix, NearestToCpG = TRUE, NearestMethod = "TSS")

genes_all <- rownames(RNAMatrix)
genes_all <- gsub("\\|.*","", genes_all)
# 20502 unqiue
RNAMatrix <- RNAMatrix[-(1:29),]
# 20502    66
rownames(RNAMatrix) <- genes_all[-(1:29)]


sig_matrix <- read.table("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/LM22.txt",header=T,sep="\t",row.names=1,check.names=F)


shared_genes_average <- intersect(intersect(rownames(sig_matrix),rownames(RNAMatrix)), rownames(BetaMatrix_average))
#502
shared_genes_average_nearest <- intersect(intersect(rownames(sig_matrix),rownames(RNAMatrix)), rownames(BetaMatrix_average_nearest))
shared_genes_average_nearest_TSS  <- intersect(intersect(rownames(sig_matrix),rownames(RNAMatrix)), rownames(BetaMatrix_average_nearest_TSS))

BetaMatrix_average <- BetaMatrix_average[shared_genes_average,]
BetaMatrix_average_nearest <- BetaMatrix_average_nearest[shared_genes_average_nearest,]
BetaMatrix_average_nearest_TSS <- BetaMatrix_average_nearest_TSS[shared_genes_average_nearest_TSS,]

save("BetaMatrix_average", "BetaMatrix_average_nearest","BetaMatrix_average_nearest_TSS","RNAMatrix",file = "KICHforGP.RData")


corr <- rep(NA, length(shared_genes))
for (i in 1:length(shared_genes)){
  gene <- shared_genes[i]
  corr[i] <- cor(RNAMatrix[gene,],as.numeric(BetaMatrix_average[gene,]))
}


corr <- rep(NA, dim(BetaMatrix)[2])
for (i in 1:dim(BetaMatrix)[2]){
  corr[i] <- cor(RNAMatrix[shared_genes,i],BetaMatrix_average[shared_genes,i],method = "spearman")
}






