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
BetaMatrix_maxvar <- CollapseCpGsGenes(BetaMatrix, method = "maxvar", include = "all")
BetaMatrix_pca <- CollapseCpGsGenes(BetaMatrix, method = "PCA", include = "all")

BetaMatrix_average_nona <- apply(BetaMatrix_average,1,function(x){x[which(is.na(x)==TRUE)] <- mean(na.omit(x)); return(x)})
BetaMatrix_average_nona <- t(BetaMatrix_average_nona)
BetaMatrix_maxvar_nona <- apply(BetaMatrix_maxvar,1,function(x){x[which(is.na(x)==TRUE)] <- mean(na.omit(x)); return(x)})
BetaMatrix_maxvar_nona <- t(BetaMatrix_maxvar_nona)
BetaMatrix_pca_nona <- apply(BetaMatrix_pca,1,function(x){x[which(is.na(x)==TRUE)] <- mean(na.omit(x)); return(x)})
BetaMatrix_pca_nona <- t(BetaMatrix_pca_nona)

sig_matrix <- read.table("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/LM22.txt",header=T,sep="\t",row.names=1,check.names=F)

source('/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/CIBERSORT_custom.R')

res_collapse_average_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_average_nona, absolute = TRUE)
res_collapse_maxvar_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_maxvar_nona, absolute = TRUE)
res_collapse_pca_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_pca_nona, absolute = TRUE)


save("res_collapse_average_cibersort", "res_collapse_maxvar_cibersort", "res_collapse_maxvar_cibersort",
     file = "/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/KICH_res_collapse_cibersort_noreverse.RData")


