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










BetaMatrix_nona <- BetaMatrix[complete.cases(BetaMatrix),]
res_methyl <- epidish(BetaMatrix_nona, centDHSbloodDMC.m,  method = "CBS")




res_RNA <- cbind(PAAD_res_RNA_cibersort[,1]+PAAD_res_RNA_cibersort[,2],
                 PAAD_res_RNA_cibersort[,11]+PAAD_res_RNA_cibersort[,12],
                 PAAD_res_RNA_cibersort[,5]+PAAD_res_RNA_cibersort[,6]+PAAD_res_RNA_cibersort[,7],
                 PAAD_res_RNA_cibersort[,4],
                 PAAD_res_RNA_cibersort[,13],
                 PAAD_res_RNA_cibersort[,22],
                 PAAD_res_RNA_cibersort[,21])

res_RNA_normalized <- apply(res_RNA,1,function(x){return(x/sum(x))})
res_RNA_normalized <- t(res_RNA_normalized)




corr <- rep(NA,183)
for (i in 1:183){
   corr[i] <- cor(res_RNA_normalized[i,],res_methyl$est[i,],method = "spearman")
 }


# # within each cell type, across samples

corr <- rep(NA,7)
for (i in 1:7){
  corr[i] <- cor(res_RNA_normalized[,i],res_methyl$est[,i], use="complete.obs", method = "spearman")
}





### collapsing corr

corr <- rep(NA,22)
for (i in 1:22){
  corr[i] <- cor(PAAD_res_RNA_cibersort[,i],res_collapse_pca_cibersort[,i], use="complete.obs", method = "spearman")
}


corr <- rep(NA, 183)
for (i in 1:183){
  corr[i] <- cor(PAAD_res_RNA_cibersort[i,],res_collapse_average_cibersort[i,],method = "spearman")
}





res_RNA <- cbind(BRCA_res_RNA_cibersort[,1]+BRCA_res_RNA_cibersort[,2],
                 BRCA_res_RNA_cibersort[,11]+BRCA_res_RNA_cibersort[,12],
                 BRCA_res_RNA_cibersort[,5]+BRCA_res_RNA_cibersort[,6]+BRCA_res_RNA_cibersort[,7],
                 BRCA_res_RNA_cibersort[,4],
                 BRCA_res_RNA_cibersort[,13],
                 BRCA_res_RNA_cibersort[,22],
                 BRCA_res_RNA_cibersort[,21])

res_RNA_normalized <- apply(res_RNA,1,function(x){return(x/sum(x))})
res_RNA_normalized <- t(res_RNA_normalized)


 #
res_collapse <- cbind(res_collapse_average_cibersort[,1]+res_collapse_average_cibersort[,2],
                      res_collapse_average_cibersort[,11]+res_collapse_average_cibersort[,12],
                      res_collapse_average_cibersort[,5]+res_collapse_average_cibersort[,6]+res_collapse_average_cibersort[,7],
                      res_collapse_average_cibersort[,4],
                      res_collapse_average_cibersort[,13],
                      res_collapse_average_cibersort[,22],
                      res_collapse_average_cibersort[,21])


#maxvar
res_collapse <- cbind(res_collapse_maxvar_cibersort[,1]+res_collapse_maxvar_cibersort[,2],
                      res_collapse_maxvar_cibersort[,11]+res_collapse_maxvar_cibersort[,12],
                      res_collapse_maxvar_cibersort[,5]+res_collapse_maxvar_cibersort[,6]+res_collapse_maxvar_cibersort[,7],
                      res_collapse_maxvar_cibersort[,4],
                      res_collapse_maxvar_cibersort[,13],
                      res_collapse_maxvar_cibersort[,22],
                      res_collapse_maxvar_cibersort[,21])

# pca
res_collapse <- cbind(res_collapse_pca_cibersort[,1]+res_collapse_pca_cibersort[,2],
                      res_collapse_pca_cibersort[,11]+res_collapse_pca_cibersort[,12],
                      res_collapse_pca_cibersort[,5]+res_collapse_pca_cibersort[,6]+res_collapse_pca_cibersort[,7],
                      res_collapse_pca_cibersort[,4],
                      res_collapse_pca_cibersort[,13],
                      res_collapse_pca_cibersort[,22],
                      res_collapse_pca_cibersort[,21])





res_collapse_normalized <- apply(res_collapse,1,function(x){return(x/sum(x))})
res_collapse_normalized <- t(res_collapse_normalized)

corr <- rep(NA,7)
for (i in 1:7){
  corr[i] <- cor(res_RNA_normalized[,i],res_collapse_normalized[,i], use="complete.obs", method = "spearman")
}



