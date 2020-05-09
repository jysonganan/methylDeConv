disease <- "BRCA"
methyl_path <- paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_methyl/", disease, "_methyl.RData")
mRNA_path <- paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mrna/", disease, "_mrna.RData")
load(methyl_path)
load(mRNA_path)
sample_shared <- intersect(substr(colnames(BetaMatrix),1,19),substr(colnames(mRNAMatrix),1,19))
colnames(BetaMatrix) <- substr(colnames(BetaMatrix),1,19)
colnames(mRNAMatrix) <- substr(colnames(mRNAMatrix),1,19)
BetaMatrix <- BetaMatrix[,sample_shared]  #485577    253
mRNAMatrix <- mRNAMatrix[,sample_shared]  #17814   253

### impute the missing values in mRNA matrix
mRNAMatrix_nona <- apply(mRNAMatrix,1,function(x){x[which(is.na(x)==TRUE)] <- mean(na.omit(x)); return(x)})
mRNAMatrix_nona <- t(mRNAMatrix_nona)


sig_matrix <- read.table("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/LM22.txt",header=T,sep="\t",row.names=1,check.names=F)

source('/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/CIBERSORT_custom.R')

## CIBERSORT on the mRNA data
res_mRNA_cibersort <- CIBERSORT(sig_matrix, mRNAMatrix_nona)



### collapsing
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

BetaMatrix_average_TSS200 <- CollapseCpGsGenes(BetaMatrix, include = "TSS200")
BetaMatrix_maxvar_TSS200 <- CollapseCpGsGenes(BetaMatrix, method = "maxvar", include = "TSS200")
BetaMatrix_pca_TSS200 <- CollapseCpGsGenes(BetaMatrix, method = "PCA", include = "TSS200")

BetaMatrix_average_TSS200_nona <- apply(BetaMatrix_average_TSS200,1,function(x){x[which(is.na(x)==TRUE)] <- mean(na.omit(x)); return(x)})
BetaMatrix_average_TSS200_nona  <- t(BetaMatrix_average_TSS200_nona)
BetaMatrix_maxvar_TSS200_nona <- apply(BetaMatrix_maxvar_TSS200,1,function(x){x[which(is.na(x)==TRUE)] <- mean(na.omit(x)); return(x)})
BetaMatrix_maxvar_TSS200_nona <- t(BetaMatrix_maxvar_TSS200_nona)
BetaMatrix_pca_TSS200_nona <- apply(BetaMatrix_pca_TSS200,1,function(x){x[which(is.na(x)==TRUE)] <- mean(na.omit(x)); return(x)})
BetaMatrix_pca_TSS200_nona <- t(BetaMatrix_pca_TSS200_nona)

BetaMatrix_average_TSS <- CollapseCpGsGenes(BetaMatrix, include = "TSS200&TSS1500")
BetaMatrix_maxvar_TSS <- CollapseCpGsGenes(BetaMatrix, method = "maxvar", include = "TSS200&TSS1500")
BetaMatrix_pca_TSS <- CollapseCpGsGenes(BetaMatrix, method = "PCA", include = "TSS200&TSS1500")

BetaMatrix_average_TSS_nona <- apply(BetaMatrix_average_TSS,1,function(x){x[which(is.na(x)==TRUE)] <- mean(na.omit(x)); return(x)})
BetaMatrix_average_TSS_nona  <- t(BetaMatrix_average_TSS_nona)
BetaMatrix_maxvar_TSS_nona <- apply(BetaMatrix_maxvar_TSS,1,function(x){x[which(is.na(x)==TRUE)] <- mean(na.omit(x)); return(x)})
BetaMatrix_maxvar_TSS_nona <- t(BetaMatrix_maxvar_TSS_nona)
BetaMatrix_pca_TSS_nona <- apply(BetaMatrix_pca_TSS,1,function(x){x[which(is.na(x)==TRUE)] <- mean(na.omit(x)); return(x)})
BetaMatrix_pca_TSS_nona <- t(BetaMatrix_pca_TSS_nona)



res_collapse_average_all_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_average_nona)
res_collapse_maxvar_all_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_maxvar_nona)
res_collapse_pca_all_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_pca_nona)

res_collapse_average_TSS200_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_average_TSS200_nona)
res_collapse_maxvar_TSS200_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_maxvar_TSS200_nona)
res_collapse_pca_TSS200_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_pca_TSS200_nona)

res_collapse_average_TSS_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_average_TSS_nona)
res_collapse_maxvar_TSS_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_maxvar_TSS_nona)
res_collapse_pca_TSS_cibersort <- CIBERSORT(sig_matrix, BetaMatrix_pca_TSS_nona)

save("res_mRNA_cibersort", "res_collapse_average_all_cibersort", "res_collapse_maxvar_all_cibersort",
     "res_collapse_pca_all_cibersort", "res_collapse_average_TSS200_cibersort", "res_collapse_maxvar_TSS200_cibersort",
     "res_collapse_pca_TSS200_cibersort", "res_collapse_average_TSS_cibersort","res_collapse_maxvar_TSS_cibersort",
     "res_collapse_pca_TSS_cibersort", file = "/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/brcaCIBERSORTmRNAvsMethylCollapsed.RData")


##### check results
## within cell type, across samples
# corr <- rep(NA,22)
# names(corr) <- colnames(res_mRNA_cibersort)[1:22]
# for (i in 1:22){
#   corr[i] <- cor(res_mRNA_cibersort[,i],res_collapse_average_TSS200_cibersort[,i], method = "spearman")
# }
# 
# ## within sample, across cell types
# corr <- rep(NA,253)
# for (i in 1:253){
#   corr[i] <- cor(res_mRNA_cibersort[i,1:22],res_collapse_pca_TSS_cibersort[i,1:22], method = "spearman")
# }
# 
# 
# 


# res_cbs <- noNa_res2_CBS[,4:10]
# res_rpc <- noNa_res2_RPC[,4:10]
# res_cp <- noNa_res2_CP[,4:10]
# res_rna_cibersort <- cbind(res_mRNA_cibersort[,1]+res_mRNA_cibersort[,2],
#                            res_mRNA_cibersort[,11]+res_mRNA_cibersort[,12],
#                            res_mRNA_cibersort[,5]+res_mRNA_cibersort[,6]+res_mRNA_cibersort[,7],
#                            res_mRNA_cibersort[,4],
#                            res_mRNA_cibersort[,13],
#                            res_mRNA_cibersort[,22],
#                            res_mRNA_cibersort[,21])
# 
# colnames(res_rna_cibersort) <- colnames(res_cbs)
# 
# res_collapse_cibersort <- cbind(res_collapse_average_all_cibersort[,1]+res_collapse_average_all_cibersort[,2],
#                            res_collapse_average_all_cibersort[,11]+res_collapse_average_all_cibersort[,12],
#                            res_collapse_average_all_cibersort[,5]+res_collapse_average_all_cibersort[,6]+res_collapse_average_all_cibersort[,7],
#                            res_collapse_average_all_cibersort[,4],
#                            res_collapse_average_all_cibersort[,13],
#                            res_collapse_average_all_cibersort[,22],
#                            res_collapse_average_all_cibersort[,21])
# 
# colnames(res_collapse_cibersort) <- colnames(res_cbs)
# 
# head(res_cbs)
# head(res_rpc)
# head(res_cp)
# head(res_rna_cibersort)
# head(res_collapse_cibersort)
# 
# res_cbs_normalized <- apply(res_cbs,1,function(x){return(x/sum(x))})
# res_cbs_normalized <- t(res_cbs_normalized)
# 
# res_rpc_normalized <- apply(res_rpc,1,function(x){return(x/sum(x))})
# res_rpc_normalized <- t(res_rpc_normalized)
# 
# res_cp_normalized <- apply(res_cp,1,function(x){return(x/sum(x))})
# res_cp_normalized <- t(res_cp_normalized)
# 
# res_rna_cibersort_normalized <- apply(res_rna_cibersort,1,function(x){return(x/sum(x))})
# res_rna_cibersort_normalized<- t(res_rna_cibersort_normalized)
# 
# res_collapse_cibersort_normalized <- apply(res_collapse_cibersort,1,function(x){return(x/sum(x))})
# res_collapse_cibersort_normalized <- t(res_collapse_cibersort_normalized)
# 
# # within each cell type, across samples
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_rna_cibersort_normalized[,i],res_cbs_normalized[,i], use="complete.obs", method = "spearman")
# }
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_rna_cibersort_normalized[,i],res_rpc_normalized[,i], use="complete.obs",method = "spearman")
# }
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_rna_cibersort_normalized[,i],res_cp_normalized[,i], use="complete.obs",method = "spearman")
# }
# 
# 
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_collapse_cibersort_normalized[,i],res_cbs_normalized[,i], use="complete.obs",method = "spearman")
# }
# corr
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_collapse_cibersort_normalized[,i],res_rpc_normalized[,i], use="complete.obs",method = "spearman")
# }
# corr
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_collapse_cibersort_normalized[,i],res_cp_normalized[,i], use="complete.obs",method = "spearman")
# }
# corr
# 
# 
# 
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_cp_normalized[,i],res_cbs_normalized[,i], use="complete.obs", method = 'spearman')
# }
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_rpc_normalized[,i],res_cbs_normalized[,i], use="complete.obs", method = 'spearman')
# }
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_cp_normalized[,i],res_rpc_normalized[,i], use="complete.obs", method = 'spearman')
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # within each cell type, across samples
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_rna_cibersort_normalized[,i],res_CBS_IC$est[,i], use="complete.obs", method = "spearman")
# }
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_rna_cibersort_normalized[,i],res_RPC_IC$est[,i], use="complete.obs",method = "spearman")
# }
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_rna_cibersort_normalized[,i],res_CP_IC$est[,i], use="complete.obs",method = "spearman")
# }
# 
# 
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_collapse_cibersort_normalized[,i],res_CBS_IC$est[,i], use="complete.obs",method = "spearman")
# }
# corr
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_collapse_cibersort_normalized[,i],res_RPC_IC$est[,i], use="complete.obs",method = "spearman")
# }
# corr
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_collapse_cibersort_normalized[,i],res_CP_IC$est[,i], use="complete.obs",method = "spearman")
# }
# corr
# 
# 
# 
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_CP_IC$est[,i],res_CBS_IC$est[,i])
# }
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_RPC_IC$est[,i],res_CBS_IC$est[,i])
# }
# 
# corr <- rep(NA,7)
# names(corr) <- colnames(res_cbs)
# for (i in 1:7){
#   corr[i] <- cor(res_CP_IC$est[,i],res_RPC_IC$est[,i])
# }
# 
# 
# # corr <- rep(NA,253)
# # for (i in 1:253){
# #   corr[i] <- cor(res_mRNA_cibersort[i,1:22],res_collapse_pca_TSS_cibersort[i,1:22], method = "spearman")
# # }






# corr <- rep(NA,253)
# for (i in 1:253){
#   corr[i] <- cor(res_rna_cibersort_normalized[i,],res_collapse_cibersort_normalized[i,])
# }
# 
# 
# corr <- rep(NA,253)
# for (i in 1:253){
#   corr[i] <- cor(res_rna_cibersort_normalized[i,],res_CBS_IC$est[i,])
# }
# mean(corr)
# 
# corr <- rep(NA,253)
# for (i in 1:253){
#   corr[i] <- cor(res_rna_cibersort_normalized[i,],res_CBS_IC$est[i,],method = "spearman")
# }
# mean(corr)
# 
# 
# corr <- rep(NA,253)
# for (i in 1:253){
#   corr[i] <- cor(res_collapse_cibersort_normalized[i,],res_RPC_IC$est[i,],method = "spearman")
# }
# mean(corr)
# 
#
