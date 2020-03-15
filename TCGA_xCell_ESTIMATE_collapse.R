disease <- "KICH"
methyl_path <- paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_methyl/", disease, "_methyl.RData")
RNA_path <- paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_rna/", disease, "_rna.RData")
load(methyl_path)
load(RNA_path)
# "BetaMatrix" "RNAMatrix" 
# check the shared samples
sample_shared <- intersect(substr(colnames(BetaMatrix),1,19),substr(colnames(RNAMatrix),1,19))
colnames(BetaMatrix) <- substr(colnames(BetaMatrix),1,19)
colnames(RNAMatrix) <- substr(colnames(RNAMatrix),1,19)
BetaMatrix <- BetaMatrix[,sample_shared]
RNAMatrix <- RNAMatrix[,sample_shared]

genes_all <- rownames(RNAMatrix)
genes_all <- gsub("\\|.*","", genes_all)
RNAMatrix <- RNAMatrix[-(1:29),]
rownames(RNAMatrix) <- genes_all[-(1:29)]

source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/CollapseCpGsGenes.R")

load("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/xCell.data.rda")
source("/sonas-hs/wigler/hpc/home/jsong/MethylDeConv/xCell_custom.R")
source("ESTIMATE_custom.R")

BetaMatrix_average <- CollapseCpGsGenes(BetaMatrix, include = "all")
BetaMatrix_maxvar <- CollapseCpGsGenes(BetaMatrix, method = "maxvar", include = "all")
BetaMatrix_pca <- CollapseCpGsGenes(BetaMatrix, method = "PCA", include = "all")

BetaMatrix_average_TSS200 <- CollapseCpGsGenes(BetaMatrix, include = "TSS200")
BetaMatrix_maxvar_TSS200 <- CollapseCpGsGenes(BetaMatrix, method = "maxvar", include = "TSS200")
BetaMatrix_pca_TSS200 <- CollapseCpGsGenes(BetaMatrix, method = "PCA", include = "TSS200")

BetaMatrix_average_TSS <- CollapseCpGsGenes(BetaMatrix, include = "TSS200&TSS1500")
BetaMatrix_maxvar_TSS <- CollapseCpGsGenes(BetaMatrix, method = "maxvar", include = "TSS200&TSS1500")
BetaMatrix_pca_TSS <- CollapseCpGsGenes(BetaMatrix, method = "PCA", include = "TSS200&TSS1500")

# xCell
tmp_rank <- apply(-BetaMatrix_average, 2, rank)
BetaMatrix_average_xCell <- xCellAnalysis(tmp_rank, rnaseq = TRUE)
print("xcell - methylation average - include all, completed!")

tmp_rank <- apply(-BetaMatrix_maxvar, 2, rank)
BetaMatrix_maxvar_xCell <- xCellAnalysis(tmp_rank, rnaseq = TRUE)
print("xcell - methylation maxvar - include all, completed!")

tmp_rank <- apply(-BetaMatrix_pca, 2, rank)
BetaMatrix_pca_xCell <- xCellAnalysis(tmp_rank, rnaseq = TRUE)
print("xcell - methylation pca - include all, completed!")

tmp_rank <- apply(-BetaMatrix_average_TSS200, 2, rank)
BetaMatrix_average_TSS200_xCell <- xCellAnalysis(tmp_rank, rnaseq = TRUE)
print("xcell - methylation average - include TSS200, completed!")

tmp_rank <- apply(-BetaMatrix_maxvar_TSS200, 2, rank)
BetaMatrix_maxvar_TSS200_xCell <- xCellAnalysis(tmp_rank, rnaseq = TRUE)
print("xcell - methylation maxvar - include TSS200, completed!")

tmp_rank <- apply(-BetaMatrix_pca_TSS200, 2, rank)
BetaMatrix_pca_TSS200_xCell <- xCellAnalysis(tmp_rank, rnaseq = TRUE)
print("xcell - methylation pca - include TSS200, completed!")

tmp_rank <- apply(-BetaMatrix_average_TSS, 2, rank)
BetaMatrix_average_TSS_xCell <- xCellAnalysis(tmp_rank, rnaseq = TRUE)
print("xcell - methylation average - include TSS200&TSS1500, completed!")

tmp_rank <- apply(-BetaMatrix_maxvar_TSS, 2, rank)
BetaMatrix_maxvar_TSS_xCell <- xCellAnalysis(tmp_rank, rnaseq = TRUE)
print("xcell - methylation maxvar - include TSS200&TSS1500, completed!")

tmp_rank <- apply(-BetaMatrix_pca_TSS, 2, rank)
BetaMatrix_pca_TSS_xCell <- xCellAnalysis(tmp_rank, rnaseq = TRUE)
print("xcell - methylation pca - include TSS200&TSS1500, completed!")

tmp_rank <- apply(RNAMatrix,2,rank)
RNAMatrix_xCell <- xCellAnalysis(tmp_rank, rnaseq = TRUE)
print("xcell - RNA, completed!")





# ESTIMATE

BetaMatrix_average_estimate <- estimate_custom(-BetaMatrix_average, platform = "illumina")
print("estimate - methylation average - include all, completed!")

BetaMatrix_maxvar_estimate <- estimate_custom(-BetaMatrix_maxvar, platform = "illumina")
print("estimate - methylation maxvar - include all, completed!")

BetaMatrix_pca_estimate <- estimate_custom(-BetaMatrix_pca, platform = "illumina")
print("estimate - methylation pca - include all, completed!")

BetaMatrix_average_TSS200_estimate <- estimate_custom(-BetaMatrix_average_TSS200, platform = "illumina")
print("estimate - methylation average - include TSS200, completed!")

BetaMatrix_maxvar_TSS200_estimate <- estimate_custom(-BetaMatrix_maxvar_TSS200, platform = "illumina")
print("estimate - methylation maxvar - include TSS200, completed!")

BetaMatrix_pca_TSS200_estimate <- estimate_custom(-BetaMatrix_pca_TSS200, platform = "illumina")
print("estimate  - methylation pca - include TSS200, completed!")

BetaMatrix_average_TSS_estimate <- estimate_custom(-BetaMatrix_average_TSS, platform = "illumina")
print("estimate  - methylation average - include TSS200&TSS1500, completed!")

BetaMatrix_maxvar_TSS_estimate <- estimate_custom(-BetaMatrix_maxvar_TSS, platform = "illumina")
print("estimate  - methylation maxvar - include TSS200&TSS1500, completed!")

BetaMatrix_pca_TSS_estimate <- estimate_custom(-BetaMatrix_pca_TSS, platform = "illumina")
print("estimate  - methylation pca - include TSS200&TSS1500, completed!")


RNAMatrix_estimate <- estimate_custom(RNAMatrix, platform = "illumina")
print("estimate  - RNA, completed!")

output_file_name <- paste0("TCGA", disease, "_xCell_ESTIMATE_collapse.RData")
save("BetaMatrix_average_xCell","BetaMatrix_maxvar_xCell","BetaMatrix_pca_xCell",
     "BetaMatrix_average_TSS200_xCell","BetaMatrix_maxvar_TSS200_xCell","BetaMatrix_pca_TSS200_xCell",
     "BetaMatrix_average_TSS_xCell", "BetaMatrix_maxvar_TSS_xCell", "BetaMatrix_pca_TSS_xCell", "RNAMatrix_xCell",
     "BetaMatrix_average_estimate", "BetaMatrix_maxvar_estimate", "BetaMatrix_pca_estimate",
     "BetaMatrix_average_TSS200_estimate", "BetaMatrix_maxvar_TSS200_estimate", "BetaMatrix_pca_TSS200_estimate",
     "BetaMatrix_average_TSS_estimate", "BetaMatrix_maxvar_TSS_estimate", "BetaMatrix_pca_TSS_estimate", "RNAMatrix_estimate",
     file = output_file_name)
