# probe selection (From compTable(mean methylation table) to ModelPars table)
source("pickCompProbes.R")
library(minfi)


#FlowSorted.Blood.450k.compTable  #FlowSorted.Blood.450k.JaffeModelPars
library(FlowSorted.Blood.450k)
dat <- pickCompProbes(preprocessQuantile(FlowSorted.Blood.450k), cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),
                      probeSelect = "both")


## DLPFC
library(FlowSorted.DLPFC.450k)
GRset_frontCortex <- preprocessQuantile(FlowSorted.DLPFC.450k)
##FlowSorted.DLPFC.450k.compTable, FlowSorted.DLPFC.450k.ModelPars
dat_frontCortex <- pickCompProbes(GRset_frontCortex, cellTypes = c("NeuN_pos", "NeuN_neg"), probeSelect = "both")
FlowSorted.DLPFC.450k.compTable <- dat_frontCortex[[1]]
save(FlowSorted.DLPFC.450k.compTable, file = "FlowSorted.DLPFC.450k.compTable.RData")
FlowSorted.DLPFC.450k.ModelPars <- dat_frontCortex[[1]][dat_frontCortex[[2]], 3:4]
save(FlowSorted.DLPFC.450k.ModelPars, file = "FlowSorted.DLPFC.450k.ModelPars.RData")

# library(FlowSorted.Blood.450k)
# GRset.normalized <- preprocessQuantile(FlowSorted.Blood.450k)
# GRset.normalized <- getBeta(GRset.normalized)
# mean_dat <- matrix(NA, dim(GRset.normalized)[1], 11)
# colnames(mean_dat) <- c("Fstat", "p.value", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "low", "high", "range")
# rownames(mean_dat) <- rownames(GRset.normalized)
# mean_dat[,3] <- apply(GRset.normalized[,which(substr(colnames(GRset.normalized),1,4) == "CD8+")],1,mean)
# mean_dat[,4] <- apply(GRset.normalized[,which(substr(colnames(GRset.normalized),1,4) == "CD4+")],1,mean)
# mean_dat[,5] <- apply(GRset.normalized[,which(substr(colnames(GRset.normalized),1,4) == "CD56")],1,mean)
# mean_dat[,6] <- apply(GRset.normalized[,which(substr(colnames(GRset.normalized),1,4) == "CD19")],1,mean)
# mean_dat[,7] <- apply(GRset.normalized[,which(substr(colnames(GRset.normalized),1,4) == "CD14")],1,mean)
# mean_dat[,8] <- apply(GRset.normalized[,which(substr(colnames(GRset.normalized),1,4) == "Gran")],1,mean)
# ### approximately close to FlowSorted.Blood.450k.compTable
# tmpMat <- GRset.normalized[,substr(colnames(GRset.normalized),1,4) %in% c("CD8+","CD4+", "CD56", "CD19", "CD14", "Gran")]
# group <- factor(substr(colnames(tmpMat),1,4)) 
# Fstat <- apply(tmpMat, 1, function(x){summary(aov(x~group))[[1]][1,4]})
# mean_dat[,1] <-Fstat
# pvalue <- apply(tmpMat, 1, function(x){summary(aov(x~group))[[1]][1,5]})
# mean_dat[,2] <-pvalue







  