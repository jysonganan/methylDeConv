# venn plot of the shared probes of different reference profiles for blood

load("FlowSorted.Blood.EPIC.IDOLModelPars.RData")
EPIC <- rownames(FlowSorted.Blood.EPIC.IDOLModelPars)

library(EpiDISH)
data("centDHSbloodDMC.m")
EpiDISH <- rownames(centDHSbloodDMC.m)

library(FlowSorted.Blood.450k)
data("FlowSorted.Blood.450k.JaffeModelPars")
Jaffe450k <- rownames(FlowSorted.Blood.450k.JaffeModelPars)

library(VennDiagram)
grid.newpage()
draw.triple.venn(area1 = 450, area2 = 333, area3 = 600, n12 = 13, n23 = 110, n13 = 54, 
                 n123 = 7, category = c("EPIC", "EpiDISH", "Jaffe450k"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))

#cbind(centDHSbloodDMC.m[intersect(EpiDISH,Jaffe450k),], FlowSorted.Blood.450k.JaffeModelPars[intersect(EpiDISH,Jaffe450k),])[1:10,]


###
load("/Users/junesong/Downloads/immuneCelltypeGenevsMethylCellProp_KICH.RData")
plot(res_collapse_average_Nearest_TSS_cibersort[,"T cells CD8"],log10(RNAMatrix["CD8A",]+1))
cor(res_collapse_average_Nearest_TSS_cibersort[,"T cells CD8"],log10(RNAMatrix["CD8A",]+1))

load("/Users/junesong/Downloads/immuneCelltypeGenevsMethylCellProp_BRCA.RData")
### best one!!
plot(res_collapse_average_Shore_cibersort[,"T cells CD8"],log10(RNAMatrix["CD8A",]+1))
cor(res_collapse_average_Shore_cibersort[,"T cells CD8"],log10(RNAMatrix["CD8A",]+1))


plot(res_collapse_average_Nearest_TSS_cibersort[,"B cells naive"]+res_collapse_average_Nearest_TSS_cibersort[,"B cells memory"],
     log10(RNAMatrix["CD19",]+1))
cor(res_collapse_average_Nearest_TSS_cibersort[,"B cells naive"]+res_collapse_average_Nearest_TSS_cibersort[,"B cells memory"],
    log10(RNAMatrix["CD19",]+1))


plot(res_collapse_average_Shore_cibersort[,"Monocytes"],
     log10(RNAMatrix["CD14",]+1))
cor(res_collapse_average_Shore_cibersort[,"Monocytes"],
    log10(RNAMatrix["CD14",]+1))


####################

load("/Users/junesong/Downloads/KICHforGP.RData")
x1 <- BetaMatrix_average[,1]
y1 <- RNAMatrix[rownames(BetaMatrix_average),1]
x2 <- BetaMatrix_average[,40]
y2 <- RNAMatrix[rownames(BetaMatrix_average),1]
dat1 <- cbind(x1,y1)
dat2 <- cbind(x2,y2)
colnames(dat1)[1] = colnames(dat2)[1] <- "methylCollaspe"
colnames(dat1)[2] = colnames(dat2)[2] <- "rna"
write.csv(dat1, "/Users/junesong/Downloads/dat1.csv")
write.csv(dat2, "/Users/junesong/Downloads/dat2.csv")

write.csv(BetaMatrix_average,"/Users/junesong/Downloads/BetaMatrix_average.csv")
write.csv(RNAMatrix[rownames(BetaMatrix_average),], "/Users/junesong/Downloads/RNAMatrix_average.csv")


### heatmap

heatmap()






corr <- rep(NA,66)
for (i in 1:66){
  corr[i] <- cor(BetaMatrix_average[,i],log2(RNAMatrix[rownames(BetaMatrix_average),i]+1))
}
#-0.3934161
corr <- rep(NA,66)
for (i in 1:66){
  corr[i] <- cor(BetaMatrix_average[,i],RNAMatrix[rownames(BetaMatrix_average),i])
}
# -0.229


corr <- rep(NA,66)
for (i in 1:66){
  corr[i] <- cor(BetaMatrix_average_nearest[,i],log2(RNAMatrix[rownames(BetaMatrix_average_nearest),i]+1))
}
#-0.3460073

corr <- rep(NA,66)
for (i in 1:66){
  corr[i] <- cor(BetaMatrix_average_nearest_TSS[,i],log2(RNAMatrix[rownames(BetaMatrix_average_nearest_TSS),i]+1))
}
#-0.3614898


plot(BetaMatrix_average[,1],log2(RNAMatrix[rownames(BetaMatrix_average),1]+1))
plot(BetaMatrix_average_nearest[,1],log2(RNAMatrix[rownames(BetaMatrix_average_nearest),1]+1))


rna <- log10(RNAMatrix[rownames(BetaMatrix_average_nearest),1]+1)
beta <- BetaMatrix_average_nearest[,1]
order <- order(rna)
rna <- rna[order]
beta <- beta[order]
matplot(1:455,cbind(rna,beta),type="l",col=c("blue","red"))





## signature gene vs cell type proportions
