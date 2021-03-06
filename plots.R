## barplots of deconvolution performance (average corr with true proportions) using different feature selection + Houseman/RPC/CBS
library(tidyr)
library(ggplot2)
Houseman <- c(0.978, 0.969, 0.950, 0.966, 0.947, 0.947, 0.831, 0.866)
RPC <- c(0.978, 0.978, 0.828, 0.983, 0.952, 0.938, 0.800, 0.886)
CBS <- c(0.956, 0.962, 0.839, 0.952, 0.942, 0.940, 0.850, 0.877)
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                       "glmnetpreselect", "RFpreselect", "OptVariables"), Houseman = Houseman, RPC = RPC, CBS = CBS)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection), 
       aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                                  "glmnetpreselect", "RFpreselect", "OptVariables")), y = Spearman_Correlation, fill = Deconvolution)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  labs(x = "Feature selection")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Deconvolution algorithms")






### matched 450k, EPIC samples of blood and saliva. Facet: EPIC / EPIC + Epithelial library
Houseman <- c(0.895, 0.838, 0.819, 0.829)
RPC <- c(0.819, 0.838, 0.752, 0.885)
CBS <- c(0.810, 0.848, 0.638, 0.867)
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","glmnetpreselect", "oneVsAllttest","glmnetpreselect"), 
                  facet = c("EPIC reference library", "EPIC reference library", "EPIC + Epithelial reference library", "EPIC + Epithelial reference library"), 
                  Houseman = Houseman, RPC = RPC, CBS = CBS)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(FeatureSelection,facet)), 
       aes(x = factor(FeatureSelection, level = c("oneVsAllttest","glmnetpreselect")), y = Spearman_Correlation, fill = Deconvolution)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")))+
  labs(x = "Feature selection")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Deconvolution algorithms")



Houseman <- c(0.619, 0.931, 0.410, 0.668)
RPC <- c(0.855, 0.962, 0.329, 0.683)
CBS <- c(0.827, 0.911, 0.372, 0.614)
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","glmnetpreselect", "oneVsAllttest","glmnetpreselect"), 
                  facet = c("EPIC reference library", "EPIC reference library", "EPIC + Epithelial reference library", "EPIC + Epithelial reference library"), 
                  Houseman = Houseman, RPC = RPC, CBS = CBS)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(FeatureSelection,facet)), 
       aes(x = factor(FeatureSelection, level = c("oneVsAllttest","glmnetpreselect")), y = Spearman_Correlation, fill = Deconvolution)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")))+
  labs(x = "Feature selection")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Deconvolution algorithms")







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


<<<<<<< HEAD
=======
load("/Users/junesong/Desktop/causal inference/Flow450kProbe_sub.RData")
length(probes_multiclassGlmnet_sub)
length(probes_pairwiseGlmnet_sub)
load("/Users/junesong/Desktop/causal inference/Flow450kProbeCompTable_1.RData")

compTable_1 <- compTable
probes_oneVsAllttest_1 <- probes_oneVsAllttest
probes_multiclassGlmnet_1 <- probes_multiclassGlmnet
probes_pairwiseGlmnet_1 <- probes_pairwiseGlmnet
probes_oneVsAllLimma_1 <- probes_oneVsAllLimma
probes_pairwiseLimma_1 <- probes_pairwiseLimma


load("/Users/junesong/Desktop/causal inference/Flow450kProbeCompTable.RData")
#1
length(probes_oneVsAllttest_1)
length(probes_oneVsAllLimma_1)
length(probes_pairwiseLimma_1)
length(probes_pairwiseGlmnet_1)
length(probes_multiclassGlmnet_1)
#2
length(probes_oneVsAllttest)
length(probes_oneVsAllLimma)
length(probes_pairwiseLimma)
length(probes_pairwiseGlmnet)
length(probes_multiclassGlmnet)
#3
length(probes_pairwiseGlmnet_sub)
length(probes_multiclassGlmnet_sub)

## compare sub with full
length(intersect(probes_pairwiseGlmnet_sub,probes_pairwiseGlmnet))
length(intersect(probes_multiclassGlmnet_sub,probes_multiclassGlmnet))

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(area1 = 9047, area2 = 7987,  cross.area = 6986, 
                   category = c("EN on all probes", "EN on 80% probes"), 
                   lty = "blank", fill = c("blue", "mediumorchid")) # area2 7916 area3 7099

grid.newpage()
draw.pairwise.venn(area1 = 910, area2 = 817,  cross.area = 709, 
                   category = c("EN on all probes", "EN on 80% probes"), 
                   lty = "blank", fill = c("yellow", "red"))  #area2 829 area3 = 709


length(intersect(probes_oneVsAllttest, probes_oneVsAllttest_1)) ## the two are same
length(intersect(probes_oneVsAllLimma, probes_oneVsAllLimma_1)) ## the two are same
length(intersect(probes_pairwiseLimma, probes_pairwiseLimma_1)) ## the two are same
length(intersect(probes_pairwiseGlmnet,probes_pairwiseGlmnet_1))
length(intersect(probes_multiclassGlmnet,probes_multiclassGlmnet_1))

grid.newpage()
draw.pairwise.venn(area1 = 9047, area2 = 9059,  cross.area = 9047, 
                   category = c("1", "2"), 
                   lty = "blank", fill = c("blue", "mediumorchid"))

grid.newpage()
draw.pairwise.venn(area1 = 910, area2 = 968,  cross.area = 910, 
                   category = c("1", "2"), 
                   lty = "blank", fill = c("yellow", "red"))


length(probes_oneVsAllttest)
length(probes_oneVsAllLimma)
length(probes_pairwiseLimma)
length(probes_pairwiseGlmnet)
length(probes_multiclassGlmnet)
grid.newpage()
draw.quintuple.venn(area1 = 600, area2= 600, area3 = 681, area4 = 9047, area5 = 910, n12 = 538, n13 = 112, n14 = 314,
                    n15 = 145, n23 = 116, n24 = 318, n25 = 153, n34 = 618, n35 = 165, n45 = 546, n123 = 111, n124 = 298,
                    n125 = 145, n134 = 99, n135 = 67, n145 = 128, n234 = 103, n235 = 68, n245 = 130, n345 = 133, n1234 = 98, n1235 = 67,
                    n1245 = 128, n1345 = 62, n2345 = 63, n12345 = 62, category = c("1", "2", "3", "4","5"), 
                    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                    cat.cex = 2,
                    margin = 0.05,
                    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
                            1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
                    ind = TRUE)

grid.newpage()
draw.triple.venn(area1 = 600, area2= 600, area3 = 681, n12 = 538, n23 = 116, n13 = 112, 
                 n123 = 111, category = c("1", "2", "3"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))



>>>>>>> 940b0d2637da5d15910cca3c2c439f9158f24513
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
