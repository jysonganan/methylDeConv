# smoking EWAS example
setwd("/Users/junesong/Desktop/causal inference")
input_methyl <- read.table("GSE50660_matrix_processed.txt", header = T, sep = "\t")
library(EpiDISH)

print('Reference matrix set as whole blood by default!')
dat <- input_methyl[,-1]
rownames(dat) <- input_methyl[,1]
cell_Prop <- epidish(beta.m = dat, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
boxplot(cell_Prop)

write.table(cell_Prop, "GSE50660_cell_prop_smoking_Epidish.txt", sep = "\t", col.names = T, quote = F)

data(centDHSbloodDMC.m)
library(RefFreeEWAS)
cell_Prop_2 = projectMix(dat, centDHSbloodDMC.m)

common_probe=intersect(row.names(centBloodSub.m),row.names(dat))
cell_Prop_2 = projectMix(dat[common_probe,], centDHSbloodDMC.m[common_probe,])

## boxplot of the estimated cell proportions
par(mfrow=c(1,2))
boxplot(cell_Prop, main = "Epidish RPC")
boxplot(cell_Prop_2, main = "Houseman")

## 
plot(cell_Prop[,1],cell_Prop_2[,1])
abline(1,1, col = "red")

library(stats)
for(i in 1:7){
  tmp <- cbind(cell_Prop[,i], cell_Prop_2[,i])
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- c("epidish", "houseman")
  reg <- lm(epidish~houseman, data=tmp)
  plot(tmp, main = colnames(cell_Prop_2)[i])
  abline(reg, col = "blue")
  legend("topleft", bty="n", legend=paste("R2 is", format(summary(reg)$adj.r.squared)), col = "red", cex = 0.6)
}



for(i in 1:7){
  tmp <- cbind(cell_Prop[,i], cell_Prop_2[,i])
  colnames(tmp) <- c("epidishRPC", "houseman")
  tmp <- as.data.frame(tmp)
  reg <- lm(houseman~epidishRPC, data = tmp)
  plot(tmp, main = colnames(cell_Prop_2)[i])
  abline(reg, col = "blue")
  legend("topleft", bty="n", legend=paste("R2 is", format(summary(reg)$adj.r.squared)), col = "red", cex = 0.6)
}

