
##### analyze Collaspe on Benchmark data (GEO with FACS)
corr <- rep(NA,18)
for (i in 1:18){
  corr[i] <- cor(res_methyl_houseman[i,], as.numeric(facs_prop[i,]),method = "spearman")
}

corr <- rep(NA,18)
for (i in 1:18){
  corr[i] <- cor(res_methyl_cbs[i,], as.numeric(facs_prop[i,]),method = "spearman")
}


corr <- rep(NA,18)
for (i in 1:18){
  corr[i] <- cor(res_methyl_rpc[i,], as.numeric(facs_prop[i,]),method = "spearman")
}




corr <- rep(NA,6)
for (i in 1:6){
  corr[i] <- cor(res_methyl_houseman[,i], as.numeric(facs_prop[,i]),method = "spearman")
}


corr <- rep(NA,6)
for (i in 1:6){
  corr[i] <- cor(res_methyl_cbs[,i], as.numeric(facs_prop[,i]),method = "spearman")
}


corr <- rep(NA,6)
for (i in 1:6){
  corr[i] <- cor(res_methyl_rpc[,i], as.numeric(facs_prop[,i]),method = "spearman")
}

res_collapse <- cbind(res_collapse_average_cibersort[,4],
                 res_collapse_average_cibersort[,5]+res_collapse_average_cibersort[,6]+res_collapse_average_cibersort[,7],
                 res_collapse_average_cibersort[,11]+res_collapse_average_cibersort[,12],
                 res_collapse_average_cibersort[,1]+res_collapse_average_cibersort[,2],
                 res_collapse_average_cibersort[,13],
                 res_collapse_average_cibersort[,21]+res_collapse_average_cibersort[,22])



res_collapse <- cbind(res_collapse_maxvar_cibersort[,4],
                 res_collapse_maxvar_cibersort[,5]+res_collapse_maxvar_cibersort[,6]+res_collapse_maxvar_cibersort[,7],
                 res_collapse_maxvar_cibersort[,11]+res_collapse_maxvar_cibersort[,12],
                 res_collapse_maxvar_cibersort[,1]+res_collapse_maxvar_cibersort[,2],
                 res_collapse_maxvar_cibersort[,13],
                 res_collapse_maxvar_cibersort[,21]+res_collapse_maxvar_cibersort[,22])




res_collapse <- cbind(res_collapse_pca_cibersort[,4],
                 res_collapse_pca_cibersort[,5]+res_collapse_pca_cibersort[,6]+res_collapse_pca_cibersort[,7],
                 res_collapse_pca_cibersort[,11]+res_collapse_pca_cibersort[,12],
                 res_collapse_pca_cibersort[,1]+res_collapse_pca_cibersort[,2],
                 res_collapse_pca_cibersort[,13],
                 res_collapse_pca_cibersort[,21]+res_collapse_pca_cibersort[,22])


res_collapse <- cbind(res_collapse_average_nearest_cibersort[,4],
                 res_collapse_average_nearest_cibersort[,5]+res_collapse_average_nearest_cibersort[,6]+res_collapse_average_nearest_cibersort[,7],
                 res_collapse_average_nearest_cibersort[,11]+res_collapse_average_nearest_cibersort[,12],
                 res_collapse_average_nearest_cibersort[,1]+res_collapse_average_nearest_cibersort[,2],
                 res_collapse_average_nearest_cibersort[,13],
                 res_collapse_average_nearest_cibersort[,21]+res_collapse_average_nearest_cibersort[,22])

res_collapse <- cbind(res_collapse_average_nearest_TSS_cibersort[,4],
                 res_collapse_average_nearest_TSS_cibersort[,5]+res_collapse_average_nearest_TSS_cibersort[,6]+res_collapse_average_nearest_TSS_cibersort[,7],
                 res_collapse_average_nearest_TSS_cibersort[,11]+res_collapse_average_nearest_TSS_cibersort[,12],
                 res_collapse_average_nearest_TSS_cibersort[,1]+res_collapse_average_nearest_TSS_cibersort[,2],
                 res_collapse_average_nearest_TSS_cibersort[,13],
                 res_collapse_average_nearest_TSS_cibersort[,21]+res_collapse_average_nearest_TSS_cibersort[,22])



res_collapse <- cbind(res_collapse_maxvar_nearest_cibersort[,4],
                 res_collapse_maxvar_nearest_cibersort[,5]+res_collapse_maxvar_nearest_cibersort[,6]+res_collapse_maxvar_nearest_cibersort[,7],
                 res_collapse_maxvar_nearest_cibersort[,11]+res_collapse_maxvar_nearest_cibersort[,12],
                 res_collapse_maxvar_nearest_cibersort[,1]+res_collapse_maxvar_nearest_cibersort[,2],
                 res_collapse_maxvar_nearest_cibersort[,13],
                 res_collapse_maxvar_nearest_cibersort[,21]+res_collapse_maxvar_nearest_cibersort[,22])

res_collapse <- cbind(res_collapse_maxvar_nearest_TSS_cibersort[,4],
                 res_collapse_maxvar_nearest_TSS_cibersort[,5]+res_collapse_maxvar_nearest_TSS_cibersort[,6]+res_collapse_maxvar_nearest_TSS_cibersort[,7],
                 res_collapse_maxvar_nearest_TSS_cibersort[,11]+res_collapse_maxvar_nearest_TSS_cibersort[,12],
                 res_collapse_maxvar_nearest_TSS_cibersort[,1]+res_collapse_maxvar_nearest_TSS_cibersort[,2],
                 res_collapse_maxvar_nearest_TSS_cibersort[,13],
                 res_collapse_maxvar_nearest_TSS_cibersort[,21]+res_collapse_maxvar_nearest_TSS_cibersort[,22])




res_collapse_normalized <- apply(res_collapse,1,function(x){return(x/sum(x))})
res_collapse_normalized <- t(res_collapse_normalized)

corr <- rep(NA,18)
for (i in 1:18){
  corr[i] <- cor(res_collapse_normalized[i,], as.numeric(facs_prop[i,]),method = "spearman")
}


corr <- rep(NA,6)
for (i in 1:6){
  corr[i] <- cor(res_collapse_normalized[,i], as.numeric(facs_prop[,i]),method = "spearman")
}

