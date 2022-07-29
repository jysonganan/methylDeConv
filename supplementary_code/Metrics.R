########################
## Spearman correlations
########################

within_sample_corr <- function(true_proportions, deconv_res){
  corr <- rep(NA, 12)
  for (i in 1:12){
    corr[i] <- cor(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]),
                   as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), method = "spearman")
  }
  print(mean(corr))
  return(mean(corr))
}


within_sample_corr_sd <- function(true_proportions, deconv_res){
  corr <- rep(NA, 12)
  for (i in 1:12){
    corr[i] <- cor(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]),
                   as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), method = "spearman")
  }
  print(sd(corr))
  return(sd(corr))
}




############
## RMSE
############
library(Metrics)
within_sample_RMSE <- function(true_proportions, deconv_res){
  RMSE <- rep(NA, 12)
  for (i in 1:12){
    RMSE[i] <- rmse(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), 
                    as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]))
  }
  print(mean(RMSE))
  return(mean(RMSE))
}


within_sample_RMSE_sd <- function(true_proportions, deconv_res){
  RMSE <- rep(NA, 12)
  for (i in 1:12){
    RMSE[i] <- rmse(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), 
                    as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]))
  }
  print(sd(RMSE))
  return(sd(RMSE))
}




############
## SMAPE
############
# https://medium.com/@davide.sarra/how-to-interpret-smape-just-like-mape-bf799ba03bdc

within_sample_SMAPE <- function(true_proportions, deconv_res){
  SMAPE <- rep(NA, 12)
  for (i in 1:12){
    SMAPE[i] <- smape(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), 
                      as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]))
  }
  print(mean(SMAPE))
  return(mean(SMAPE))
}


within_sample_SMAPE_sd <- function(true_proportions, deconv_res){
  SMAPE <- rep(NA, 12)
  for (i in 1:12){
    SMAPE[i] <- smape(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), 
                      as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]))
  }
  print(sd(SMAPE))
  return(sd(SMAPE))
}


