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


## SMAPE is undefined if any |actual|+|predicted| is 0.
## custom the function: if both actual and predicted are 0, the term is 0.


SMAPE_custom <- function(actual, predicted){
  id = which((actual == 0) & (predicted == 0) == TRUE)
  if (length(id) > 0){
    actual_noZero <- actual[-id]
    predicted_noZero <- predicted[-id]
  }
  else{
    actual_noZero <- actual
    predicted_noZero <- predicted
  }
  nom = ae(actual_noZero, predicted_noZero)
  denom = abs(actual_noZero) + abs(predicted_noZero)
  mean_tmp = sum(nom/denom)/length(actual)
  return(2*mean_tmp)
}


within_sample_SMAPE_custom <- function(true_proportions, deconv_res){
  SMAPE <- rep(NA, 12)
  for (i in 1:12){
    SMAPE[i] <- SMAPE_custom(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), 
                             as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]))
  }
  print(mean(SMAPE))
  return(mean(SMAPE))
}


within_sample_SMAPE_sd_custom <- function(true_proportions, deconv_res){
  SMAPE <- rep(NA, 12)
  for (i in 1:12){
    SMAPE[i] <- SMAPE_custom(as.numeric(true_proportions[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]), 
                             as.numeric(deconv_res[i,c("Bcell", "CD4T","CD8T","Mono","Neu","NK")]))
  }
  print(sd(SMAPE))
  return(sd(SMAPE))
}
