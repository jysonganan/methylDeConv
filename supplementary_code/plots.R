## plots

### Figure 1
## barplots of deconvolution performance (average corr with true proportions) using different feature selection + deconvolution algorithms
library(tidyr)
library(ggplot2)
cor_Houseman <- c(0.978, 0.969, 0.950, 0.966, 0.947, 0.947, 0.831, 0.955)
cor_RPC <- c(0.978, 0.978, 0.828, 0.983, 0.952, 0.938, 0.800, 0.950)
cor_CBS <- c(0.956, 0.962, 0.839, 0.952, 0.942, 0.950, 0.846, 0.932)
cor_MethylResolver <- c(0.978, 0.983, 0.868, 0.978, 0.956, 0.920, 0.788, 0.959)

#####  or load results from featureSelection_Comparison.R

# cor_Houseman <- c(cor_Houseman_1, cor_Houseman_2, cor_Houseman_3, cor_Houseman_4, cor_Houseman_5,
#                   cor_Houseman_6, cor_Houseman_7, cor_Houseman_8, cor_Houseman_9)
# cor_RPC <- c(cor_RPC_1, cor_RPC_2, cor_RPC_3, cor_RPC_4, cor_RPC_5, cor_RPC_6, cor_RPC_7, cor_RPC_8, cor_RPC_9)
# cor_CBS <- c(cor_CBS_1, cor_CBS_2, cor_CBS_3, cor_CBS_4, cor_CBS_5, cor_CBS_6, cor_CBS_7, cor_CBS_8, cor_CBS_9)
# cor_MethylResolver <- c(cor_MethylResolver_1, cor_MethylResolver_2, cor_MethylResolver_3, cor_MethylResolver_4,
#                         cor_MethylResolver_5, cor_MethylResolver_6, cor_MethylResolver_7, cor_MethylResolver_8,
#                         cor_MethylResolver_9)

df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                       "glmnetpreselect", "RFpreselect", "Rfepreselect"),
                  Houseman = cor_Houseman, RPC = cor_RPC, CBS = cor_CBS, MethylResolver = cor_MethylResolver)

print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
       aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                                  "glmnetpreselect", "RFpreselect", "Rfepreselect")), y = Spearman_Correlation, fill = Deconvolution)) +
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  labs(x = "Feature selection")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Deconvolution algorithms"))
