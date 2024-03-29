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



### Supplementary Figure 1
## barplots of deconvolution performance (average corr with true proportions) using OneVsAllttest with top 50, 100, 150, 200 probes

library(tidyr)
library(ggplot2)
cor_Houseman <- c(0.971, 0.978, 0.964, 0.964)
cor_RPC <- c(0.973, 0.978, 0.969, 0.964)
cor_CBS <- c(0.952, 0.956, 0.966, 0.968)
cor_MethylResolver <- c(0.971, 0.978, 0.964, 0.964)

df1 <- data.frame(FeatureSelection = c("top 50","top 100","top 150","top 200"),
                  Houseman = cor_Houseman, RPC = cor_RPC, CBS = cor_CBS, MethylResolver = cor_MethylResolver)

print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
             aes(x = factor(FeatureSelection, level = c("top 50","top 100","top 150","top 200")), y = Spearman_Correlation, fill = Deconvolution)) +
        geom_bar(stat = 'identity', position = 'dodge')+
        geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
        labs(x = "Feature selection with oneVsAllttest")+
        labs(y = "Average Spearman correlation")+
        labs(fill = "Deconvolution algorithms"))


cor_Houseman <- c(0.9710931, 0.9783402, 0.9640545, 0.9640545)
cor_RPC <- c(0.9734393, 0.9783402, 0.9688164, 0.9640545)
cor_CBS <- c(0.9515192, 0.9562116, 0.9663312, 0.9682206)
cor_MethylResolver <- c(0.9711626, 0.9783402, 0.9640545, 0.9640545)
df1 <- data.frame(FeatureSelection = c("top 50","top 100","top 150","top 200"),
                  Houseman = cor_Houseman, RPC = cor_RPC, CBS = cor_CBS, MethylResolver = cor_MethylResolver)
ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
       aes(x = factor(FeatureSelection, level = c("top 50","top 100","top 150","top 200")),
           y = Spearman_Correlation, group = Deconvolution)) +
  geom_point(aes(col = Deconvolution))+
  geom_line(aes(colour = Deconvolution, alpha= Deconvolution, linetype= Deconvolution, size = Deconvolution))+
  #geom_line(aes(linetype= Deconvolution, colour = Deconvolution), size = 1)+
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "solid"))+
  scale_alpha_manual(values=c(1,0.4,0.6,1))+
  scale_color_manual(values=c('violet','green','red', 'cyan3'))+
  scale_size_manual(values=c(1,1.4,1.4,1))+
  labs(x = "oneVsAllttest")+
  labs(y = "Average Spearman correlation")+
  labs(group = "Deconvolution algorithms") +
  theme_grey(base_size = 15)



