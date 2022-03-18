## plots
#####################
### Figure 1
#####################
## barplots of deconvolution performance (average corr with true proportions) using different feature selection + deconvolution algorithms
library(tidyr)
library(ggplot2)
cor_Houseman <- c(0.978, 0.969, 0.950, 0.966, 0.947, 0.947, 0.831, 0.955)
cor_RPC <- c(0.978, 0.978, 0.828, 0.983, 0.952, 0.938, 0.800, 0.950)
cor_CBS <- c(0.956, 0.962, 0.839, 0.952, 0.942, 0.950, 0.846, 0.932)
cor_MethylResolver <- c(0.978, 0.983, 0.868, 0.978, 0.956, 0.920, 0.788, 0.959)
#### or load results from featureSelection_Comparison.R
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
  labs(fill = "Deconvolution algorithms")+
  coord_cartesian(ylim=c(0.5,1)))



#####################
### Figure 2
#####################
Houseman <- c(0.978, 0.947, 0.95, 0.94)
RPC <- c(0.978, 0.938, 0.94, 0.94)
CBS <- c(0.956, 0.950, 0.95, 0.92)
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","glmnetpreselect", "oneVsAllttest","glmnetpreselect"), 
                  facet = c("A. EPIC reference library", "A. EPIC reference library", "B. EPIC + Epithelial reference library", "B. EPIC + Epithelial reference library"), 
                  Houseman = Houseman, RPC = RPC, CBS = CBS)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(FeatureSelection,facet)), 
       aes(x = factor(FeatureSelection, level = c("oneVsAllttest","glmnetpreselect")), y = Spearman_Correlation, fill = Deconvolution)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  facet_wrap(~factor(facet,level = c("A. EPIC reference library","B. EPIC + Epithelial reference library")))+
  labs(x = "Feature selection")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Deconvolution algorithms")+
  coord_cartesian(ylim=c(0.5,1))+
  theme(strip.text.x = element_text(size = 15))




#####################
### Figure 6
#####################
Houseman <- c(0.989, 0.970, 0.901, 0.639,0.258, 0.990, 0.983,0.973,0.894,0.711)
RPC <- c(0.991,0.970,0.828,0.542,0.236,0.992,0.986,0.982,0.939,0.829)
CBS <- c(0.974, 0.955,0.875,0.643,0.324,0.980,0.981,0.978,0.919,0.772)
MethylResolver <- c(0.992,0.974,0.864,0.651,0.368,0.991,0.985,0.980,0.939,0.819)
df1 <- data.frame(NonImmuneProportions = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9", "0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9"),
                  facet = c(rep("EPIC reference library", 5), rep("EPIC + Epithelial reference library",5)),
                  Houseman = Houseman, RPC = RPC, CBS = CBS, MethylResolver = MethylResolver)
ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(NonImmuneProportions,facet)),
       aes(x = factor(NonImmuneProportions, level = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9")), 
           y = Spearman_Correlation, group =Deconvolution)) +
  geom_point(aes(col = Deconvolution))+ 
  geom_line(aes(colour = Deconvolution), size = 1)+
  #geom_line(group = 4, aes(col = factor(Deconvolution)))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")))+
  labs(x = "Non-immune (Epithelial) proportions")+
  labs(y = "Average Spearman correlation")+
  labs(group= "Deconvolution algorithms")

se = c(0.00116, 0.00185, 0.00404, 0.01235, 0.01702, 0.00111, 0.00156, 0.00158, 0.00557, 0.01085,  
  0.00098, 0.00189, 0.00736, 0.01340, 0.01790, 0.00098, 0.00133, 0.00121, 0.00351, 0.00707,
  0.00201, 0.00255, 0.00540, 0.01171, 0.01696, 0.00183, 0.00146, 0.00151, 0.00472, 0.01010,
  0.00093, 0.00174, 0.00560, 0.01210, 0.01675, 0.00107, 0.00150, 0.00132, 0.00353, 0.00772)
#se <- se * sqrt(600)
se <- se * 1.96
ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(NonImmuneProportions,facet)),
       aes(x = factor(NonImmuneProportions, level = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9")), 
           y = Spearman_Correlation, group =Deconvolution)) +
  geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se, col = Deconvolution), width=.1) +
  geom_point(aes(col = Deconvolution))+ 
  geom_line(aes(col = Deconvolution), size = 1)+
  #geom_line(group = 4, aes(col = factor(Deconvolution)))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")))+
  labs(x = "Non-immune (Epithelial) proportions")+
  labs(y = "Average Spearman correlation")+
  labs(group= "Deconvolution algorithms")





#####################
### Figure 7
#####################
Houseman <- c(0.707,0.753,0.662,0.406,0.161,0.538,0.561,0.392,0.066,-0.016)
RPC <- c(0.926,0.951,0.840,0.528,0.191,0.940,0.927,0.922,0.825,0.601)
CBS <- c(0.936,0.947,0.812,0.504,0.229,0.938,0.933,0.925,0.819,0.586)
MethylResolver <- c(0.923,0.948,0.874,0.606,0.266,0.938,0.932,0.924,0.836,0.600)
df1 <- data.frame(NonImmuneProportions = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9", "0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9"),
                  facet = c(rep("EPIC reference library", 5), rep("EPIC + Epithelial reference library",5)),
                  Houseman = Houseman, RPC = RPC, CBS = CBS, MethylResolver = MethylResolver)

se = c(0.00884,0.00819,0.01214,0.01586,0.01673,	0.01382,0.01143,0.01245,0.01676,0.01761,
       0.00324,0.00258,0.00706,0.01374,0.01759, 0.00282,0.00382,0.00293,0.00612,0.01403,
       0.00297,0.00242,0.00905,0.01483,0.01764,	0.00299,0.00356,0.00304,0.00754,0.01440,
       0.00330,0.00288,0.00612,0.01297,0.01739,0.00286,0.00373,0.00299,0.00702,0.01409)

se <- se * 1.96
ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(NonImmuneProportions,facet)),
       aes(x = factor(NonImmuneProportions, level = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9")), 
           y = Spearman_Correlation, group =Deconvolution)) +
  geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se, col = Deconvolution), width=.1) +
  geom_point(aes(col = Deconvolution))+ 
  geom_line(aes(col = Deconvolution), size = 1)+
  #geom_line(group = 4, aes(col = factor(Deconvolution)))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")))+
  labs(x = "Non-immune (Epithelial) proportions")+
  labs(y = "Average Spearman correlation")+
  labs(group= "Deconvolution algorithms")






#####################
### Figure 8
######################## 
##### within cell type MethylResolver
Bcell <- c(0.999,0.992, 0.982, 0.943, 0.727, 0.999,0.997,0.997,0.993,0.930)
CD4T <- c(0.996,0.991,0.933,0.576,0.406,0.996,0.992,0.981,0.941,0.873)
CD8T <- c(0.993,0.972,0.903,0.585,0.207,0.993,0.989,0.970,0.949,0.833)
Mono <- c(0.997,0.991,0.966,0.827,0.570,0.997,0.997,0.994,0.968,0.920)
Neu <- c(0.998,0.994,0.966,0.851,0.544,0.998,0.998,0.994,0.978,0.773)
NK <- c(0.997,0.991,0.961,0.901,0.735,0.997,0.996,0.995,0.981,0.946)
Epithelial <- c(NA,NA,NA,NA,NA,NA,0.996,0.997,0.996,0.973)
df1 <- data.frame(NonImmuneProportions = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9", "0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9"),
                  facet = c(rep("EPIC reference library", 5), rep("EPIC + Epithelial reference library",5)),
                  Bcell = Bcell, CD4T = CD4T, CD8T = CD8T, Mono = Mono, Neu = Neu, NK = NK, Nonimmune_Epithelial = Epithelial)
ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(NonImmuneProportions,facet)),
       aes(x = factor(NonImmuneProportions, level = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9")), y = Spearman_Correlation, fill = Deconvolution)) +
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")),dir = "v")+
  labs(x = "Non-immune (Epithelial) proportions")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Cell Types")





#####################
### Figure 9
#####################
### within cell type MethylResolver
Bcell <- c(0.996,0.990,0.982,0.937,0.704,0.996,0.990,0.991,0.981,0.824)
CD4T <- c(0.979,0.983,0.930,0.591,0.480,0.981,0.987,0.951,0.864,0.725)
CD8T <- c(0.979,0.963,0.880,0.607,0.326,0.983,0.979,0.943,0.876,0.654)
Mono <- c(0.989,0.985,0.950,0.776,0.479,0.991,0.993,0.985,0.925,0.822)
Neu <- c(0.989,0.989,0.946,0.854,0.525,0.989,0.988,0.975,0.961,0.632)
NK <- c(0.976,0.984,0.959,0.863,0.626,0.972,0.986,0.980,0.940,0.903)
Epithelial <- c(NA,NA,NA,NA,NA,NA,0.979,0.992,0.988,0.919)
df1 <- data.frame(NonImmuneProportions = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9", "0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9"),
                  facet = c(rep("EPIC reference library", 5), rep("EPIC + Epithelial reference library",5)),
                  Bcell = Bcell, CD4T = CD4T, CD8T = CD8T, Mono = Mono, Neu = Neu, NK = NK, Nonimmune_Epithelial = Epithelial)
ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(NonImmuneProportions,facet)),
       aes(x = factor(NonImmuneProportions, level = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9")), y = Spearman_Correlation, fill = Deconvolution)) +
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")),dir = "v")+
  labs(x = "Non-immune (Epithelial) proportions")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Cell Types")










#####################
### Figure 11, 12
#####################

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


###################
## Figure 13
##################

load("res_140038.RData")
library(ggplot2)
library(tidyr)
pdf(file = "Plot2.pdf",  
    width = 10,
    height = 8) 
df <- gather(gender_dat_blood, series,value,-group)
ggplot(df) + geom_boxplot(aes(series ,value,color=group)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE140038-RPC-onevsAllttest")
dev.off()


###################
## Figure 14
###################

load("res_111631_111165.RData")
library(ggplot2)
library(tidyr)
pdf(file = "Plot3.pdf",  
    width = 10,
    height = 8) 
df <- gather(gender_dat_blood, series,value,-group)
ggplot(df) + geom_boxplot(aes(series ,value,color=group)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE111631+GSE111165-RPC-onevsAllttest")
dev.off()




##################
## Figure 16
##################
load("res_112308.RData")
pdf(file = "Plot4.pdf",  
    width = 10,
    height = 8) 
df <- gather(gender_dat_blood, series,value,-group)
ggplot(df) + geom_boxplot(aes(series ,value,color=group)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE112308-Melanoma_Houseman-RPC-onevsAllttest")
dev.off()



#####################
### Supplementary Figure 1
#####################

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
