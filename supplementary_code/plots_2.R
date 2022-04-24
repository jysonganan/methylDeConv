## plots
#####################
### Figure 1
#####################
## barplots of deconvolution performance (average corr with true proportions) using different feature selection + deconvolution algorithms
library(tidyr)
library(ggplot2)

pdf(file = "Figure1.pdf",  
    width = 16,
    height = 8) 

cor_Houseman <- c(0.978, 0.969, 0.950, 0.966, 0.947, 0.947, 0.831, 0.964)
cor_RPC <- c(0.978, 0.978, 0.828, 0.983, 0.952, 0.938, 0.800, 0.959)
cor_CBS <- c(0.956, 0.962, 0.839, 0.952, 0.942, 0.950, 0.846, 0.952)
cor_MethylResolver <- c(0.978, 0.983, 0.868, 0.978, 0.956, 0.920, 0.788, 0.955)

df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                       "glmnetpreselect", "RFpreselect", "Rfepreselect"),
                  Houseman = cor_Houseman, RPC = cor_RPC, CBS = cor_CBS, MethylResolver = cor_MethylResolver)


sd_1 <- c(0.0299, 0.0521, 0.0530, 0.0548, 0.0662, 0.0658, 0.1200, 0.0516,
          0.0299, 0.0299, 0.1766, 0.0252, 0.0679, 0.0743, 0.1670, 0.0659,
          0.0639, 0.0471, 0.1831, 0.0501, 0.0812, 0.0470, 0.1028, 0.0501,
          0.0299, 0.0282, 0.1804, 0.0299, 0.0692, 0.1005, 0.1517, 0.0648)

se = sd_1/sqrt(12)


print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
             aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                                        "glmnetpreselect", "RFpreselect", "Rfepreselect")), y = Spearman_Correlation, fill = Deconvolution)) +
        geom_bar(stat = 'identity', position = 'dodge')+
        #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
        geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
        labs(x = "Feature selection")+
        labs(y = "Average Spearman correlation")+
        labs(fill = "Deconvolution algorithms")+
        coord_cartesian(ylim=c(0.5,1)))

dev.off()




#####################
### Figure 2
#####################
library(VennDiagram)
grid.newpage()
draw.triple.venn(area1 = 600, area2= 600, area3 = 946, n12 = 576, n23 = 350, n13 = 347, 
                 n123 = 345, category = c("oneVsAllttest", "oneVsAllLimma", "glmnetpreselect"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))






#####################
### Figure 3
#####################
pdf(file = "Figure3.pdf",  
    width = 16,
    height = 8) 


Houseman <- c(0.978, 0.947, 0.947, 0.942)
RPC <- c(0.978, 0.938, 0.942, 0.937)
CBS <- c(0.956, 0.950, 0.952, 0.918)
MethylResolver <- c(0.978, 0.920, 0.928, 0.936)
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","glmnetpreselect", "oneVsAllttest","glmnetpreselect"), 
                  facet = c("A. EPIC reference library", "A. EPIC reference library", "B. 450k reference library", "B. 450k reference library"), 
                  Houseman = Houseman, RPC = RPC, CBS = CBS, MethylResolver = MethylResolver)


sd_1 <- c(0.0299, 0.0658, 0.0662, 0.0685,
          0.0299, 0.0743, 0.0640, 0.0660,
          0.0639, 0.0470, 0.0679, 0.1150,
          0.0299, 0.1005, 0.1004, 0.1181)

se <- sd_1/sqrt(12)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(FeatureSelection,facet)), 
       aes(x = factor(FeatureSelection, level = c("oneVsAllttest","glmnetpreselect")), y = Spearman_Correlation, fill = Deconvolution)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
  facet_wrap(~factor(facet,level = c("A. EPIC reference library","B. 450k reference library")))+
  labs(x = "Feature selection")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Deconvolution algorithms")+
  coord_cartesian(ylim=c(0.5,1))+
  theme(strip.text.x = element_text(size = 15))




dev.off()






#####################
### Figure 5
#####################

pdf(file = "Figure5.pdf",  
    width = 16,
    height = 8) 


cor_Houseman <- c(0.969, 0.969, 0.911, 0.957, 0.952, 0.937)
cor_RPC <- c(0.983, 0.974, 0.911, 0.983, 0.971, 0.932)
cor_CBS <- c(0.964, 0.974, 0.920, 0.971, 0.952, 0.952)
cor_MethylResolver <- c(0.978, 0.962, 0.877, 0.978, 0.971, 0.944)

df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet","glmnetpreselect"),
                  Houseman = cor_Houseman, RPC = cor_RPC, CBS = cor_CBS, MethylResolver = cor_MethylResolver)

sd_1 = c(0.0521, 0.0521, 0.1116, 0.0688, 0.0679, 0.0755,
         0.0252, 0.0306, 0.1193, 0.0252, 0.0354, 0.0729,
         0.0516, 0.0306, 0.0872, 0.0482, 0.0664, 0.0675, 
         0.0502, 0.0555, 0.1323, 0.0299, 0.0354, 0.0760)

se = sd_1/sqrt(12)

print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
             aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                                        "glmnetpreselect")), y = Spearman_Correlation, fill = Deconvolution)) +
        geom_bar(stat = 'identity', position = 'dodge')+
        #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
        geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
        labs(x = "Feature selection")+
        labs(y = "Average Spearman correlation")+
        labs(fill = "Deconvolution algorithms")+
        coord_cartesian(ylim=c(0.5,1.0)))



dev.off()



#####################
### Figure 6
#####################
Houseman <- c(0.989, 0.970, 0.901, 0.639,0.258, 0.990, 0.983,0.973,0.894,0.711)
RPC <- c(0.991,0.970,0.828,0.542,0.236,0.992,0.986,0.982,0.939,0.829)
CBS <- c(0.974, 0.955,0.875,0.643,0.324,0.980,0.981,0.978,0.919,0.772)
MethylResolver <- c(0.992,0.974,0.864,0.651,0.368,0.991,0.985,0.980,0.939,0.819)
df1 <- data.frame(NonImmuneProportions = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9", "0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9"),
                  facet = c(rep("EPIC reference library", 5), rep("EPIC + Epithelial reference library", 5)),
                  Houseman = Houseman, RPC = RPC, CBS = CBS, MethylResolver = MethylResolver)


se = c(0.00116, 0.00185, 0.00404, 0.01235, 0.01702, 0.00111, 0.00156, 0.00158, 0.00557, 0.01085,  
       0.00098, 0.00189, 0.00736, 0.01340, 0.01790, 0.00098, 0.00133, 0.00121, 0.00351, 0.00707,
       0.00201, 0.00255, 0.00540, 0.01171, 0.01696, 0.00183, 0.00146, 0.00151, 0.00472, 0.01010,
       0.00093, 0.00174, 0.00560, 0.01210, 0.01675, 0.00107, 0.00150, 0.00132, 0.00353, 0.00772)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(NonImmuneProportions,facet)),
       aes(x = factor(NonImmuneProportions, level = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9")), 
           y = Spearman_Correlation, fill =Deconvolution)) +
  
  #geom_point(aes(col = Deconvolution))+ 
  #geom_line(aes(col = Deconvolution), size = 1)+
  #geom_line(group = 4, aes(col = factor(Deconvolution)))+
  geom_bar(stat = 'identity', position = 'dodge')+
  #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
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
                  facet = c(rep("EPIC reference library", 5), rep("EPIC + Epithelial reference library", 5)),
                  Houseman = Houseman, RPC = RPC, CBS = CBS, MethylResolver = MethylResolver)


se = c(0.00884,0.00819,0.01214,0.01586,0.01673,	0.01382,0.01143,0.01245,0.01676,0.01761,
       0.00324,0.00258,0.00706,0.01374,0.01759, 0.00282,0.00382,0.00293,0.00612,0.01403,
       0.00297,0.00242,0.00905,0.01483,0.01764,	0.00299,0.00356,0.00304,0.00754,0.01440,
       0.00330,0.00288,0.00612,0.01297,0.01739,0.00286,0.00373,0.00299,0.00702,0.01409)


ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(NonImmuneProportions,facet)),
       aes(x = factor(NonImmuneProportions, level = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9")), 
           y = Spearman_Correlation, fill =Deconvolution)) +
  
  #geom_point(aes(col = Deconvolution))+ 
  #geom_line(aes(col = Deconvolution), size = 1)+
  #geom_line(group = 4, aes(col = factor(Deconvolution)))+
  geom_bar(stat = 'identity', position = 'dodge')+
  #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
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
  labs(y = "Spearman correlation")+
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
  labs(y = "Spearman correlation")+
  labs(fill = "Cell Types")


######################
### Figure 10
######################
pdf(file = "Figure10.pdf",  
    width = 12,
    height = 8) 

cell_types <- unique(reference_EPIC_extend$ref_phenotype)

## coefficient of variation
cv_mat <- matrix(NA, nrow = dim(reference_EPIC_extend$ref_betamatrix)[1], ncol = 7)
for (i in 1:7){
  ids <- which(reference_EPIC_extend$ref_phenotype == cell_types[i])
  dat <- reference_EPIC_extend$ref_betamatrix[,ids]
  cv_mat[,i] <- apply(dat, 1, function(x){return(sd(x)/mean(x))})
}


df <- data.frame(cell_type = rep(cell_types, rep(866091,7)), 
                 cv = c(cv_mat[,1], cv_mat[,2], cv_mat[,3], cv_mat[,4], 
                        cv_mat[,5], cv_mat[,6], cv_mat[,7]))
ggplot(df, aes(x=cv, fill=cell_type, color = cell_type)) +
  geom_histogram(position="identity", alpha=0.4, bins = 100) + 
  labs(x = "coefficients of variation")+
  labs(y = "count of CpGs")+
  coord_cartesian(xlim=c(0,1.5))

dev.off()







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





####################
## Figure 11  (new)
####################
pdf(file = "Figure11.pdf",  
    width = 16,
    height = 8) 


Houseman <- c(0.876, 0.752, 0.781, 0.838)
RPC <- c(0.781, 0.848, 0.762, 0.876)
CBS <- c(0.743, 0.876, 0.648, 0.886)
MethylResolver <- c(0.8, 0.810, 0.848, 0.8)
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","glmnetpreselect", "oneVsAllttest","glmnetpreselect"), 
                  facet = c("EPIC reference library", "EPIC reference library", "EPIC + Epithelial reference library", "EPIC + Epithelial reference library"), 
                  Houseman = Houseman, RPC = RPC, CBS = CBS, MethylResolver = MethylResolver)

sd_1 <- c(0.0668, 0.3306, 0.2298, 0.1892, 0.2151, 0.2906, 0.2089, 0.1927,
          0.2549, 0.1221, 0.3011, 0.2012, 0.2130, 0.2884, 0.1384, 0.1869)
se = sd_1/sqrt(6)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(FeatureSelection,facet)), 
       aes(x = factor(FeatureSelection, level = c("oneVsAllttest","glmnetpreselect")), y = Spearman_Correlation, fill = Deconvolution)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")))+
  labs(x = "Feature selection")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Deconvolution algorithms")+
  coord_cartesian(ylim=c(0.5,1))

dev.off()


####################
## Figure 12  (new)
####################
pdf(file = "Figure12.pdf",  
    width = 16,
    height = 8) 


Houseman <- c(0.590, 0.781, 0.270, 0.362)
RPC <- c(0.793, 0.957, 0.224, -0.005)
CBS <- c(0.674, 0.850, 0.231, 0.187)
MethylResolver <- c(0.581, 0.895, 0.392, -0.019)
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","glmnetpreselect", "oneVsAllttest","glmnetpreselect"), 
                  facet = c("EPIC reference library", "EPIC reference library", "EPIC + Epithelial reference library", "EPIC + Epithelial reference library"), 
                  Houseman = Houseman, RPC = RPC, CBS = CBS, MethylResolver = MethylResolver)

sd_1 <- c(0.1711, 0.1749, 0.3077, 0.3139, 0.1110, 0.0689, 0.1320, 0.0743,
          0.2429, 0.1177, 0.1825, 0.1519, 0.3226, 0.1221, 0.1775, 0.0812)
se = sd_1/sqrt(6)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(FeatureSelection,facet)), 
       aes(x = factor(FeatureSelection, level = c("oneVsAllttest","glmnetpreselect")), y = Spearman_Correlation, fill = Deconvolution)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + Epithelial reference library")))+
  labs(x = "Feature selection")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Deconvolution algorithms")+
  coord_cartesian(ylim=c(-0.1,1))

dev.off()



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



###################
## Figure 15
###################

load("res_133395.RData")
library(ggplot2)
library(tidyr)
pdf(file = "Plot3.pdf",  
    width = 10,
    height = 8) 
df <- gather(gender_dat_blood, series,value,-group)
ggplot(df) + geom_boxplot(aes(series ,value,color=group)) +
  xlab('cell types')+
  ylab('proportions') +
  ggtitle("GSE133395-RPC-onevsAllttest")
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
  ggtitle("GSE112308-RPC-onevsAllttest")
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





####################
## Supp Figure 8
####################
Houseman <- c(0.991,0.917,0.734,0.504,0.201,0.986,0.984,0.975,0.910,0.722)
RPC <- c(0.992,0.914,0.710,0.472,0.154,0.991,0.985,0.976,0.908,0.728)
CBS <- c(0.976,0.885,0.686,0.461,0.126,0.969,0.971,0.963,0.892,0.703)
df1 <- data.frame(NonImmuneProportions = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9", "0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9"),
                  facet = c(rep("EPIC reference library", 5), rep("EPIC + cfDNA reference library",5)),
                  Houseman = Houseman, RPC = RPC, CBS = CBS)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(NonImmuneProportions,facet)), 
       aes(x = factor(NonImmuneProportions, level = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9")), y = Spearman_Correlation, fill = Deconvolution)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + cfDNA reference library")))+
  labs(x = "Non-Immune (cfDNA) Proportions")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Deconvolution algorithms")



####################
## Supp Figure 9
####################

### within cell type Houseman
Bcell <- c(0.998,0.995,0.987,0.968,0.837,0.998,0.997,0.996,0.986,0.866)
CD4T <- c(0.996,0.994,0.982,0.900,0.827,0.994,0.994,0.982,0.924,0.829)
CD8T <- c(0.995,0.992,0.970,0.943,0.839,0.992,0.988,0.972,0.934,0.764)
Mono <- c(0.996,0.993,0.975,0.912,0.788,0.996,0.997,0.990,0.932,0.814)
Neu <- c(0.996,0.974,0.805,0.465,0.451,0.997,0.995,0.983,0.928,0.480)
NK <- c(0.997,0.994,0.978,0.957,0.865,0.997,0.997,0.995,0.977,0.923)
cfDNA <- c(NA,NA,NA,NA,NA,NA,0.975,0.978,0.956,0.796)
df1 <- data.frame(NonImmuneProportions = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9", "0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9"),
                  facet = c(rep("EPIC reference library", 5), rep("EPIC + cfDNA reference library",5)),
                  Bcell = Bcell, CD4T = CD4T, CD8T = CD8T, Mono = Mono, Neu = Neu, NK = NK, Nonimmune_cfDNA = cfDNA)
ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(NonImmuneProportions,facet)),
       aes(x = factor(NonImmuneProportions, level = c("0","0.1-0.2", "0.2-0.5","0.5-0.8","0.8-0.9")), y = Spearman_Correlation, fill = Deconvolution)) +
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  facet_wrap(~factor(facet,level = c("EPIC reference library","EPIC + cfDNA reference library")),dir = "v")+
  labs(x = "Non-immune (cfDNA) proportions")+
  labs(y = "Average Spearman correlation")+
  labs(fill = "Cell Types")

