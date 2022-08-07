library(tidyr)
library(ggplot2)

setwd("./Downloads/")
pdf(file = "Figure1.pdf",  
    width = 16,
    height = 8) 


cor_Houseman <- c(0.978, 0.969, 0.950, 0.966, 0.947, 0.947, 0.831, 0.964)
cor_RPC <- c(0.978, 0.978, 0.828, 0.983, 0.952, 0.938, 0.800, 0.959)
cor_CBS <- c(0.956, 0.962, 0.839, 0.952, 0.942, 0.950, 0.846, 0.952)
cor_MethylResolver <- c(0.978, 0.983, 0.868, 0.978, 0.956, 0.920, 0.788, 0.955)
cor_ARIC <- c(0.993, 0.988, 0.948, 0.988, 0.971, 0.969, 0.831, 0.939)
cor_TOAST_1 <- c(0.652, 0.621, 0.461, 0.564, 0.646, 0.665, 0.452, 0.549)
cor_TOAST_2 <- c(0.712, 0.677, 0.411, 0.564, 0.646, 0.665, 0.452, 0.624)
  
  
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                       "glmnetpreselect", "RFpreselect", "Rfepreselect"),
                  Houseman = cor_Houseman, RPC = cor_RPC, CBS = cor_CBS, MethylResolver = cor_MethylResolver, ARIC = cor_ARIC, 
                  TOAST_RefFree_1 = cor_TOAST_1, TOAST_RefFree_2 = cor_TOAST_2)


sd_1 <- c(0.0299, 0.0521, 0.0530, 0.0548, 0.0662, 0.0658, 0.1200, 0.0516,
          0.0299, 0.0299, 0.1766, 0.0252, 0.0679, 0.0743, 0.1670, 0.0659,
          0.0639, 0.0471, 0.1831, 0.0501, 0.0812, 0.0470, 0.1028, 0.0501,
          0.0299, 0.0282, 0.1804, 0.0299, 0.0692, 0.1005, 0.1517, 0.0648,
          0.0172, 0.0255, 0.0885, 0.0222, 0.0354, 0.0506, 0.1592, 0.0902,
          0.3488, 0.3619, 0.3764, 0.3385, 0.3349, 0.2944, 0.3611, 0.3014,
          0.2586, 0.2844, 0.4214, 0.3385, 0.3349, 0.2944, 0.3611, 0.2639)

se = sd_1/sqrt(12)


print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
             aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                                        "glmnetpreselect", "RFpreselect", "Rfepreselect")), y = Spearman_Correlation, fill = Deconvolution)) +
        geom_bar(stat = 'identity', position = 'dodge')+
        #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
        geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
        labs(x = "Feature selection")+
        labs(y = "Average Spearman correlation")+
        labs(fill = "Deconvolution algorithms"))
        
dev.off()











pdf(file = "Figure1_RMSE.pdf",  
    width = 16,
    height = 8) 


RMSE_Houseman <- c(0.0130, 0.0129, 0.0195, 0.0141, 0.0146, 0.0175, 0.0219, 0.0139)
RMSE_RPC <- c(0.0147, 0.0142, 0.0288, 0.0115, 0.0134, 0.0194, 0.0256, 0.0144)
RMSE_CBS <- c(0.0163, 0.0157, 0.0262, 0.0144, 0.0201, 0.0191, 0.0291, 0.0162)
RMSE_MethylResolver <- c(0.0148, 0.0131, 0.0286, 0.0123, 0.0121, 0.0185, 0.0275, 0.0155)
RMSE_ARIC <- c(0.0168, 0.0168, 0.0216, 0.0158, 0.0130, 0.0183, 0.0342, 0.0204)
RMSE_TOAST_1 <- c(0.1732, 0.1718, 0.2223, 0.2269, 0.1962, 0.1802, 0.2312, 0.1904)
RMSE_TOAST_2 <- c(0.1888, 0.1961, 0.2418, 0.2269, 0.1962, 0.1802, 0.2312, 0.2060)
  
  
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                       "glmnetpreselect", "RFpreselect", "Rfepreselect"),
                  Houseman = RMSE_Houseman, RPC = RMSE_RPC, CBS = RMSE_CBS, MethylResolver = RMSE_MethylResolver, ARIC = RMSE_ARIC, 
                  TOAST_RefFree_1 = RMSE_TOAST_1, TOAST_RefFree_2 = RMSE_TOAST_2)


sd_1 <- c(0.0031, 0.0031, 0.0037, 0.0036, 0.0049, 0.0051, 0.0071, 0.0024,
          0.0058, 0.0060, 0.0045, 0.0024, 0.0043, 0.0085, 0.0092, 0.0020,
          0.0058, 0.0062, 0.0057, 0.0083, 0.0065, 0.0044, 0.0061, 0.0081,
          0.0078, 0.0059, 0.0060, 0.0031, 0.0028, 0.0065, 0.0109, 0.0024,
          0.0064, 0.0062, 0.0084, 0.0065, 0.0025, 0.0084, 0.0109, 0.0058,
          0.0845, 0.0863, 0.0888, 0.0768, 0.0844, 0.0722, 0.0952, 0.0912,
          0.0909, 0.0969, 0.1186, 0.0768, 0.0844, 0.0722, 0.0952, 0.0835)

se = sd_1/sqrt(12)

print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
             aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                                        "glmnetpreselect", "RFpreselect", "Rfepreselect")), y = Spearman_Correlation, fill = Deconvolution)) +
        geom_bar(stat = 'identity', position = 'dodge')+
        #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
        geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
        labs(x = "Feature selection")+
        labs(y = "Average RMSE")+
        labs(fill = "Deconvolution algorithms"))
        
dev.off()







pdf(file = "Figure1_SMAPE.pdf",  
    width = 16,
    height = 8) 


SMAPE_Houseman <- c(0.1657, 0.1667, 0.2348, 0.1628, 0.1360, 0.2078, 0.2693, 0.1482)
SMAPE_RPC <- c(0.1509, 0.1497, 0.3088, 0.1185, 0.1237, 0.2040, 0.2951, 0.1367)
SMAPE_CBS <- c(0.1591, 0.1560, 0.2806, 0.1485, 0.1674, 0.2141, 0.3111, 0.1642)
SMAPE_MethylResolver <- c(0.1566, 0.1409, 0.2964, 0.1123, 0.1181, 0.1808, 0.3273, 0.1457)
SMAPE_ARIC <- c(0.1398, 0.1418, 0.1858, 0.1255, 0.1254, 0.1590, 0.2505, 0.1374)
SMAPE_TOAST_1 <- c(1.2084, 1.2098, 1.3597, 1.3964, 1.3637, 1.2809, 1.4500, 1.2372)
SMAPE_TOAST_2 <- c(1.3131, 1.3176, 1.4966, 1.3964, 1.3637, 1.2809, 1.4500, 1.3643)
  
  
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                       "glmnetpreselect", "RFpreselect", "Rfepreselect"),
                  Houseman = SMAPE_Houseman, RPC = SMAPE_RPC, CBS = SMAPE_CBS, MethylResolver = SMAPE_MethylResolver, ARIC = SMAPE_ARIC, 
                  TOAST_RefFree_1 = SMAPE_TOAST_1, TOAST_RefFree_2 = SMAPE_TOAST_2)


sd_1 <- c(0.1022, 0.0993, 0.1320, 0.1150, 0.1074, 0.1464, 0.2127, 0.0924,
          0.1174, 0.1168, 0.1180, 0.0918, 0.1056, 0.1486, 0.2243, 0.0952,
          0.0964, 0.0950, 0.1229, 0.0879, 0.1005, 0.1346, 0.1898, 0.1087,
          0.1656, 0.1255, 0.1198, 0.0947, 0.0991, 0.1363, 0.2557, 0.0917,
          0.1215, 0.1208, 0.1213, 0.0984, 0.0945, 0.1267, 0.0991, 0.0813,
          0.2141, 0.2080, 0.1634, 0.1469, 0.2086, 0.1986, 0.1926, 0.1669,
          0.2160, 0.2393, 0.2096, 0.1469, 0.2086, 0.1986, 0.1926, 0.1762)

se = sd_1/sqrt(12)

print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
             aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                                        "glmnetpreselect", "RFpreselect", "Rfepreselect")), y = Spearman_Correlation, fill = Deconvolution)) +
        geom_bar(stat = 'identity', position = 'dodge')+
        #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
        geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
        labs(x = "Feature selection")+
        labs(y = "Average SMAPE")+
        labs(fill = "Deconvolution algorithms"))
        
dev.off()









## within cell-type

Houseman <- c(0.978, 0.947, 0.947, 0.942)
RPC <- c(0.978, 0.938, 0.942, 0.937)
CBS <- c(0.956, 0.950, 0.952, 0.918)
MethylResolver <- c(0.978, 0.920, 0.928, 0.936)
ARIC <- c(0.993, 0.969, 0.937, 0.885)
TOAST_1 <- c(0.652, 0.665, 0.671, 0.603)
TOAST_2 <- c(0.712, 0.665, 0.635, 0.603)


pdf(file = "Figure3.pdf",  
    width = 16,
    height = 8) 

df1 <- data.frame(FeatureSelection = c("oneVsAllttest","glmnetpreselect", "oneVsAllttest","glmnetpreselect"), 
                  facet = c("A. EPIC reference library", "A. EPIC reference library", "B. 450k reference library", "B. 450k reference library"), 
                  Houseman = Houseman, RPC = RPC, CBS = CBS, MethylResolver = MethylResolver,
                  ARIC = ARIC, TOAST_1 = TOAST_1, TOAST_2 = TOAST_2)


sd_1 <- c(0.0299, 0.0658, 0.0662, 0.0685,
          0.0299, 0.0743, 0.0640, 0.0660,
          0.0639, 0.0470, 0.0679, 0.1150,
          0.0299, 0.1005, 0.1004, 0.1181,
          0.0172, 0.0506, 0.0949, 0.1161,
          0.3488, 0.2944, 0.2235, 0.3179,
          0.2586, 0.2944, 0.3360, 0.3179)


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







pdf(file = "Figure3_RMSE.pdf",  
    width = 16,
    height = 8) 


Houseman <- c(0.0130, 0.0175, 0.0180, 0.0198)
RPC <- c(0.0147, 0.0194, 0.0204, 0.0216)
CBS <- c(0.0163, 0.0191, 0.0223, 0.0194)
MethylResolver <- c(0.0148, 0.0185, 0.0209, 0.0239)
ARIC <- c(0.0168, 0.0183, 0.0250, 0.0295)
TOAST_1 <- c(0.1732, 0.1802, 0.1827, 0.1764)
TOAST_2 <- c(0.1888, 0.1802, 0.2069, 0.1764)




df1 <- data.frame(FeatureSelection = c("oneVsAllttest","glmnetpreselect", "oneVsAllttest","glmnetpreselect"), 
                  facet = c("A. EPIC reference library", "A. EPIC reference library", "B. 450k reference library", "B. 450k reference library"), 
                  Houseman = Houseman, RPC = RPC, CBS = CBS, MethylResolver = MethylResolver,
                  ARIC = ARIC, TOAST_1 = TOAST_1, TOAST_2 = TOAST_2)


sd_1 <- c(0.0031, 0.0051, 0.0039, 0.0044,
          0.0058, 0.0085, 0.0061, 0.0062,
          0.0058, 0.0044, 0.0091, 0.0064,
          0.0078, 0.0065, 0.0070, 0.0095,
          0.0064, 0.0084, 0.0074, 0.0089,
          0.0845, 0.0722, 0.0664, 0.0538,
          0.0909, 0.0722, 0.0824, 0.0538)


se <- sd_1/sqrt(12)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(FeatureSelection,facet)), 
       aes(x = factor(FeatureSelection, level = c("oneVsAllttest","glmnetpreselect")), y = Spearman_Correlation, fill = Deconvolution)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
  facet_wrap(~factor(facet,level = c("A. EPIC reference library","B. 450k reference library")))+
  labs(x = "Feature selection")+
  labs(y = "Average RMSE")+
  labs(fill = "Deconvolution algorithms")+
  theme(strip.text.x = element_text(size = 15))

dev.off()








pdf(file = "Figure3_SMAPE.pdf",  
    width = 16,
    height = 8) 


Houseman <- c(0.1657, 0.2078, 0.1595, 0.1561)
RPC <- c(0.1509, 0.2040, 0.1745, 0.1457)
CBS <- c(0.1591, 0.2141, 0.1573, 0.1401)
MethylResolver <- c(0.1566, 0.1808, 0.1635, 0.1483)
ARIC <- c(0.1398, 0.1590, 0.1993, 0.1640)
TOAST_1 <- c(1.2084, 1.2809, 1.2845, 1.2498)
TOAST_2 <- c(1.3131, 1.2809, 1.3696, 1.2498)




df1 <- data.frame(FeatureSelection = c("oneVsAllttest","glmnetpreselect", "oneVsAllttest","glmnetpreselect"), 
                  facet = c("A. EPIC reference library", "A. EPIC reference library", "B. 450k reference library", "B. 450k reference library"), 
                  Houseman = Houseman, RPC = RPC, CBS = CBS, MethylResolver = MethylResolver,
                  ARIC = ARIC, TOAST_1 = TOAST_1, TOAST_2 = TOAST_2)


sd_1 <- c(0.1022, 0.1464, 0.1097, 0.1026,
          0.1174, 0.1486, 0.0948, 0.0540,
          0.0964, 0.1346, 0.0427, 0.0457,
          0.1656, 0.1363, 0.0756, 0.0350,
          0.1215, 0.1267, 0.1021, 0.0479,
          0.2141, 0.1986, 0.2190, 0.2259,
          0.2160, 0.1986, 0.2282, 0.2259)


se <- sd_1/sqrt(12)

ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -c(FeatureSelection,facet)), 
       aes(x = factor(FeatureSelection, level = c("oneVsAllttest","glmnetpreselect")), y = Spearman_Correlation, fill = Deconvolution)) + 
  geom_bar(stat = 'identity', position = 'dodge')+
  #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
  facet_wrap(~factor(facet,level = c("A. EPIC reference library","B. 450k reference library")))+
  labs(x = "Feature selection")+
  labs(y = "Average SMAPE")+
  labs(fill = "Deconvolution algorithms")+
  theme(strip.text.x = element_text(size = 15))

dev.off()
















pdf(file = "Figure4.pdf",  
    width = 16,
    height = 8) 


cor_Houseman <- c(0.969, 0.969, 0.911, 0.957, 0.952, 0.937, 0.907, 0.899)
cor_RPC <- c(0.983, 0.974, 0.911, 0.983, 0.971, 0.932, 0.899, 0.885)
cor_CBS <- c(0.964, 0.974, 0.920, 0.971, 0.952, 0.952, 0.868, 0.911)
cor_MethylResolver <- c(0.978, 0.962, 0.877, 0.978, 0.971, 0.944, 0.885, 0.911)
cor_ARIC <- c(0.981, 0.983, 0.939, 0.974, 0.952, 0.976, 0.861, 0.868)
cor_TOAST_1 <- c(0.366, 0.391, 0.372, 0.281, 0.383, 0.508, 0.414, 0.502)
cor_TOAST_2 <- c(0.366, 0.391, 0.325, 0.281, 0.383, 0.508, 0.312, 0.502)
  
  
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                       "glmnetpreselect", "RFpreselect", "Rfepreselect"),
                  Houseman = cor_Houseman, RPC = cor_RPC, CBS = cor_CBS, MethylResolver = cor_MethylResolver, ARIC = cor_ARIC, 
                  TOAST_RefFree_1 = cor_TOAST_1, TOAST_RefFree_2 = cor_TOAST_2)


sd_1 <- c(0.0521, 0.0521, 0.1116, 0.0688, 0.0679, 0.0755, 0.0752, 0.0965,
          0.0252, 0.0306, 0.1193, 0.0252, 0.0354, 0.0729, 0.0969, 0.0855,
          0.0516, 0.0306, 0.0872, 0.0482, 0.0664, 0.0675, 0.1488, 0.0774,
          0.0502, 0.0555, 0.1323, 0.0299, 0.0354, 0.0760, 0.1327, 0.0836,
          0.0340, 0.0252, 0.0961, 0.0506, 0.0679, 0.0351, 0.1922, 0.2302,
          0.3469, 0.3382, 0.2564, 0.4497, 0.3937, 0.3972, 0.2827, 0.2950,
          0.3469, 0.3382, 0.4171, 0.4497, 0.3937, 0.3972, 0.4793, 0.2950)

se = sd_1/sqrt(12)


print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
             aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                                        "glmnetpreselect", "RFpreselect", "Rfepreselect")), y = Spearman_Correlation, fill = Deconvolution)) +
        geom_bar(stat = 'identity', position = 'dodge')+
        #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
        geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
        labs(x = "Feature selection")+
        labs(y = "Average Spearman correlation")+
        labs(fill = "Deconvolution algorithms"))
        
dev.off()











pdf(file = "Figure4_RMSE.pdf",  
    width = 16,
    height = 8) 


RMSE_Houseman <- c(0.0148, 0.0147, 0.0179, 0.0150, 0.0140, 0.0189, 0.0171, 0.0207)
RMSE_RPC <- c(0.0146, 0.0152, 0.0177, 0.0124, 0.0121, 0.0199, 0.0184, 0.0204)
RMSE_CBS <- c(0.0142, 0.0152, 0.0169, 0.0113, 0.0108, 0.0153, 0.0257, 0.0169)
RMSE_MethylResolver <- c(0.0143, 0.0152, 0.0192, 0.0139, 0.0110, 0.0183, 0.0192, 0.0193)
RMSE_ARIC <- c(0.0192, 0.0196, 0.0209, 0.0162, 0.0124, 0.0180, 0.0434, 0.0372)
RMSE_TOAST_1 <- c(0.2187, 0.2146, 0.2062, 0.2133, 0.2121, 0.2134, 0.2496, 0.2141)
RMSE_TOAST_2 <- c(0.2187, 0.2146, 0.2239, 0.2133, 0.2121, 0.2134, 0.2499, 0.2141)
  
  
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                       "glmnetpreselect", "RFpreselect", "Rfepreselect"),
                  Houseman = RMSE_Houseman, RPC = RMSE_RPC, CBS = RMSE_CBS, MethylResolver = RMSE_MethylResolver, ARIC = RMSE_ARIC, 
                  TOAST_RefFree_1 = RMSE_TOAST_1, TOAST_RefFree_2 = RMSE_TOAST_2)


sd_1 <- c(0.0033, 0.0038, 0.0032, 0.0033, 0.0028, 0.0050, 0.0033, 0.0024,
          0.0062, 0.0061, 0.0042, 0.0030, 0.0023, 0.0081, 0.0040, 0.0021,
          0.0037, 0.0041, 0.0035, 0.0062, 0.0025, 0.0062, 0.0086, 0.0035,
          0.0067, 0.0073, 0.0042, 0.0043, 0.0047, 0.0064, 0.0042, 0.0044,
          0.0079, 0.0082, 0.0066, 0.0053, 0.0027, 0.0084, 0.0238, 0.0210,
          0.0864, 0.0875, 0.0958, 0.0958, 0.1010, 0.0958, 0.0963, 0.0894,
          0.0864, 0.0875, 0.1084, 0.0958, 0.1010, 0.0958, 0.0945, 0.0894)

se = sd_1/sqrt(12)

print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
             aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                                        "glmnetpreselect", "RFpreselect", "Rfepreselect")), y = Spearman_Correlation, fill = Deconvolution)) +
        geom_bar(stat = 'identity', position = 'dodge')+
        #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
        geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
        labs(x = "Feature selection")+
        labs(y = "Average RMSE")+
        labs(fill = "Deconvolution algorithms"))
        
dev.off()







pdf(file = "Figure4_SMAPE.pdf",  
    width = 16,
    height = 8) 


SMAPE_Houseman <- c(0.1540, 0.1532, 0.2393, 0.1566, 0.1236, 0.1967, 0.1843, 0.1941)
SMAPE_RPC <- c(0.1279, 0.1350, 0.2300, 0.1179, 0.1132, 0.1715, 0.1602, 0.1961)
SMAPE_CBS <- c(0.1353, 0.1429, 0.2187, 0.1286, 0.1128, 0.1527, 0.2009, 0.1972)
SMAPE_MethylResolver <- c(0.1233, 0.1353, 0.2385, 0.1134, 0.1407, 0.1360, 0.1784, 0.2045)
SMAPE_ARIC <- c(0.1649, 0.1562, 0.2053, 0.1264, 0.1441, 0.1452, 0.2762, 0.2414)
SMAPE_TOAST_1 <- c(1.4729, 1.4618, 1.4758, 1.5140, 1.4847, 1.5348, 1.5889, 1.4619)
SMAPE_TOAST_2 <- c(1.4729, 1.4618, 1.6317, 1.5140, 1.4847, 1.5348, 1.5887, 1.4619)
  
  
df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                       "glmnetpreselect", "RFpreselect", "Rfepreselect"),
                  Houseman = SMAPE_Houseman, RPC = SMAPE_RPC, CBS = SMAPE_CBS, MethylResolver = SMAPE_MethylResolver, ARIC = SMAPE_ARIC, 
                  TOAST_RefFree_1 = SMAPE_TOAST_1, TOAST_RefFree_2 = SMAPE_TOAST_2)


sd_1 <- c(0.1071, 0.1119, 0.1291, 0.0999, 0.1048, 0.1330, 0.1218, 0.1047,
          0.1087, 0.1103, 0.1457, 0.0942, 0.0998, 0.1309, 0.0431, 0.1135,
          0.1117, 0.1132, 0.1283, 0.0874, 0.1036, 0.1232, 0.0419, 0.1175,
          0.1187, 0.1215, 0.1577, 0.0953, 0.1261, 0.1017, 0.1058, 0.1101,
          0.1311, 0.1112, 0.1243, 0.0959, 0.1029, 0.1153, 0.1654, 0.1543,
          0.2013, 0.2074, 0.2390, 0.2028, 0.2250, 0.2172, 0.2169, 0.2231,
          0.2013, 0.2074, 0.2209, 0.2028, 0.2250, 0.2172, 0.2672, 0.2231)

se = sd_1/sqrt(12)

print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
             aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet","multiGlmnet",
                                                        "glmnetpreselect", "RFpreselect", "Rfepreselect")), y = Spearman_Correlation, fill = Deconvolution)) +
        geom_bar(stat = 'identity', position = 'dodge')+
        #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
        geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
        labs(x = "Feature selection")+
        labs(y = "Average SMAPE")+
        labs(fill = "Deconvolution algorithms"))
        
dev.off()







#################### DL feature selection results



pdf(file = "AddDL_feature_selection.pdf",  
    width = 16,
    height = 8) 
library(tidyr)
library(ggplot2)

cor_Houseman <- c(0.969, 0.969, 0.911, 0.938, 0.937)
cor_RPC <- c(0.983, 0.974, 0.911, 0.983, 0.932)
cor_CBS <- c(0.964, 0.974, 0.920, 0.974, 0.952)
cor_MethylResolver <- c(0.978, 0.962, 0.877, 0.952, 0.944)


df1 <- data.frame(FeatureSelection = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet + DualNet", "glmnetpreselect"),
                  Houseman = cor_Houseman, RPC = cor_RPC, CBS = cor_CBS, MethylResolver = cor_MethylResolver)

sd_1 <- c(0.0521, 0.0521, 0.1116, 0.0598, 0.0755,
          0.0252, 0.0306, 0.1193, 0.0252, 0.0729, 
          0.0516, 0.0306, 0.0872, 0.0280, 0.0675, 
          0.0502, 0.0555, 0.1323, 0.0608, 0.0760)

se = sd_1/sqrt(12)

print(ggplot(data = df1 %>% gather(Deconvolution, Spearman_Correlation, -FeatureSelection),
             aes(x = factor(FeatureSelection, level = c("oneVsAllttest","oneVsAllLimma","pairwiseLimma","pairwiseGlmnet + DualNet",
                                                        "glmnetpreselect")), y = Spearman_Correlation, fill = Deconvolution)) +
        geom_bar(stat = 'identity', position = 'dodge')+
        #geom_text(aes(label= round(Spearman_Correlation,2)), position = position_dodge(0.9))+
        geom_errorbar(aes(ymin=Spearman_Correlation-se, ymax=Spearman_Correlation+se), width=.5, position=position_dodge(.9)) +
        labs(x = "Feature selection")+
        labs(y = "Average Spearman correlation")+
        labs(fill = "Deconvolution algorithms"))

dev.off()





pdf(file = "AddDL_feature_selection_vennPlot.pdf",  
    width = 16,
    height = 8)

library(VennDiagram)
grid.newpage()


draw.quad.venn(area1 = 700, area2= 700, area3 = 1000, area4 = 1049, n12 = 647, n13 = 10, n14 = 389, n23 = 10, n24 = 402,
               n34 = 16, n123 = 10, n124 = 389, n134 = 8, n234 = 8, n1234 = 8, category = c("oneVsAllttest", "oneVsAllLimma", 
                                                                "pairwiseLimma + DualNet", "glmnetpreselect"), lty = "blank", 
               fill = c("skyblue", "pink1", "mediumorchid", "red")) 

dev.off()



