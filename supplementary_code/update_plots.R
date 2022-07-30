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
                  Houseman = cor_Houseman, RPC = cor_RPC, CBS = cor_CBS, MethylResolver = cor_MethylResolver, ARIC = cor_ARIC, 
                  TOAST_RefFree_1 = cor_TOAST_1, TOAST_RefFree_2 = cor_TOAST_2)


sd_1 <- c(0.0031, 0.0031, 0.0037, 0.0036, 0.0049, 0.0051, 0.0071, 0.0024,
          0.0058, 0.0060, 0.0045, 0.0024, 0.0043, 0.0085, 0.0092, 0.0020,
          0.0058, 0.0062, 0.0057, 0.0083, 0.0065, 0.0044, 0.0061, 0.0081,
          0.0078, 0.0059, 0.0060, 0.0031, 0.0028, 0.0065, 0.0109, 0.0024,
          0.0064, 0.0062, 0.0084, 0.0065, 0.0025, 0.0084, 0.0109, 0.0058,0.1888, 
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
        labs(y = "Average Spearman correlation")+
        labs(fill = "Deconvolution algorithms"))
        
dev.off()
