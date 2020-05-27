## matched 450k and EPIC samples GSE111165

library(GEOquery)
library(minfi)
getGEOSuppFiles("GSE111165")
untar("GSE111165/GSE111165_RAW.tar", exdir = "GSE111165/idat")
head(list.files("GSE111165/idat", pattern = "idat"))

idatFiles <- list.files("GSE111165/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

geoMat <- getGEO("GSE111165")



## for 450k 
pD.all <- pData(geoMat[[1]])
