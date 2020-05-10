### GSE122126

#Purified pancreatic acinar cells, pancreatic duct cells, 
# pancreatic beta cells, vascular endothelial cells and colon epithelial cells of 450k and EPIC data.
library(GEOquery)
library(minfi)
getGEOSuppFiles("GSE122126")
untar("GSE122126/GSE122126_RAW.tar", exdir = "GSE122126/idat")
head(list.files("GSE122126/idat", pattern = "idat"))

idatFiles <- list.files("GSE122126/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)
rgSet <- read.metharray.exp("GSE122126/idat", force = TRUE)
grSet <- preprocessNoob(rgSet, dyeMethod = "single")
betaMat <- getBeta(grSet)
#866091     90

geoMat <- getGEO("GSE122126")
pD.all <- pData(geoMat[[2]])
pD <- pD.all[, c("title", "geo_accession", "age:ch1", "disease state:ch1", "gender:ch1", "sample type:ch1")]
# add phenotype data
sampleNames(rgSet) <- substr(sampleNames(rgSet), 1, 10)
pD <- pD[sampleNames(rgSet),]
pD <- as(pD, "DataFrame")
pData(rgSet) <- pD



