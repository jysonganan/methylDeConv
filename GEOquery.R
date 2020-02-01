## 40 samples from individuals hospitalized with acute mania as well as unaffected controls.

library(GEOquery)

# supplemental files to GEO series (GSE), GEO platforms (GPL) and GEO samples (GSM)
getGEOSuppFiles("GSE68777")
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
head(list.files("GSE68777/idat", pattern = "idat"))

idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

library(minfi)
#targets <- read.metharray.sheet("GSE68777/idat") #however, no sample sheet!
rgSet <- read.metharray.exp("GSE68777/idat")
# 450k, (622399,40)




### phenotype data stored in GEO
geoMat <- getGEO("GSE68777")
pD.all <- pData(geoMat[[1]])
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)

names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis: ", "", pD$group)
pD$sex <- sub("^Sex: ", "", pD$sex)

## add phenotype data
sampleNames(rgSet) <- sub(".*_5", "5", sampleNames(rgSet))

rownames(pD) <- pD$title
pD <- pD[sampleNames(rgSet),]
pD <- as(pD, "DataFrame")
pData(rgSet) <- pD
rgSet







##
grSet <- preprocessQuantile(rgSet)
getBeta(grSet)[1:3,1:3]
head(getIslandStatus(grSet))
