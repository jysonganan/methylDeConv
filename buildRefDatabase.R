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




### 
## Reformat FlowSorted.450k; with preprocessNoob

library(FlowSorted.Blood.450k)
CellLines.matrix = NULL
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
keep <- which(ref_phenotype %in% cellTypes)
ref_betamatrix <- ref_betamatrix[,keep]
ref_phenotype <- ref_phenotype[keep]

source("refCompTableProbeSelection.R")
compTable <- ref_compTable(ref_betamatrix, ref_phenotype)
probes_oneVsAllttest <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype, probeSelect = "both", MaxDMRs = 100)
probes_oneVsAllLimma <- ref_probe_selection_oneVsAllLimma(ref_betamatrix, ref_phenotype, probeSelect = "both")
## the above two methods: intersection 538/600 probes
probes_pairwiseLimma <- ref_probe_selection_pairwiseLimma(ref_betamatrix, ref_phenotype)  #681 probes selected
probes_multiclassGlmnet <- ref_probe_selection_multiclassGlmnet(ref_betamatrix, ref_phenotype)

# dat_450k <- pickCompProbes(preprocessQuantile(FlowSorted.Blood.450k), cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),
#                            probeSelect = "both")
