library(FlowSorted.Blood.450k)
wh.WBC <- which(FlowSorted.Blood.450k$CellType == "WBC")
wh.PBMC <- which(FlowSorted.Blood.450k$CellType == "PBMC")
RGset <- FlowSorted.Blood.450k[, c(wh.WBC, wh.PBMC)]
sampleNames(RGset) <- paste(RGset$CellType, c(seq(along = wh.WBC), seq(along = wh.PBMC)), sep = "_")

counts <- estimateCellCounts(RGset, meanPlot = FALSE)
round(counts, 2)

res <- estimateCellCounts(RGset, meanPlot = FALSE, returnAll = TRUE)
## res[[1]] component table 12*6
## res[[2]] normalized reference dataset 485512*11 (Fstat, p.value, CD8T, CD4T, NK, ...)
## res[[3]] normalized input user data 485512*12

source("projectCellType.R")
## without the preselection of probes, using all 485512 probes
res2 <- projectCellType(getBeta(res[[3]]), as.matrix(res[[2]][,3:8]))
sum((res[[1]]-res2)^2)/72 #0.00429997
library(hydroGOF)
mse(res[[1]],res2)
#CD8T         CD4T           NK        Bcell         Mono         Gran 
#0.0131918574 0.0092794638 0.0010565830 0.0003860891 0.0015351611 0.0003506676 
mean(mse(res[[1]],res2)) #0.00429997



## if use the FlowSorted.Blood.450k.compTable (456655 probes) 
res3 <- projectCellType(getBeta(res[[3]])[rownames(FlowSorted.Blood.450k.compTable),], 
                        as.matrix(FlowSorted.Blood.450k.compTable[,3:8]))
## if use FlowSorted.Blood.450k.JaffeModelPars (600 probes)
res4 <- projectCellType(getBeta(res[[3]])[rownames(FlowSorted.Blood.450k.JaffeModelPars),], 
                        as.matrix(FlowSorted.Blood.450k.JaffeModelPars))
mse(res[[1]],res4)
#CD8T         CD4T           NK        Bcell         Mono         Gran 
#2.489272e-03 9.331895e-04 1.389416e-04 1.089519e-03 5.481756e-05 4.514935e-04 
mean(mse(res[[1]],res4))  # 0.0008595389


## preprocess RGset with minfi
library(minfi)
manifest <- getManifest(RGset)
manifest
#IlluminaMethylationManifest object
#Annotation
#array: IlluminaHumanMethylation450k
#Number of type I probes: 135476 
#Number of type II probes: 350036 
#Number of control probes: 850 
#Number of SNP type I probes: 25 
#Number of SNP type II probes: 40 
head(getProbeInfo(manifest))

MSet <- preprocessRaw(RGset)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
beta <- getBeta(GRset)
sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
annotation <- getAnnotation(GRset)
names(annotation)

MSet.illumina <- preprocessIllumina(RGset, bg.correct = TRUE,
                                    normalize = "controls")
MSet.swan <- preprocessSWAN(RGset)
GRset.quantile <- preprocessQuantile(RGset, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE, 
                                     mergeManifest = FALSE, sex = NULL)

#### approximately close to res[[3]], normalized input user data 485512*12
MSet.noob <- preprocessNoob(RGset)
GRset.funnorm <- preprocessFunnorm(RGset)


res5 <- projectCellType(getBeta(GRset.quantile)[rownames(FlowSorted.Blood.450k.JaffeModelPars),], 
                        as.matrix(FlowSorted.Blood.450k.JaffeModelPars))
