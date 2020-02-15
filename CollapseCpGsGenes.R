# collapse CpGs into genes
library(IlluminaHumanMethylation450kmanifest)
data(IlluminaHumanMethylation450kmanifest)
library(minfi)
getManifest(IlluminaHumanMethylation450kmanifest)
df <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = c("I", "II", "Control",
                              "I-Green", "I-Red", "SnpI", "SnpII"))
class(df)
getManifestInfo(IlluminaHumanMethylation450kmanifest, type = c("nLoci", "locusNames"))

head(getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "I"))
head(IlluminaHumanMethylation450kmanifest@data$TypeI)
head(IlluminaHumanMethylation450kmanifest@data$TypeII)
head(IlluminaHumanMethylation450kmanifest@data$TypeControl)
