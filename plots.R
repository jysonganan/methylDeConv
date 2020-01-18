# venn plot of the shared probes of different reference profiles for blood

load("FlowSorted.Blood.EPIC.IDOLModelPars.RData")
EPIC <- rownames(FlowSorted.Blood.EPIC.IDOLModelPars)

library(EpiDISH)
data("centDHSbloodDMC.m")
EpiDISH <- rownames(centDHSbloodDMC.m)

library(FlowSorted.Blood.450k)
data("FlowSorted.Blood.450k.JaffeModelPars")
Jaffe450k <- rownames(FlowSorted.Blood.450k.JaffeModelPars)
