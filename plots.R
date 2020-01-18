# venn plot of the shared probes of different reference profiles for blood

load("FlowSorted.Blood.EPIC.IDOLModelPars.RData")
EPIC <- rownames(FlowSorted.Blood.EPIC.IDOLModelPars)

library(EpiDISH)
data("centDHSbloodDMC.m")
EpiDISH <- rownames(centDHSbloodDMC.m)

library(FlowSorted.Blood.450k)
data("FlowSorted.Blood.450k.JaffeModelPars")
Jaffe450k <- rownames(FlowSorted.Blood.450k.JaffeModelPars)

library(VennDiagram)
grid.newpage()
draw.triple.venn(area1 = 450, area2 = 333, area3 = 600, n12 = 13, n23 = 110, n13 = 54, 
                 n123 = 7, category = c("EPIC", "EpiDISH", "Jaffe450k"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))

#cbind(centDHSbloodDMC.m[intersect(EpiDISH,Jaffe450k),], FlowSorted.Blood.450k.JaffeModelPars[intersect(EpiDISH,Jaffe450k),])[1:10,]
