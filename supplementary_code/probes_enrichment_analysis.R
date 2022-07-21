## enrichment analysis
library(methylDeConv)
#library(methylGSA)
library(missMethyl)

reference_EPIC <- build_reference_EPIC(extend = FALSE)
compTable_EPIC <- ref_compTable(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype)
compTable_EPIC <- compTable_EPIC[,3:8]
probes_EPIC_oneVsAllttest <- ref_probe_selection_oneVsAllttest(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype) ##600
set.seed(2)
probes_EPIC_glmnetpreselect <- ref_probe_selection_twoStage(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype,
                                        preselect = 300, ml_model = "elastic net")
probes_EPIC_glmnetpreselect <- probes_EPIC_glmnetpreselect[-1] ## 946
probes_EPIC_overlapping <- intersect(probes_EPIC_oneVsAllttest, probes_EPIC_glmnetpreselect) ## 347


reference_Extended <- build_reference_EPIC(extend = TRUE, include = "Epithelial")
compTable_Extended <- ref_compTable(reference_Extended$ref_betamatrix, reference_Extended$ref_phenotype)
compTable_Extended <- compTable_Extended[,3:9]
probes_Extended_oneVsAllttest <- ref_probe_selection_oneVsAllttest(reference_Extended$ref_betamatrix, reference_Extended$ref_phenotype) ##700
set.seed(2)
probes_Extended_glmnetpreselect <- ref_probe_selection_twoStage(reference_Extended$ref_betamatrix, reference_Extended$ref_phenotype,
                                                            preselect = 300, ml_model = "elastic net")
probes_Extended_glmnetpreselect <- probes_Extended_glmnetpreselect[-1] ## 1048
probes_Extended_overlapping <- intersect(probes_Extended_oneVsAllttest, probes_Extended_glmnetpreselect) ## 389


#### missMethyl, hypergeometric test (taking into account the number of CpG sites per gene on the 450K array)
gst_EPIC_oneVsAllttest_KEGG <- gometh(sig.cpg=probes_EPIC_oneVsAllttest, all.cpg=rownames(compTable_EPIC), collection="KEGG", array.type="EPIC")
gst_EPIC_oneVsAllttest_GO <- gometh(sig.cpg=probes_EPIC_oneVsAllttest, all.cpg=rownames(compTable_EPIC), collection="GO", array.type="EPIC")
gst_EPIC_glmnetpreselect_KEGG <- gometh(sig.cpg=probes_EPIC_glmnetpreselect, all.cpg=rownames(compTable_EPIC), collection="KEGG", array.type="EPIC")
gst_EPIC_glmnetpreselect_GO <- gometh(sig.cpg=probes_EPIC_glmnetpreselect, all.cpg=rownames(compTable_EPIC), collection="GO", array.type="EPIC")
gst_EPIC_overlapping_KEGG <- gometh(sig.cpg=probes_EPIC_overlapping, all.cpg=rownames(compTable_EPIC), collection="KEGG", array.type="EPIC")
gst_EPIC_overlapping_GO <- gometh(sig.cpg=probes_EPIC_overlapping, all.cpg=rownames(compTable_EPIC), collection="GO", array.type="EPIC")


gst_Extended_oneVsAllttest_KEGG <- gometh(sig.cpg=probes_Extended_oneVsAllttest, all.cpg=rownames(compTable_EPIC), collection="KEGG", array.type="EPIC")
gst_Extended_oneVsAllttest_GO <- gometh(sig.cpg=probes_Extended_oneVsAllttest, all.cpg=rownames(compTable_EPIC), collection="GO", array.type="EPIC")
gst_Extended_glmnetpreselect_KEGG <- gometh(sig.cpg=probes_Extended_glmnetpreselect, all.cpg=rownames(compTable_EPIC), collection="KEGG", array.type="EPIC")
gst_Extended_glmnetpreselect_GO <- gometh(sig.cpg=probes_Extended_glmnetpreselect,all.cpg=rownames(compTable_EPIC), collection="GO", array.type="EPIC")
gst_Extended_overlapping_KEGG <- gometh(sig.cpg=probes_Extended_overlapping, all.cpg=rownames(compTable_EPIC), collection="KEGG", array.type="EPIC")
gst_Extended_overlapping_GO <- gometh(sig.cpg=probes_Extended_overlapping, all.cpg=rownames(compTable_EPIC), collection="GO", array.type="EPIC")


gst_EPIC_oneVsAllttest_KEGG <- gst_EPIC_oneVsAllttest_KEGG[order(gst_EPIC_oneVsAllttest_KEGG[,"FDR"]),]
gst_EPIC_oneVsAllttest_GO <- gst_EPIC_oneVsAllttest_GO[order(gst_EPIC_oneVsAllttest_GO[,"FDR"]),]
gst_EPIC_glmnetpreselect_KEGG <- gst_EPIC_glmnetpreselect_KEGG[order(gst_EPIC_glmnetpreselect_KEGG[,"FDR"]),]
gst_EPIC_glmnetpreselect_GO <- gst_EPIC_glmnetpreselect_GO[order(gst_EPIC_glmnetpreselect_GO[,"FDR"]),]
gst_EPIC_overlapping_KEGG <- gst_EPIC_overlapping_KEGG[order(gst_EPIC_overlapping_KEGG[,"FDR"]),]
gst_EPIC_overlapping_GO <- gst_EPIC_overlapping_GO[order(gst_EPIC_overlapping_GO[,"FDR"]),]


gst_Extended_oneVsAllttest_KEGG <- gst_Extended_oneVsAllttest_KEGG[order(gst_Extended_oneVsAllttest_KEGG[,"FDR"]),]
gst_Extended_oneVsAllttest_GO <- gst_Extended_oneVsAllttest_GO[order(gst_Extended_oneVsAllttest_GO[,"FDR"]),]
gst_Extended_glmnetpreselect_KEGG <- gst_Extended_glmnetpreselect_KEGG[order(gst_Extended_glmnetpreselect_KEGG[,"FDR"]),]
gst_Extended_glmnetpreselect_GO <- gst_Extended_glmnetpreselect_GO[order(gst_Extended_glmnetpreselect_GO[,"FDR"]),]
gst_Extended_overlapping_KEGG <- gst_Extended_overlapping_KEGG[order(gst_Extended_overlapping_KEGG[,"FDR"]),]
gst_Extended_overlapping_GO <- gst_Extended_overlapping_GO[order(gst_Extended_overlapping_GO[,"FDR"]),]

# save("gst_EPIC_oneVsAllttest_KEGG", "gst_EPIC_oneVsAllttest_GO", "gst_EPIC_glmnetpreselect_KEGG", "gst_EPIC_glmnetpreselect_GO",
#      "gst_EPIC_overlapping_KEGG", "gst_EPIC_overlapping_GO", "gst_Extended_oneVsAllttest_KEGG", "gst_Extended_oneVsAllttest_GO",
#      "gst_Extended_glmnetpreselect_KEGG", "gst_Extended_glmnetpreselect_GO", "gst_Extended_overlapping_KEGG", "gst_Extended_overlapping_GO",
#      file = "enrichment_analysis_missMethyl.RData")




###
## use package methylGSA

library(methylGSA)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

cpg.mypval <- rep(0.5, nrow(compTable_EPIC))
names(cpg.mypval) <- rownames(compTable_EPIC)
cpg.mypval[which(names(cpg.mypval) %in% probes_EPIC_oneVsAllttest)] <- 0.0005
EPIC_oneVsAllttest_KEGG = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "KEGG", array.type = "EPIC", minsize = 15, maxsize = 500)
EPIC_oneVsAllttest_GO = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "GO", array.type = "EPIC", minsize = 15, maxsize = 500)

#cpg.mypval_1 <- runif(nrow(compTable_EPIC), 0.1,0.9)
#names(cpg.mypval_1) <- rownames(compTable_EPIC)
#cpg.mypval_1[which(names(cpg.mypval_1) %in% probes_EPIC_oneVsAllttest)] <- runif(600, 0.0001, 0.0005)


cpg.mypval <- rep(0.5, nrow(compTable_EPIC))
names(cpg.mypval) <- rownames(compTable_EPIC)
cpg.mypval[which(names(cpg.mypval) %in% probes_EPIC_glmnetpreselect)] <- 0.0005
EPIC_glmnetpreselect_KEGG = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "KEGG", array.type = "EPIC", minsize = 15, maxsize = 500)
EPIC_glmnetpreselect_GO = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "GO", array.type = "EPIC", minsize = 15, maxsize = 500)


cpg.mypval <- rep(0.5, nrow(compTable_EPIC))
names(cpg.mypval) <- rownames(compTable_EPIC)
cpg.mypval[which(names(cpg.mypval) %in% probes_EPIC_overlapping)] <- 0.0005
EPIC_overlapping_KEGG = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "KEGG", array.type = "EPIC", minsize = 15, maxsize = 500)
EPIC_overlapping_GO = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "GO", array.type = "EPIC", minsize = 15, maxsize = 500)






cpg.mypval <- rep(0.5, nrow(compTable_Extended))
names(cpg.mypval) <- rownames(compTable_Extended)
cpg.mypval[which(names(cpg.mypval) %in% probes_Extended_oneVsAllttest)] <- 0.0005
Extended_oneVsAllttest_KEGG = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "KEGG", array.type = "EPIC", minsize = 15, maxsize = 500)
Extended_oneVsAllttest_GO = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "GO", array.type = "EPIC", minsize = 15, maxsize = 500)


cpg.mypval <- rep(0.5, nrow(compTable_Extended))
names(cpg.mypval) <- rownames(compTable_Extended)
cpg.mypval[which(names(cpg.mypval) %in% probes_Extended_glmnetpreselect)] <- 0.0005
Extended_glmnetpreselect_KEGG = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "KEGG", array.type = "EPIC", minsize = 15, maxsize = 500)
Extended_glmnetpreselect_GO = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "GO", array.type = "EPIC", minsize = 15, maxsize = 500)


cpg.mypval <- rep(0.5, nrow(compTable_Extended))
names(cpg.mypval) <- rownames(compTable_Extended)
cpg.mypval[which(names(cpg.mypval) %in% probes_Extended_overlapping)] <- 0.0005
Extended_overlapping_KEGG = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "KEGG", array.type = "EPIC", minsize = 15, maxsize = 500)
Extended_overlapping_GO = methylgometh(cpg.pval = cpg.mypval, sig.cut = 0.001, topDE = NULL, GS.type = "GO", array.type = "EPIC", minsize = 15, maxsize = 500)
