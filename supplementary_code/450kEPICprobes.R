## for each cell type, test statistics for one-vs-all t test on EPIC & 450k arrays (common probes)
library(methylDeConv)
reference_EPIC <- build_reference_EPIC(extend = FALSE)
compTable_EPIC <- ref_compTable(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype)
compTable_EPIC <- compTable_EPIC[,3:8]

### 450k
reference_450k <- build_reference_450k(extend = FALSE)
compTable_450k <- ref_compTable(reference_450k$ref_betamatrix, reference_450k$ref_phenotype)
compTable_450k <- compTable_450k[,3:8]


splitit <- function(x) {split(seq_along(x), x)}

ref_phenotype <- as.factor(reference_EPIC$ref_phenotype)
tIndexes <- splitit(ref_phenotype)

tstatList_EPIC <- lapply(tIndexes, function(i) {
  x <- rep(0,ncol(reference_EPIC$ref_betamatrix))
  x[i] <- 1
  return(genefilter::rowttests(reference_EPIC$ref_betamatrix, factor(x)))
})


ref_phenotype <- as.factor(reference_450k$ref_phenotype)
tIndexes <- splitit(ref_phenotype)

tstatList_450k <- lapply(tIndexes, function(i) {
  x <- rep(0,ncol(reference_450k$ref_betamatrix))
  x[i] <- 1
  return(genefilter::rowttests(reference_450k$ref_betamatrix, factor(x)))
})


common_probes = intersect(rownames(compTable_450k),rownames(compTable_EPIC))
## B cell vs all else
Bcell_stat <- cbind(tstatList_450k[[1]][common_probes, c(2,3)], tstatList_EPIC[[1]][common_probes, c(2,3)])

CD4T_stat <- cbind(tstatList_450k[[2]][common_probes, c(2,3)], tstatList_EPIC[[2]][common_probes, c(2,3)])

CD8T_stat <- cbind(tstatList_450k[[3]][common_probes, c(2,3)], tstatList_EPIC[[3]][common_probes, c(2,3)])

Mono_stat <- cbind(tstatList_450k[[4]][common_probes, c(2,3)], tstatList_EPIC[[5]][common_probes, c(2,3)])

Neu_stat <- cbind(tstatList_450k[[5]][common_probes, c(2,3)], tstatList_EPIC[[4]][common_probes, c(2,3)])

NK_stat <- cbind(tstatList_450k[[6]][common_probes, c(2,3)], tstatList_EPIC[[6]][common_probes, c(2,3)])

save("Bcell_stat", "CD4T_stat", "CD8T_stat", "Mono_stat", "Neu_stat", "NK_stat", file = "450kEPIC_commonProbes_OneVsAllttest_statistics.RData")




## 2. hypergeometric test 


probe_1_EPIC <- ref_probe_selection_oneVsAllttest(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype, MaxDMRs = 100)
probe_2_EPIC <- ref_probe_selection_oneVsAllttest(reference_EPIC$ref_betamatrix, reference_EPIC$ref_phenotype, MaxDMRs = 300)

probe_1_450k <- ref_probe_selection_oneVsAllttest(reference_450k$ref_betamatrix, reference_450k$ref_phenotype, MaxDMRs = 100)
probe_2_450k <- ref_probe_selection_oneVsAllttest(reference_450k$ref_betamatrix, reference_450k$ref_phenotype, MaxDMRs = 300)


phyper(length(intersect(common_probes, intersect(probe_1_450k, probe_1_EPIC))),
       length(intersect(common_probes, probe_1_EPIC)),
       length(common_probes) - length(intersect(common_probes, probe_1_EPIC)),
       length(intersect(common_probes, probe_1_450k)), lower.tail = FALSE,  log.p = FALSE)
       #length(intersect(common_probes, intersect(probe_1_450k, probe_1_EPIC))),
       #length(common_probes[(common_probes%in%probe_1_450k)&(!common_probes%in%probe_1_EPIC)]),
       #length(common_probes[(!common_probes%in%probe_1_450k)&(common_probes%in%probe_1_EPIC)]),
       #length(common_probes[(!common_probes%in%probe_1_450k)&(!common_probes%in%probe_1_EPIC)])
       
phyper(length(intersect(common_probes, intersect(probe_2_450k, probe_2_EPIC))),
       length(intersect(common_probes, probe_2_EPIC)),
       length(common_probes) - length(intersect(common_probes, probe_2_EPIC)),
       length(intersect(common_probes, probe_2_450k)), lower.tail = FALSE,  log.p = FALSE)
