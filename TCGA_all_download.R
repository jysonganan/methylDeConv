library(RTCGA)

disease_all <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","COADREAD","DLBC","ESCA","GBM",
                 "GBMLGG","HNSC","KICH","KIPAN","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC",
                 "MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","STES","TGCT","THCA",
                 "THYM","UCEC","UCS","UVM")

# Download methylation data
for (i in 1:length(disease_all)){
  disease <- disease_all[i]
  downloadTCGA(cancerTypes = disease, dataSet = "Merge_methylation__humanmethylation450", destDir = "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_methyl", date = "2016-01-28")
  folder_name <-paste0("gdac.broadinstitute.org_", disease, ".Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0")
  folder_name <- paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_methyl/", folder_name)
  sub_folder_name <- paste0(folder_name,"/MethylDataset")
  system(paste("mkdir",sub_folder_name))
  system(paste("mv",paste0(folder_name, "/",disease,"*.txt"),sub_folder_name))
  
  path_disease <- list.files(path = sub_folder_name, full.names = TRUE, recursive = TRUE)
  BetaMatrix <- readTCGA(path_disease, dataType = "methylation")
  system(paste("rm -rf", folder_name))
  rownames(BetaMatrix) <- BetaMatrix[,1]
  BetaMatrix <- BetaMatrix[,-1]
  BetaMatrix <- t(BetaMatrix)
  save(BetaMatrix, file = paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_methyl/",disease,"_methyl.RData"))
  rm(BetaMatrix)
}

# Download RNAseq data
for (i in 1:length(disease_all)){
  disease <- disease_all[i]
  downloadTCGA(cancerTypes = disease, dataSet = "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level", destDir = "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_rna", date = "2016-01-28")
  folder_name <- paste0("gdac.broadinstitute.org_",disease,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0")
  folder_name <- paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_rna/", folder_name)
  sub_folder_name <- paste0(folder_name,"/RNADataset")
  system(paste("mkdir",sub_folder_name))
  system(paste("mv",paste0(folder_name, "/",disease,"*.txt"),sub_folder_name))
  
  path_disease <- list.files(path = sub_folder_name, full.names = TRUE, recursive = TRUE)
  RNAMatrix <- readTCGA(path_disease, dataType = "rnaseq")
  rownames(RNAMatrix) <- RNAMatrix[,1]
  RNAMatrix <- RNAMatrix[,-1]
  RNAMatrix <- t(RNAMatrix)
  save(RNAMatrix, file = paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_rna/",disease,"_rna.RData"))
}


# Download clinical data
for (i in 1:length(disease_all)){
  disease <- disease_all[i]
  downloadTCGA(cancerTypes = disease, destDir = "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_clinical", date = "2016-01-28")
}



# Download mRNA data
disease_all <- c("BRCA","COAD","COADREAD","GBM","GBMLGG","KIPAN","KIRC","KIRP","LGG","LUAD","LUSC",
                 "OV","READ","UCEC")

for (i in 1:length(disease_all)){
  disease <- disease_all[i]
  downloadTCGA(cancerTypes = disease, dataSet = "Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.Level_3", destDir = "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mrna", date = "2016-01-28")
  folder_name <- paste0("gdac.broadinstitute.org_",disease,".Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.Level_3.2016012800.0.0")
  folder_name <- paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mrna/", folder_name)
  sub_folder_name <- paste0(folder_name,"/mRNADataset")
  system(paste("mkdir",sub_folder_name))
  system(paste("mv",paste0(folder_name, "/",disease,"*.txt"),sub_folder_name))
  
  path_disease <- list.files(path = sub_folder_name, full.names = TRUE, recursive = TRUE)
  RNAMatrix <- readTCGA(path_disease, dataType = "mRNA")
  rownames(RNAMatrix) <- RNAMatrix[,1]
  RNAMatrix <- RNAMatrix[,-1]
  RNAMatrix <- t(RNAMatrix)
  save(RNAMatrix, file = paste0("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mrna/",disease,"_mrna.RData"))
}

#### instead mRNA
library(magrittr)
library(dplyr)
library(devtools)
library(RTCGA)
(cohorts <- infoTCGA() %>% 
   rownames() %>% 
   sub("-counts", "", x=.))


# dir.create( "data2" )
releaseDate <- "2016-01-28"
sapply( cohorts, function(element){
tryCatch({
downloadTCGA( cancerTypes = element, 
              dataSet = "Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.Level_3",
              destDir = "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mRNA", 
              date = releaseDate )},
error = function(cond){
   cat("Error: Maybe there weren't mutations data for ", element, " cancer.\n")
}
)
})

list.files( "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mRNA") %>% 
   file.path( "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mRNA", .) %>%
   file.rename( to = substr(.,start=1,stop=50))

list.files( "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mRNA") %>%
   file.path( "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mRNA", .) %>%
   sapply(function(x){
      if (x == "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mRNA/NA")
         file.remove(x)      
   })

list.files( "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mRNA") %>% 
   file.path( "/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mRNA", .) %>%
   sapply(function(x){
      file.path(x, list.files(x)) %>%
         grep(pattern = "MANIFEST.txt", x = ., value=TRUE) %>%
         file.remove()
      })

list.files("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mRNA") %>%
   file.path("/sonas-hs/krasnitz/hpc/data/pfproj/tcga_data/tcga_mRNA", .) %>%
   sapply(function(y){
      file.path(y, list.files(y)) %>%
      assign( value = .,
              x = paste0(list.files(y) %>%
                            gsub(x = .,
                                 pattern = "\\..*",
                                 replacement = "") %>%
                            gsub(x=., pattern="-", replacement = "_"),
                         ".mRNA.path"),
              envir = .GlobalEnv)
   })

ls() %>%
   grep("mRNA\\.path", x = ., value = TRUE) %>% 
   sapply(function(element){
      tryCatch({
         readTCGA(get(element, envir = .GlobalEnv),
               dataType = "mRNA") %>%
         assign(value = .,
                x = sub("\\.path", "", x = element),
                envir = .GlobalEnv )
      }, error = function(cond){
         cat(element)
      }) 
     invisible(NULL)
    }    
)

grep( "mRNA", ls(), value = TRUE) %>%
   grep("path", x=., value = TRUE, invert = TRUE) %>%
   cat( sep="," ) #can one to id better? as from use_data documentation:
   # ...    Unquoted names of existing objects to save
   devtools::use_data(BRCA.mRNA,COAD.mRNA,COADREAD.mRNA,GBMLGG.mRNA,
                      KIPAN.mRNA,KIRC.mRNA,KIRP.mRNA,LGG.mRNA,LUAD.mRNA,
                      LUSC.mRNA,OV.mRNA,READ.mRNA,UCEC.mRNA,
                      overwrite = TRUE,
                      compress="xz")

