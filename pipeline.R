## pipeline
# sd1=apply(y1, 1, sd, na.rm = F)
# include = which(sd1>=quantile(sd1,0.05))
# O = y1[include,]
# edata1 = edata[include,]

preprocess <- function(input_methyl, input_phenotype, input_covariate = NULL){
  
}
## Generate reference profiles
MethylDeconv <- function(input_methyl, input_phenotype, input_covariate = NULL, input_reference = NULL,method = "RUV"){
  
  x1 <- as.numeric(input_phenotype[,2])
  
  if(!is.null(input_covariate)){
    x2list <- list()
    for (i in 1:ncol(input_covariate)){
      x2list[[i]] <- as.numeric(input_covariate[,i])
    }
  }
  
  if(method == "RefFreeEWAS"){
    library(RefFreeEWAS)
    tmp1 <- lm(t(input_methyl)~x1)
    dim3 <- EstDimRMT(cbind(t(coef(tmp1)),t(resid(tmp1))), FALSE)$dim
    if(is.null(input_covariate)){
      test1 <- RefFreeEwasModel(input_methyl, cbind(1,x1), dim3)
    }
    else{
      tmp2 <- cbind(1,x1)
      for (i in 1:length(x2list)){
        tmp2 <- cbind(tmp2, x2list[[i]])
      }
      test1 <- RefFreeEwasModel(input_methyl,tmp2,dim3) 
    }
    testBoot1 = BootRefFreeEwasModel(test1, 50)
  }
  
  if(method == "SVA"){
    lib = c("MASS","sva","limma")
    lapply(lib, require, character.only = TRUE)
    mod1 <- model.matrix(~.,cbind(input_phenotype, input_covariate))
    mod2 <- model.matrix(~.,input_covariate)
    svobj1 <- sva(input_methyl,mod1,mod2, n.sv = NULL, method = "two-step")
    modSv1 <- cbind(mod1, svobj1$sv)
    fit1 <- lmFit(input_methyl, modSv1, method = "robust")
    fite1 <- eBayes(fit1)
    #formdf1 <- as.formula(paste("~", paste(names(df1)[2:3],collapse="+")))
  }
  
  if(method == "ReFACTor"){  #####
    ## def include before preprocess
    orign_methyl <- orign_methyl[include,]
    input_methyl <- input_methyl[include,]
   
    for (site in 1:nrow(orign_methyl))
    {
      model <- lm(orign_methyl[site,]~x1+x2+x3+x4+x5+x6)
      orign_methyl_adj[site,] = residuals(model)
    }
    orign_methyl = orign_methyl_adj
    print('Running a standard PCA...')
    pcs = prcomp(scale(t(orign_methyl)))
    coeff = pcs$rotation
    score = pcs$x
    print('Compute a low rank approximation of input data and rank sites...')
    x = score[,1:7]%*%t(coeff[,1:7])
    An = scale(t(orign_methyl), center = T, scale = F)
    Bn = scale(x, center = T, scale = F)
    An = t(t(An)*(1/sqrt(apply(An^2,2,sum))))
    Bn = t(t(Bn)*(1/sqrt(apply(Bn^2,2,sum))))
    
    #find the distance of each site from its low rank approximation
    distances = apply((An-Bn)^2,2,sum)^0.5
    dsort = sort(distances, index.return = T)
    ranked_list = dsort$ix
    
    print('Compute ReFACTor components...')
    sites = ranked_list[1:500]
    pcs = prcomp(scale(t(orign_methyl[sites,])))
    first_score <- score[,1:7]
    score = pcs$x
    mod = model.matrix(~x1+x2+x3+x4+x5+x6+first_score)
    fit1 <- lmFit(input_methyl, mod, method = "ls")
    fite1 <- eBayes(fit1)
  }
  
  if(method == "RefFreeCellMix"){
    library(RefFreeEWAS)
    cell <- RefFreeCellMix(y2, mu0 = NULL, K = 7, iters = 5, Yfinal = NULL, verbose = TRUE)
    mod1 <- model.matrix(~x1+x2+x3+x4+x5+x6+cell$Omega)
    fit1 <- lmFit(y2, mod1, method = "robust")
    fite1 <- eBayes(fit1)
  }
  
  
  if(method == "HousemanRef_based"){
    library(RefFreeEWAS)
    if(is.null(input_reference)){
      print('Reference matrix set as whole blood by default!')
    }
    else{
      #1. estimate cell proportions
      w = projectMix(input_methyl, input_reference)
      print("Cell proportions estimated!")
      #2. adjust for cell proportions
      mod1 <- model.matrix(~., cbind(input_phenotype, input_covariate, w))
      fit1 <- lmFit(input_methyl, mod1, method = "robust")
      fite1 <- eBayes(fit1)
      tab1  = topTable(fite1, coef = "x1",number=length(input_methyl[,1]), p.val=0.05,adjust = "fdr")
      write.table(tab1,file="Ewasher_data_houseman_CpG.csv",sep=",")
      return(w, fite1, tab1)
    }
  }
    
  #   beta=t(as.matrix(read.table(".../input_data.txt",sep="\t",header=T)))
  #   source("Rcodes_Cell_mixture.R")
  #   #This program will create a file by the name “Ewasher_Houseman_cell.csv”, 
  #   # which will have cell proportions for the 204 subjects 
  #   library(limma)
  #   library(MASS)
  #   cell = read.csv("Ewasher_data_housecell.csv",header=T)
  #   cell1 = cell[,2]
  #   cell2 = cell[,3]
  #   cell3 = cell[,4]
  #   cell4 = cell[,5]
  #   cell5 = cell[,6]
  #   cell6 = cell[,7]
  #   
  #   mod11 = model.matrix(~x1+x2+x3+x4+x5+x6+cell2+cell3+cell4+cell5+cell6)
  #   fit1  = lmFit(y12,mod11,method="robust")
  #   fite1 = eBayes(fit1)
  #   tab1  = topTable(fite1, coef = "x1",number=length(y1[,1]), p.val=0.05,adjust = "fdr")
  #   write.table(tab1,file="Ewasher_data_houseman_CpG.csv",sep=",")
  # 
  # 
  
  
  
  
  if(method == "450k_Ref_based"){
    lib = c("minfi","quadprog","FlowSorted.Blood.450k",
            "IlluminaHumanMethylation450kmanifest",
            "IlluminaHumanMethylation450kanno.ilmn12.hg19");
    lapply(lib, require, character.only = TRUE);
    
    #input_methyl
    # grSet1=read.table("input_data.txt",header=T); grSet=grSet1[,-1];
    # rownames(grSet)=grSet1[,1]; grSet=data.matrix(grSet)
    # 
    if(is.null(input_reference)){
      referenceMset = get('FlowSorted.Blood.450k.compTable')
      cell = c("CD8T","CD4T", "NK","Bcell","Mono","Gran","Eos")
      compData = minfi:::pickCompProbes(referenceMset, cellTypes=cell)
      coefs = compData$coefEsts
      coefs = coefs[intersect(rownames(input_methyl),rownames(coefs)),]
      rm(referenceMset)
      counts = minfi:::projectCellType(input_methyl[rownames(coefs), ], coefs)
      rownames(counts) = colnames(input_methyl)
      retun(counts)
    }
  
    ## adjust for cell proportions
    # cell =  counts
    # cell1 = cell[,1];cell2 = cell[,3];cell3 = cell[,4];cell5 = cell[,6];
    # cell6 = cell[,7];cell7 = cell[,8];
    # mod11 = model.matrix(~x1+x2+x3+x4+x5+x6+cell2+cell3+cell4+cell5+cell6+cell7);
    # fit1  = lmFit(y12,mod11,method="robust");
    # fite1 = eBayes(fit1);
    # tab1  = topTable(fite1, coef = "x1",number=length(y1[,1]), p.val=0.05,adjust = "fdr");
    # write.table(tab1,file="Ewasher_data_minfi_CpG.csv",sep=",")
    # 
  }
  
  if(method == "RUV"){
    
  }
}

setwd("/Users/junesong/Desktop/causal inference/FastLmm.Py/doc/FastLmm-EWASher/demo")

y <- read.table("input_data.txt", sep = "\t", header = T)
y1 <- y[,-1]
rownames(y1) = y[,1]
miss = which(is.na(y1),arr.ind=TRUE)[,1]
y1 = y1[-miss,]
edata = as.matrix(log2(y1/(1-y1)))


covar1 = read.table("input_phenotype.txt",header=FALSE);
covar2<-read.table("covariates.txt",header=FALSE)
x1 = as.numeric(covar1[,2]);
x2<-as.numeric(covar2[,3])
x3<-as.numeric(covar2[,4])
x4<-as.numeric(covar2[,5])
x5<-as.numeric(covar2[,6])
x6<-as.numeric(covar2[,7])

tmp1 = lm(t(edata)~x1);
dim3 = EstDimRMT(cbind(t(coef(tmp1)), t(resid(tmp1))), FALSE)$dim;
test1 = RefFreeEwasModel(y2,cbind(1,x1,x2,x3,x4,x5,x6),dim3);
testBoot1 = BootRefFreeEwasModel(test1, 50)
              