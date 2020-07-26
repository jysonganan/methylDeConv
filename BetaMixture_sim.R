



#https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

estBetaParams_datRow <- function(datRow){
  mu = mean(datRow, na.rm = TRUE)
  var = var(datRow, na.rm = TRUE)
  return(estBetaParams(mu, var))
}

betaSim <- function(ref_betamatrix, ref_phenotype){
  phenotype <- levels(as.factor(ref_phenotype))
  dat_sim <- matrix(NA, nrow(ref_betamatrix), length(phenotype))
  for (i in 1:length(phenotype)){
    dat <- ref_betamatrix[,ref_phenotype == phenotype[i]]
    dat_sim[,i] <- apply(dat, 1, function(x){
      para <- estBetaParams_datRow(x)
      return(rbeta(1, para[[1]], para[[2]]))
    })
  }
  rownames(dat_sim) <- rownames(ref_betamatrix)
  colnames(dat_sim) <- phenotype
  return(dat_sim)
}


sim1 <- betaSim(ref_betamatrix, ref_phenotype)

set.seed(2)
purified_datasets_sim <- list()
for (i in 1:20){
  purified_datasets_sim[[i]] <-  betaSim(ref_betamatrix, ref_phenotype)
}

save(purified_datasets_sim, file = "purified_datasets_sim.RData")
## impute NAs with zeros

## Within in each of three groups (no epithelial, low epithelial and high epithelial)
## generate 600 mixture samples:

## combined the 20 ramdomly generated purified datasets (beta distribution) 
## with 30 randomly generated probablities to generate 600 mixture samples.



## 1. no epithelial
set.seed(3)
proportions_sim <- rdirichlet(30, c(1,1,1,1,1,1))
proportions_sim <- cbind(proportions_sim, 0)
true_proportions_sim <- proportions_sim

for (i in 1:19){
  true_proportions_sim <- rbind(true_proportions_sim, proportions_sim)
}

mixture_sim <- list()
for (i in 1:20){
  for (j in 1:30){
    res = true_proportions_sim[j,] * t(purified_datasets_sim[[i]])
    res = apply(res, 2, sum)
    mixture_sim <- c(mixture_sim, res)
  }
}

mixture_sim_mat <- matrix(NA, nrow(ref_betamatrix), 600)
for (i in 1:600){
  mixture_sim_mat[,i] <- mixture_sim[[i]]
}

true_proportions_sim_noEpithelial <- true_proportions_sim
mixture_sim_mat_noEpithelial <- mixture_sim_mat



## 2. low epithelial
set.seed(4)
proportions_sim <- rdirichlet(30, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(30, 0.1, 0.2)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
true_proportions_sim <- proportions_sim

for (i in 1:19){
  true_proportions_sim <- rbind(true_proportions_sim, proportions_sim)
}

mixture_sim <- list()
for (i in 1:20){
  for (j in 1:30){
    res = true_proportions_sim[j,] * t(purified_datasets_sim[[i]])
    res = apply(res, 2, sum)
    mixture_sim <- c(mixture_sim, res)
  }
}

mixture_sim_mat <- matrix(NA, nrow(ref_betamatrix), 600)
for (i in 1:600){
  mixture_sim_mat[,i] <- mixture_sim[[i]]
}

true_proportions_sim_lowEpithelial <- true_proportions_sim
mixture_sim_mat_lowEpithelial <- mixture_sim_mat

## 3. high epithelial

set.seed(5)
proportions_sim <- rdirichlet(30, c(1,1,1,1,1,1))
proportions_sim_epithelial <- runif(30, 0.8, 1.0)
proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
true_proportions_sim <- proportions_sim

for (i in 1:19){
  true_proportions_sim <- rbind(true_proportions_sim, proportions_sim)
}

mixture_sim <- list()
for (i in 1:20){
  for (j in 1:30){
    res = true_proportions_sim[j,] * t(purified_datasets_sim[[i]])
    res = apply(res, 2, sum)
    mixture_sim <- c(mixture_sim, res)
  }
}

mixture_sim_mat <- matrix(NA, nrow(ref_betamatrix), 600)
for (i in 1:600){
  mixture_sim_mat[,i] <- mixture_sim[[i]]
}

true_proportions_sim_highEpithelial <- true_proportions_sim
mixture_sim_mat_highEpithelial <- mixture_sim_mat
