
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


#'Purified profiles sampling in Beta Mixture Simulation.
#'
#'Purified profiles sampling in Beta Mixture Simulation.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@param n The number of simulated purfied cell profiles. Default value: 20.
#'@return The simulated purified profiles.
#'@export

betaSim_purifiedProfiles <- function(ref_betamatrix, ref_phenotype, n = 20){
  set.seed(2)
  purified_datasets_sim <- list()
  for (i in 1:n){
    purified_datasets_sim[[i]] <-  betaSim(ref_betamatrix, ref_phenotype)
  }

  for (i in 1:n){
    purified_datasets_sim[[i]][which(is.na(purified_datasets_sim[[i]])==TRUE)] <- 0
  }
  return(purified_datasets_sim)
}



#'Beta Mixture Simulation
#'
#'Beta Mixture Simulation to generate the mixture profiles with true/known proportions.
#'@param nonimmune_level The levels are 1,2,3,4,5. Default value is 1. 1: No non-immune component; level 2: the non-immune proportion is 0.1-0.2;
#'level 3: the non-immune proportion is 0.2-0.5; level 4: the non-immune proportion is 0.5-0.8; level 5: the non-immune proportion is 0.8-0.9.
#'@param purified_datasets_sim The reference matrix ref_betamatrix.
#'@param n The number of simulated proportions. Default value: 30.
#'@return The simulated purified profiles.
#'@export

betaSim_mixtureProfiles <- function(nonimmune_level = 1, purified_datasets_sim, n = 30){
  n_purified = length(purified_datasets_sim)
  if (nonimmune_level == 1){
    set.seed(3)
    proportions_sim <- MCMCpack::rdirichlet(n, c(1,1,1,1,1,1))
    proportions_sim <- cbind(proportions_sim, 0)
    proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]
    true_proportions_sim <- proportions_sim
      }
  else if (nonimmune_level == 2){
    set.seed(4)
    proportions_sim <- MCMCpack::rdirichlet(n, c(1,1,1,1,1,1))
    proportions_sim_epithelial <- stats::runif(n, 0.1, 0.2)
    proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
    proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
    proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]
    true_proportions_sim <- proportions_sim
  }
  else if (nonimmune_level == 3){
    set.seed(11)
    proportions_sim <- rdirichlet(n, c(1,1,1,1,1,1))
    proportions_sim_epithelial <- runif(n, 0.2, 0.5)
    proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
    proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
    proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]
    true_proportions_sim <- proportions_sim
  }
  else if (nonimmune_level == 4){
    set.seed(7)
    proportions_sim <- rdirichlet(n, c(1,1,1,1,1,1))
    proportions_sim_epithelial <- runif(n, 0.5, 0.8)
    proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
    proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
    proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]
    true_proportions_sim <- proportions_sim
  }
  else if (nonimmune_level == 5){
    set.seed(5)
    proportions_sim <- rdirichlet(n, c(1,1,1,1,1,1))
    proportions_sim_epithelial <- runif(n, 0.8, 0.9)
    proportions_sim <- proportions_sim * (1-proportions_sim_epithelial)
    proportions_sim <- cbind(proportions_sim, proportions_sim_epithelial)
    proportions_sim <- proportions_sim[,c(1,2,3,7,4,5,6)]
    true_proportions_sim <- proportions_sim
  }


  for (i in 1:(n_purified-1)){
    true_proportions_sim <- rbind(true_proportions_sim, proportions_sim)
  }

  mixture_sim_mat <- NULL
  for (i in 1:n_purified){
    res <- matrix(NA, nrow(purified_datasets_sim[[i]]), n)
    for (m in 1:nrow(res)){
      res[m,] <- purified_datasets_sim[[i]][m,] %*% t(proportions_sim)
    }

    mixture_sim_mat <- cbind(mixture_sim_mat, res)
  }

  rownames(mixture_sim_mat) <- rownames(purified_datasets_sim[[1]])
  colnames(true_proportions_sim) <- colnames(purified_datasets_sim[[1]])
  return(list(mixture_sim_mat = mixture_sim_mat, true_proportions_sim = true_proportions_sim))
}

