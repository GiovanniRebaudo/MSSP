#functions description: 
#generate_zipf :  returns the pmfs of a number J of zipf distributions 
#                 with overlapping support
#generate_zipf_reorder :  returns the pmfs of a number J of zipf distributions 
#                 with overlapping support and "similar" p.m.f.
#sample_from_pop: returns a sample with replacement from one pop given pmfs
#sample_from_pop_all: returns J samples with replacement from J pop given pmfs
#sample_alpha_DP: returns a sample of the concentration parameter of a DP
#                 from its full conditional in a marginal sampler 
#prob_new_group_DP: ratio between EPPF of DP corresponding to a group 
#                   of equal obs with new value (needed for +SSM)
#prob_new_group_PY: ratio between EPPF of PY corresponding to a group 
#                   of equal obs with new value (needed for +SSM)
#prob_new_group_DM: ratio between EPPF of DM corresponding to a group 
#                   of equal obs with new value (needed for +SSM)
#prob_group_DP: ratio between EPPF of DP corresponding to a group 
#               of equal obs with old value (needed for +SSM)
#prob_group_PY: ratio between EPPF of PY corresponding to a group 
#               of equal obs with old value (needed for +SSM)
#prob_group_DM: ratio between EPPF of PY corresponding to a group (needed for +SSM)


generate_zipf<- function(param = c(rep(1.3, 2), rep(2, 6)), 
                         tot_species = 3000, j_species = 2500, seed = 0){
  ## The function returns the pmfs of a number J of zipf distributions
  ## J is implicitly defined by the length of the input param
  ##inputs: 
  ## param[j] = zipf parameter of pop j 
  ## tot_species = the total number of species jointly in the J populations
  ## j_species = the number of species in each population
  ## seed = different seed generate different distributions
  ##output: 
  ## truth = a list of pmfs (UP TO THE NORM COST), truth[[j]][c] = pr(X_{ji} = c)  
  
  set.seed(seed)
  J = length(param) 
  truth = vector("list", J)
  
  for(j in 1:length(param)){
    truth[[j]] = rep(0, tot_species)
    species = sample(tot_species, j_species)
    truth[[j]][species] = (1 / c(1:j_species))**(param[j])
  }
  
  return(truth)
}

generate_zipf_reorder<- function(param = c(rep(1.3, 2), rep(2, 6)), 
                                 tot_species = 3000, j_species = 2500, seed = 0){
  ## The function returns the pmfs of a number J of zipf distributions
  ## J is implicitly defined by the length of the input param
  ##inputs: 
  ## param[j] = zipf parameter of pop j 
  ## tot_species = the total number of species jointly in the J populations
  ## j_species = the number of species in each population
  ## seed = different seed generate different distributions
  ##output: 
  ## truth = a list of pmfs (UP TO THE NORM COST), truth[[j]][c] = pr(X_{ji} = c)  
  
  set.seed(seed)
  J = length(param) 
  truth = vector("list", J)
  
  for(j in 1:length(param)){
    truth[[j]] = rep(0, tot_species)
    species = sort(sample(tot_species, j_species))
    truth[[j]][species] = (1 / c(1:j_species))**(param[j])
  }
  
  return(truth)
}

sample_from_pop <- function(j, truth, size = 1){
  ## The function returns a sample with replacement 
  ## of size ``size" from population j 
  ##inputs: 
  ## j = label of the population to sample
  ## truth = a list of pmfs, truth[[j]][c] = pr(X_{ji} = c)  
  ## size = how many samples
  
  return(sample(length(truth[[j]]), size, replace = TRUE, prob = truth[[j]])) 
}

sample_from_pop_all <- function(truth, size = NULL, seed = 0, verbose = TRUE){
  ## The function returns a sample with replacement 
  ## of size ``size" from population j 
  ##inputs: 
  ## truth = a list of pmfs, truth[[j]][c] = pr(X_{ji} = c)  
  ## size = vector of how many samples (default is one obs for each pop)
  ## seed 
  set.seed(seed)
  J = length(truth)
  
  if(is.null(size)){
    size = rep(1, J)
  } else if(J!=length(size)){
    if(verbose){
      cat("Size length differs from num. of pops.")
      cat("\nGenerating samples of size", size[1], "from the", J, "pmfs.\n")
    }
    size = rep(size[1],J)
  }
  
  X = matrix(NA, nrow = J, ncol = max(size))
  for (j in 1:J){
    X[j,1:size[j]] = sample(length(truth[[j]]), size[j], 
                            replace = TRUE, prob = truth[[j]])
  }
  return(X) 
}


sample_alpha_DP <- function(alpha, K, n, a = 1, b = 1){
  ## sample concentration parameter of a DP from its full conditional 
  ## via variable augmentation
  ##inputs: 
  ## alpha = current value of the concentration parameter
  ## K = number of unique species
  ## n = customer in the restaurant
  ## a,b = hyperparameters of the gamma prior on the concentration param
  ##output: 
  ## a sample from the full conditional of the concentration param
  
  eta = rbeta(1, alpha, n)
  alpha = rgamma(1, shape = a + K, rate = b - log(eta))
  return(alpha)
}


prob_new_group_DP <- function(n, alpha, ntot, log = FALSE){
  ## returns the probability of observing a future group 
  ## of equal obs with new value in the DP
  ##inputs: 
  ##  n = dimension of the group
  ##  alpha = concentration parameter
  ##  ntot = number of observations observed up to this point
  ##  log = whether returning the log prob
  
  if(log){
    return((lfactorial(alpha + ntot - 1) + 
              log(alpha) + lfactorial(n - 1) - 
              lfactorial(alpha + ntot + n - 1)) )
  }else{
    return(exp(lfactorial(alpha + ntot - 1) + 
                 log(alpha) + lfactorial(n - 1) - 
                 lfactorial(alpha + ntot + n - 1)) ) 
  }
}


prob_new_group_PY <- function(n, alpha, sigma, ntot, D, log = FALSE){
  ## returns the probability of observing a future group 
  ## of equal obs with new value in the PY
  ##inputs: 
  ##  n = dimension of the group
  ##  alpha = concentration parameter
  ##  sigma = discount parameter
  ##  ntot = number of observations observed up to this point
  ##  D = number of unique observations observed up to this point
  ##  log = whether returning the log prob
  
  if(log){
    return((lgamma(alpha + ntot) - 
              lgamma(alpha + ntot + n) + 
              log(alpha + D * sigma) - lgamma(1 - sigma) +
              lgamma(n - sigma) ) )
  }else{
    return(exp(lgamma(alpha + ntot) - 
                 lgamma(alpha + ntot + n) + 
                 log(alpha + D * sigma) - lgamma(1 - sigma) +
                 lgamma(n - sigma) ) ) 
  }
}

prob_new_group_DM <- function(n, rho, M, ntot, D, log = FALSE){
  ## returns the probability of observing a future group 
  ## of equal obs with new value in the DM
  ##inputs: 
  ##  n = dimension of the group
  ##  rho = parameter
  ##  M = num of support points
  ##  ntot = number of observations observed up to this point
  ##  D = number of unique observations observed up to this point
  ##  log = whether returning the log prob
  
  if(log){
    return((lgamma(M*rho + ntot) - 
              lgamma(M*rho + ntot + n) + 
              log(M - D) - lgamma(rho) +
              lgamma(n +1) ) )
  }else{
    return(exp(lgamma(M*rho + ntot) - 
                 lgamma(M*rho + ntot + n) + 
                 log(M - D) - lgamma(rho) +
                 lgamma(n +1) ) ) 
  }
}

prob_group_DP <- function(n, nc, alpha, ntot, log = FALSE){
  ## returns the probability of observing a future group 
  ## of equal obs with old value in the DP
  ##inputs: 
  ##  n = dimension of the group
  ##  nc = number of obs with the same species already observed
  ##  alpha = concentration parameter
  ##  ntot = number of observations observed up to this point
  ##  log = whether returning the log prob
  
  if(log){
    return((lfactorial(nc + n - 1) -
              lfactorial(nc - 1) + 
              lfactorial(alpha + ntot - 1) - 
              lfactorial(alpha + ntot + n - 1)) )
  }else{
    return(exp(lfactorial(nc + n - 1) -
                 lfactorial(nc - 1) + 
                 lfactorial(alpha + ntot - 1) - 
                 lfactorial(alpha + ntot + n - 1)) ) 
  }
}

prob_group_PY <- function(n, nc, alpha, sigma, ntot, log = FALSE){
  ## returns the probability of observing a future group 
  ## of equal obs with old value in the PY
  ##inputs: 
  ##  n = dimension of the group
  ##  nc = number of obs with the same species already observed
  ##  alpha = concentration parameter
  ##  sigma = discount parameter
  ##  ntot = number of observations observed up to this point
  ##  log = whether returning the log prob
  
  if(log){
    return((lgamma(alpha + ntot) - 
              lgamma(alpha + ntot + n) - 
              lgamma(nc - sigma) +
              lgamma(nc + n - sigma) ) )
  }else{
    return((lgamma(alpha + ntot) - 
              lgamma(alpha + ntot + n) - 
              lgamma(nc - sigma) +
              lgamma(nc + n - sigma) ) ) 
  }
}

prob_group_DM <- function(n, nc, rho, M, ntot, log = FALSE){
  ## returns the probability of observing a future group 
  ## of equal obs with old value in the DM
  ##inputs: 
  ##  n = dimension of the group
  ##  nc = number of obs with the same species already observed
  ##  rho = parameter
  ##  M = num of support points
  ##  ntot = number of observations observed up to this point
  ##  log = whether returning the log prob
  
  if(log){
    return((lgamma(M*rho + ntot) - 
              lgamma(M*rho + ntot + n) - 
              lgamma(nc + rho) +
              lgamma(nc + n + rho) ) )
  }else{
    return(exp(lgamma(M*rho + ntot) - 
                 lgamma(M*rho + ntot + n) - 
                 lgamma(nc + rho) +
                 lgamma(nc + n + rho) ) ) 
  }
}