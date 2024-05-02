#multiarmed bandit for species discovery via mSSP

#functions description: 
#generate_zipf :  returns the pmfs of a number J of zipf distributions 
#                 with overlapping support
#generate_zipf_reorder :  returns the pmfs of a number J of zipf distributions 
#                 with overlapping support and "similar" p.m.f.
#sample_from_pop: returns a sample with replacement from a given pmf
#sample_alpha_DP: returns a sample of the concentration parameter of a DP
#                 from its full conditional in a marginal sampler 
#prob_new_group_DP: ratio between EPPF of DP corresponding to a group 
#                   of equal obs with new value
#prob_new_group_PY: ratio between EPPF of PY corresponding to a group 
#                   of equal obs with new value
#prob_new_group_DM: ratio between EPPF of DM corresponding to a group 
#                   of equal obs with new value
#prob_group_DP: ratio between EPPF of DP corresponding to a group 
#               of equal obs with old value
#prob_group_PY: ratio between EPPF of PY corresponding to a group 
#               of equal obs with old value
#prob_group_DM: ratio between EPPF of PY corresponding to a group 
#               of equal obs with old value
#plusDP_MAB:  Multi arm bandit with the +DP, 
#             returns the cumulative number of species discovered
#uniform_MAB: Multi arm bandit with Uniform random choices,
#             returns the cumulative number of species discovered
#plusPY_MAB:  Multi arm bandit with the +DP, 
#             returns the cumulative number of species discovered
#indepDP_MAB: Multi arm bandit with independent DP, 
#             returns the cumulative number of species discovered
#Oracle_MAP: Multi arm bandit with oracle 

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

#+DP multi arm bandit problem 
plusDP_MAB<- function(data, a = 1, b = 1, p = rep(1/3, 3),
                       init_samples = 30, new_samples = 300, 
                       burnin = 10, iters = 200, seed = 0){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = list of J pmfs
  ##  a,b = hyperparameters of the gamma prior on the concentration params
  ##  p = hyperparameters of the mixture prior of the mixing proportions
  ##  init_samples = number of starting observations for the MAB
  ##  new_samples = number of sampling step of the MAB
  ##  burnin = length of butnin of each MCMC
  ##  iter = number of iter after burnin of each MCMC
  ##  seed 

  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = new_samples, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  species_discovered = rep(0, new_samples) #vector to save the num of discoveries

  tot_iter = burnin + iters #per each MCMC
  
  set.seed(seed)
  J = length(data) #tot number of populations
    
  X = matrix(NA,nrow = J, 
    ncol = init_samples + new_samples) #matrix of observations X[j,i]is X_{j,i}
    
  #sample initial observations
  for(j in 1:J){
    X[j, 1:init_samples] = sample_from_pop(j, data, size = init_samples) 
  }
    
  I = rep(init_samples, J) #initial sample sizes
    
  #initializing the hyperparam for the first run of the MAB
  alpha = rgamma(J + 1, shape = a, rate = b) 
  epsilon = rep(0.5, J)
  
  #initializing the restaurants for the first run of the MAB
  #r[iter,j,i] is r_{j,i}
  r = array(NA, dim = c(tot_iter+1, J, init_samples + new_samples) )
  
  for(j in 1:J){ #a
    unique_tab = tabulate(X[j,])
    for(x in which(unique_tab>0)){ #b
      if( x %in% X[-j,] ){ #c
        
        r[1, j, which(X[j,]==x)] = 0  
        
      }else{ #c
        
        r[1, j, which(X[j,]==x)] = rep(sample(c(0,1),1, 
                                  prob = c(epsilon[j], 1- epsilon[j]))
                                       ,length(which(X[j,]==x)))
        
      } #c
    } #b
  } #a
    
  for(newobs in 1:new_samples){ #MAB for cycle
    #MCMC for latent allocation to restaurants
    
    #vector for prob of new species for each mcmc iteration
    prob_new_species = matrix(0, nrow = tot_iter, ncol = J)
    
    #argument of pEPPF (i.e., matrix of counts)
    unique_tab = matrix(0, ncol = max(X, na.rm = TRUE) , nrow = J)
    for(j in 1:J){ #a
      unique_tab[j,] = tabulate(X[j,], nbins = max(X, na.rm = TRUE))
    } #a
    #customer counts in each restaurant
    #l0[j] = num of customers in common restaurant from pop j
    #lj[j] = num of customers in exclusive restraurant from pop j
    l0 = rep(0, J); lj = rep(0, J)
      
    for(iter in 2:(tot_iter+1)){
      #assigned customer to restaurants 
      r[iter,,] = r[iter-1,,]
        for(j in 1:J){ #a
          for(x in which(unique_tab[j,]>0)){ #b
            r[iter, j, which(X[j,]==x)] = NA 
            
            if( x %in% X[which(r[iter,,]==0)] ){ #c
              
              r[iter, j, which(X[j,]==x)] = 0  
              
            }else{ #c
              
              p0 = (unique_tab[j,x])*log(epsilon[j]) +
                prob_new_group_DP(unique_tab[j,x], alpha[1], 
                                  sum( 1- r[iter, ,], na.rm = TRUE), log = TRUE)

              p1 = (unique_tab[j,x])*log((1-epsilon[j]))
              temp = X[j,which(r[iter,j,]==1)]
              n_temp = length(temp)
              if(x %in% temp){ #d
                nn_temp = sum(temp == x)
                p1 = p1  + prob_group_DP(unique_tab[j,x], nn_temp, alpha[j+1], 
                                             n_temp, log = TRUE)
              }else{ #d
                p1 = p1 + prob_new_group_DP(unique_tab[j,x], alpha[j+1], 
                              n_temp, log = TRUE)
              } #d
              
              r[iter, j, which(X[j,]==x)] = sample(c(0,1), 1, 
                                     prob = exp(c(p0, p1)) )
              
            } #c
          } #b
        } #a
        
        for(j in 1:J){
          #sample the mixing proportions
          lj[j] = sum(r[iter, j, ], na.rm = TRUE)
          l0[j] = I[j] - lj[j]
          if(l0[j] == 0){ #a
            z = sample(c(0,1), 1, prob = c(p[1], p[3]/ (I[j]+1)) )
            epsilon[j] = z*rbeta(1, 1,  1 + I[j])
          }else if(l0[j] == I[j]){ #a
            z = sample(c(1,0), 1, prob = c(p[2], p[3]/ (I[j]+1)) )
            epsilon[j] = max(z, rbeta(1, I[j]+1, 1))
          }else{ #a
            epsilon[j] = rbeta(1, l0[j] + 1, lj[j] + 1)
          } #a
        }
      
        #sample alpha 
        K = length(unique(X[which(r[iter,,]==0)]))
        alpha[1] = sample_alpha_DP( alpha[1], K, sum(l0), a, b) #common restaurant
        for(j in 1:J){
          K = length(unique(X[j,which(r[iter,j,]==1)]))
          alpha[j+1] = sample_alpha_DP( alpha[j+1], K, lj[j], a , b )
        }
        
        #compute the prediction probabilities
        l_by_rest = c (sum(l0), lj)
        prob_new_species[iter-1,] = (epsilon * alpha[1]) / (alpha[1] + l_by_rest[1]) + 
          (1 - epsilon) * alpha[-1] / (alpha[-1] + l_by_rest[-1])
      }#mcmc
    
      #estimate prediction prob
      est_prob = colMeans(prob_new_species[(burnin+1):tot_iter,])
      where_vec = which(est_prob == max(est_prob))
      where = ifelse(length(where_vec)>1, sample(where_vec,1), where_vec)
      print(est_prob)
      
      #sample a new observation
      x = sample_from_pop(where, data)

      #check if species is new
      species_discovered[newobs] = !(x %in% X)
      print(c(where,species_discovered[newobs]))
      #add it to the sample
      X[where, I[where]+1] = x
      I[where] = I[where] + 1

      #initializing the restaurants for each run of the MAB
      #r[iter,j,i] is r_{j,i}
      #set to the previous
      temp = r[iter,,]
      if( x %in% X[-where,] ){ #c
        
        temp[which(X==x)] = 0  
        r[1, ,] = temp 
        
      }else{ #c
        
        r[1, where, which(X[where,]==x)] = rep(sample(c(0,1),1, 
                                prob = c(epsilon[where], 1- epsilon[where]))
                                    ,length(which(X[where,]==x)))
        
      } #c
      setTxtProgressBar(pb, newobs)
    }#for MAB
    
  return(cumsum(species_discovered))
}#function



uniform_MAB <- function(data,
                        init_samples = 30, new_samples = 300, 
                        seed = 0){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = list of J pmfs
  ##  init_samples = number of starting observations for the MAB
  ##  new_samples = number of sampling step of the MAB
  ##  seed 
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = new_samples, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  species_discovered = rep(0, new_samples) #vector to save the num of discoveries
  
  set.seed(seed)
  J = length(data) #tot number of populations
  
  X = matrix(NA,nrow = J, 
             ncol = init_samples + new_samples) #matrix of observations X[j,i]is X_{j,i}
  
  #sample initial observations
  for(j in 1:J){
    X[j, 1:init_samples] = sample_from_pop(j, data, size = init_samples) 
  }
  
  I = rep(init_samples, J) #initial sample sizes
  
  for(newobs in 1:new_samples){ #MAB for cycle
   
    #estimate prediction prob
    where = sample(J,1)
    
    #sample a new observation
    x = sample_from_pop(where, data)
    
    #check if species is new
    species_discovered[newobs] = !(x %in% X)

    X[where, I[where]+1] = x
    I[where] = I[where] + 1
    

    setTxtProgressBar(pb, newobs)
  }#for MAB
  
  return(cumsum(species_discovered))
}#function



#+PY multi arm bandit problem 
plusPY_MAB<- function(data, a_alpha = 1, b_alpha = 1,
                      a_sigma = 1, b_sigma = 2, p = rep(1/3, 3),
                      init_samples = 30, new_samples = 300, 
                      burnin = 10, iters = 200, seed = 0, 
                      niter_MH = 1, ada_step = 10,
                      ada_thresh = 0.44,
                      r_ada_input = 0){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = list of J pmfs
  ##  a_alpha, b_alpha = hyperparameters of the gamma prior on the concent params
  ##  a_sigma, b_sigma = hyperparameters of the beta prior on the discount params
  ##  p = hyperparameters of the mixture prior of the mixing proportions
  ##  init_samples = number of starting observations for the MAB
  ##  new_samples = number of sampling step of the MAB
  ##  burnin = length of butnin of each MCMC
  ##  iter = number of iter after burnin of each MCMC
  ##  seed 
  ##  niter_MH = num iter of adaptive metropolis for hyper param
  ##  ada_step, ada_thresh, r_ada_input adaptive metropolis quantitites
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = new_samples, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  species_discovered = rep(0, new_samples) #vector to save the num of discoveries
  
  tot_iter = burnin + iters #per each MCMC
  
  set.seed(seed)
  J = length(data) #tot number of populations
  
  X = matrix(NA,nrow = J, 
             ncol = init_samples + new_samples) #matrix of observations X[j,i] is X_{j,i}
  
  #sample initial observations
  for(j in 1:J){
    X[j, 1:init_samples] = sample_from_pop(j, data, size = init_samples) 
  }
  
  I = rep(init_samples, J) #initial sample sizes
  
  #initializing the hyperparam for the first run of the MAB
  alpha = rgamma(J + 1, shape = a_alpha, rate = b_alpha)
  sigma = rbeta(J + 1, a_sigma, b_sigma)
  epsilon = rep(0.5, J)
  
  # Quantities for adaptive Metropolis for hyper PY
  Prop_sd_logit_sig_j = rep(1, J+1)
  Move_sigma_j_out    = matrix(nrow=J+1, ncol=tot_iter)
  Prop_sd_log_alpha_j = rep(1, J+1)
  Move_alpha_j_out    = matrix(nrow=J+1, ncol=tot_iter)
  
  # Numerically 0 lowerbound hyperpar
  epsilon_MH = 1e-5
  # Numerically infinite upperbound hyperpar
  Max_val = 1e10
  
  #initializing the restaurants for the first run of the MAB
  #r[iter,j,i] is r_{j,i}
  r = array(NA, dim = c(tot_iter+1, J, init_samples + new_samples) )
  
  for(j in 1:J){ #a
    unique_tab = tabulate(X[j,])
    for(x in which(unique_tab>0)){ #b
      if( x %in% X[-j,] ){ #c
        
        r[1, j, which(X[j,]==x)] = 0  
        
      }else{ #c
        
        r[1, j, which(X[j,]==x)] = rep(sample(c(0,1),1, 
                                              prob = c(epsilon[j], 1- epsilon[j]))
                                       ,length(which(X[j,]==x)))
        
      } #c
    } #b
  } #a
  
  
  for(newobs in 1:new_samples){ #MAB for cycle
    
    r_ada = r_ada_input
    
    #MCMC for latent allocation to restaurants
    
    #vector for prob of new species for each mcmc iteration
    prob_new_species = matrix(0, nrow = tot_iter, ncol = J)
    
    #argument of pEPPF (i.e., matrix of counts)
    unique_tab = matrix(0, ncol = max(X, na.rm = TRUE) , nrow = J)
    for(j in 1:J){ #a
      unique_tab[j,] = tabulate(X[j,], nbins = max(X, na.rm = TRUE))
    } #a
    #customer counts in each restaurant
    #l0[j] = num of customers in common restaurant from pop j
    #lj[j] = num of customers in exclusive restraurant from pop j
    l0 = rep(0, J); lj = rep(0, J)
    
    for(iter in 2:(tot_iter+1)){
      #assigned customer to restaurants 
      r[iter,,] = r[iter-1,,]
      for(j in 1:J){ #a
        for(x in which(unique_tab[j,]>0)){ #b
          r[iter, j, which(X[j,]==x)] = NA 
          
          if( x %in% X[which(r[iter,,]==0)] ){ #c
            
            r[iter, j, which(X[j,]==x)] = 0  
            
          }else{ #c
            D = length(unique(X[which(r[iter,,]==0)]))
            p0 = (unique_tab[j,x])*log(epsilon[j]) +
              prob_new_group_PY(unique_tab[j,x], alpha[1], sigma[1], 
                                sum( 1- r[iter, ,], na.rm = TRUE), D, log = TRUE)
            
            p1 = (unique_tab[j,x])*log((1-epsilon[j]))
            temp = X[j,which(r[iter,j,]==1)]
            n_temp = length(temp)
            if(x %in% temp){ #d
              nn_temp = sum(temp == x)
              p1 = p1  + 
                prob_group_PY(unique_tab[j,x], nn_temp, alpha[j+1], sigma[j+1], 
                              n_temp, log = TRUE)
            }else{ #d
              D = length(unique(temp))
              p1 = p1 + 
                prob_new_group_PY(unique_tab[j,x], alpha[j+1], sigma[j+1],
                                  n_temp, D, log = TRUE)
            } #d
            
            r[iter, j, which(X[j,]==x)] = sample(c(0,1), 1, 
                                                 prob = exp(c(p0, p1)) )
            
          } #c
        } #b
      } #a
      
      for(j in 1:J){
        #sample the mixing proportions
        lj[j] = sum(r[iter, j, ], na.rm = TRUE)
        l0[j] = I[j] - lj[j]
        if(l0[j] == 0){ #a
          z = sample(c(0,1), 1, prob = c(p[1], p[3]/ (I[j]+1)) )
          epsilon[j] = z*rbeta(1, 1,  1 + I[j])
        }else if(l0[j] == I[j]){ #a
          z = sample(c(1,0), 1, prob = c(p[2], p[3]/ (I[j]+1)) )
          epsilon[j] = max(z, rbeta(1, I[j]+1, 1))
        }else{ #a
          epsilon[j] = rbeta(1, l0[j] + 1, lj[j] + 1)
        } #a
      }
      
      #sample hyperparam
      n_temp = sum(l0)
      temp = X[which(r[iter,,]==0)]
      nDishes = length(unique(temp))
      Dish_by_rest = nDishes
      
      if(n_temp>0){ #rest empty?
        
        vec_1_to_D_1 = ifelse(nDishes>1, (1:(nDishes-1)), 1)
        ell_d_vec = table(temp)
        
        for (iter_MH in 1:niter_MH){
          # Update parameters \alpha
          sigma_old      = sigma[1]
          alpha_old      = alpha[1]
          log_alpha_old  = log(alpha_old)
          # Propose \alpha
          log_alpha_prop = rnorm(1, mean = log_alpha_old, 
                                 sd=Prop_sd_log_alpha_j[1])
          alpha_prop     = exp(log_alpha_prop)
          
          if(epsilon_MH<alpha_prop && alpha_prop<Max_val){
            # Acc_prob_alpha is on the logarithmic scale (consider Jacobian)
            # Prior and Jacobian part
            Acc_prob_alpha = a_alpha*(log_alpha_prop - log_alpha_old)+
              b_alpha*(alpha_old - alpha_prop)
            
            # Likelihood part
            
            Acc_prob_alpha = Acc_prob_alpha +
              sum(log(alpha_prop + vec_1_to_D_1 * sigma_old) - 
                    log(alpha_old  + vec_1_to_D_1 * sigma_old)) +
              lgamma(alpha_old + n_temp) - lgamma(alpha_prop + n_temp) +
              lgamma(alpha_prop +1) - lgamma(alpha_old +1)
            
            move_alpha         = (log(runif(1)) < Acc_prob_alpha)
            alpha_old          = ifelse(move_alpha, alpha_prop, alpha_old)
            
          } else {
            
            move_alpha = FALSE
            
          }
          
          alpha[1] = alpha_old
          
          # Update parameters \sigma
          log_sigma_old   = log(sigma_old)
          logit_sigma_old = qlogis(sigma_old) # logit function
          # Propose \sigma
          logit_sig_prop = rnorm(1, mean = logit_sigma_old, 
                                 sd=Prop_sd_logit_sig_j[1]) 
          
          log_sigma_prop   = plogis(logit_sig_prop, log=T) # inv logistic function
          sigma_prop       = exp(log_sigma_prop) 
          
          if(epsilon_MH<sigma_prop && sigma_prop<1-epsilon_MH){
            # If we propose something numerically out the parameter space 
            #  we have to reject
            
            
            # Acc_prob_sigma is on the logarithmic scale (consider Jacobian)
            # Prior and Jacobian part
            Acc_prob_sigma = a_sigma*(log_sigma_prop - log_sigma_old)+
              b_sigma*(log(1-sigma_prop) - log(1-sigma_old))
            
            # Likelihood part (it can be made slightly more effiecient TBD)
            Acc_prob_sigma = Acc_prob_sigma + 
              nDishes * (lgamma(1-sigma_old) - 
                           lgamma(1-sigma_prop)) +
              sum(log(alpha[1] + vec_1_to_D_1 * sigma_prop) - 
                    log(alpha[1] + vec_1_to_D_1 * sigma_old)) +
              sum(lgamma(ell_d_vec - sigma_prop) - 
                    lgamma(ell_d_vec - sigma_old))
            
            # End Likelihood part
            
            move_sigma = (log(runif(1)) < Acc_prob_sigma)
            sigma_old  = ifelse(move_sigma, sigma_prop, sigma_old)
            
          } else {
            move_sigma = FALSE
          }
          
          sigma[1] = sigma_old
          
        }#end iter MH
        
        # Save acceptance
        Move_alpha_j_out[1, iter-1] = move_alpha
        Move_sigma_j_out[1, iter-1] = move_sigma
        
      }else{ #rest empty?
        
        alpha[1] = rgamma(1, shape = a_alpha, rate = b_alpha)
        sigma[1] = rbeta(1, a_sigma, b_sigma)
        Move_alpha_j_out[1, iter-1] = TRUE
        Move_sigma_j_out[1, iter-1] = TRUE
        
      }
      
      
      for(j in 1:(J)){ #jr
        temp = X[j,which(r[iter,j,]==1)]
        n_temp = lj[j]
        nDishes = length(unique(temp))
        Dish_by_rest = c(Dish_by_rest, nDishes)
        
        if(lj[j]>0){ #rest empty?
          vec_1_to_D_1 = ifelse(nDishes>1, (1:(nDishes-1)), 1)
          ell_d_vec = table(temp)
          
          for (iter_MH in 1:niter_MH){
            # Update parameters \alpha
            sigma_old      = sigma[j+1]
            alpha_old      = alpha[j+1]
            log_alpha_old  = log(alpha_old)
            # Propose \alpha
            log_alpha_prop = rnorm(1, mean = log_alpha_old, 
                                   sd=Prop_sd_log_alpha_j[j+1])
            alpha_prop     = exp(log_alpha_prop)
            
            if(epsilon_MH<alpha_prop && alpha_prop<Max_val){
              # Acc_prob_alpha is on the logarithmic scale (consider Jacobian)
              # Prior and Jacobian part
              Acc_prob_alpha = a_alpha*(log_alpha_prop - log_alpha_old)+
                b_alpha*(alpha_old-alpha_prop)
              
              # Likelihood part
              
              Acc_prob_alpha = Acc_prob_alpha +
                sum(log(alpha_prop + vec_1_to_D_1 * sigma_old) - 
                      log(alpha_old  + vec_1_to_D_1 * sigma_old)) +
                lgamma(alpha_old + n_temp) - lgamma(alpha_prop + n_temp) +
                lgamma(alpha_prop + 1) - lgamma(alpha_old +1)
              
              move_alpha         = (log(runif(1)) < Acc_prob_alpha)
              alpha_old          = ifelse(move_alpha, alpha_prop, alpha_old)
              
            } else {
              
              move_alpha = FALSE
            }
            
            alpha[j+1] = alpha_old
            
            # Update parameters \sigma
            log_sigma_old   = log(sigma_old)
            logit_sigma_old = qlogis(sigma_old) # logit function
            # Propose \sigma
            logit_sig_prop = rnorm(1, mean = logit_sigma_old, 
                                   sd=Prop_sd_logit_sig_j[j+1]) 
            #
            log_sigma_prop   = plogis(logit_sig_prop, log=T) # inv logistic function
            sigma_prop       = exp(log_sigma_prop) 
            
            if(epsilon_MH<sigma_prop && sigma_prop<1-epsilon_MH){
              # If we propose something numerically out the parameter space 
              #  we have to reject
              
              
              # Acc_prob_sigma is on the logarithmic scale (consider Jacobian)
              # Prior and Jacobian part
              Acc_prob_sigma = a_sigma*(log_sigma_prop - log_sigma_old)+
                b_sigma*(log(1-sigma_prop) - log(1-sigma_old))
              
              # Likelihood part (it can be made slightly more effiecient TBD)
              Acc_prob_sigma = Acc_prob_sigma + 
                nDishes * (lgamma(1-sigma_old) - 
                             lgamma(1-sigma_prop)) +
                sum(log(alpha[j+1] + vec_1_to_D_1 * sigma_prop) - 
                      log(alpha[j+1] + vec_1_to_D_1 * sigma_old)) +
                sum(lgamma(ell_d_vec - sigma_prop) - 
                      lgamma(ell_d_vec - sigma_old))
              
              # End Likelihood part
              
              move_sigma = (log(runif(1)) < Acc_prob_sigma)
              sigma_old  = ifelse(move_sigma, sigma_prop, sigma_old)
              
            } else {
              move_sigma = FALSE
            }
            
            sigma[j+1] = sigma_old
            
          }#end iter MH
          
          Move_alpha_j_out[j+1, iter-1] = move_alpha
          Move_sigma_j_out[j+1, iter-1] = move_sigma
          
        }else{ #rest empty?
          
          alpha[j+1] = rgamma(1, shape = a_alpha, rate = b_alpha)
          sigma[j+1] = rbeta(1, a_sigma, b_sigma)
          Move_alpha_j_out[j+1, iter-1] = TRUE
          Move_sigma_j_out[j+1, iter-1] = TRUE
        }
      } #jr
      
      if((iter-1)%%ada_step == 0){
        r_ada                    = r_ada + ada_step
        ada_delta                = min(0.01, 1/sqrt(iter-1))
        seq_ada_step             = (r_ada-ada_step):(r_ada)
        
        # (Ada)
        # Update proposal for \sigma_j, j = 0, 1, ..., J
        Accept_sigma_j      = apply(Move_sigma_j_out[,seq_ada_step], 1, mean)
        Dec_which_sigma_j   = Accept_sigma_j < ada_thresh
        Prop_sd_logit_sig_j = ifelse(Dec_which_sigma_j, 
                                     exp(log(Prop_sd_logit_sig_j) - ada_delta), 
                                     exp(log(Prop_sd_logit_sig_j) + ada_delta))
        
        # Update proposal for \alpha_j, j = 0, 1, ..., J
        Accept_alpha_j      = apply(Move_alpha_j_out[,seq_ada_step], 1, mean)
        Dec_which_alpha_j   = Accept_alpha_j < ada_thresh
        Prop_sd_log_alpha_j = ifelse(Dec_which_alpha_j, 
                                     exp(log(Prop_sd_log_alpha_j) - ada_delta), 
                                     exp(log(Prop_sd_log_alpha_j) + ada_delta))
      }
      # End Update proposal adaptive MH steps
      
      #compute the prediction probabilities
      l_by_rest = c (sum(l0), lj)
      prob_new_species[iter-1,] = epsilon * 
        (alpha[1] + sigma[1]*Dish_by_rest[1]) / (alpha[1] + l_by_rest[1]) + 
        (1 - epsilon) * 
        (alpha[-1] + sigma[-1]*Dish_by_rest[-1]) / (alpha[-1] + l_by_rest[-1])
    }#mcmc
    
    #estimate prediction prob
    est_prob = colMeans(prob_new_species[(burnin+1):tot_iter,])
    where_vec = which(est_prob == max(est_prob))
    where = ifelse(length(where_vec)>1, sample(where_vec,1), where_vec)
    print(est_prob)
    print(sigma)
    print(epsilon)
    print(l0)
    print(lj)
    #sample a new observation
    x = sample_from_pop(where, data)
    
    #check if species is new
    species_discovered[newobs] = !(x %in% X)
    
    print(c(where,species_discovered[newobs]))
    #add it to the sample
    X[where, I[where]+1] = x
    I[where] = I[where] + 1
    
    #initializing the restaurants for each run of the MAB
    #r[iter,j,i] is r_{j,i}
    #set to the previous
    temp = r[iter,,]
    if( x %in% X[-where,] ){ #c
      
      temp[which(X==x)] = 0  
      r[1, ,] = temp 
      
    }else{ #c
      
      r[1, where, which(X[where,]==x)] = rep(sample(c(0,1),1, 
                                  prob = c(epsilon[where], 1- epsilon[where]))
                                  ,length(which(X[where,]==x)))
      
    } #c
    setTxtProgressBar(pb, newobs)
  }#for MAB
  
  return(cumsum(species_discovered))
}#function

#DP indep multi arm bandit problem 
indepDP_MAB<- function(data, a = 1, b = 1,
                       init_samples = 30, new_samples = 300, 
                       burnin = 10, iters = 200, seed = 0){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = list of J pmfs
  ##  a,b = hyperparameters of the gamma prior on the concentration params
  ##  init_samples = number of starting observations for the MAB
  ##  new_samples = number of sampling step of the MAB
  ##  burnin = length of butnin of each MCMC
  ##  iter = number of iter after burnin of each MCMC
  ##  seed 
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = new_samples, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  species_discovered = rep(0, new_samples) #vector to save the num of discoveries
  
  tot_iter = burnin + iters #per each MCMC
  
  set.seed(seed)
  J = length(data) #tot number of populations
  
  X = matrix(NA,nrow = J, 
             ncol = init_samples + new_samples) #matrix of observations X[j,i]is X_{j,i}
  
  #sample initial observations
  for(j in 1:J){
    X[j, 1:init_samples] = sample_from_pop(j, data, size = init_samples) 
  }
  
  I = rep(init_samples, J) #initial sample sizes
  
  #num of initial unique obs
  K = apply(X, 1, function(x) sum(!is.na(unique(x))))
  
  #initializing the hyperparam for the first run of the MAB
  alpha = rgamma(J, shape = a, rate = b) 
  
  for(newobs in 1:new_samples){ #MAB for cycle
    #MCMC for latent allocation to restaurants
    
    #vector for prob of new species for each mcmc iteration
    prob_new_species = matrix(0, nrow = tot_iter, ncol = J)
    
    for(iter in 2:(tot_iter+1)){
      
      #sample alpha 
      for(j in 1:J){
        alpha[j] = sample_alpha_DP( alpha[j], K[j], I[j], a , b )
      }
      
      #compute the prediction probabilities
      prob_new_species[iter-1,] =  alpha / (alpha + I[j])
    }#mcmc
    
    #estimate prediction prob
    est_prob = colMeans(prob_new_species[(burnin+1):tot_iter,])
    where_vec = which(est_prob == max(est_prob))
    where = ifelse(length(where_vec)>1, sample(where_vec,1), where_vec)
    print(est_prob)
    
    #sample a new observation
    x = sample_from_pop(where, data)
    
    #check if species is new
    species_discovered[newobs] = !(x %in% X)
    if(!(x %in% X)){
      K[where] = K[where]+1
    }
    print(c(where,species_discovered[newobs]))
    #add it to the sample
    X[where, I[where]+1] = x
    I[where] = I[where] + 1
    
    setTxtProgressBar(pb, newobs)
  }#for MAB
  
  return(cumsum(species_discovered))
}#function

#DP indep multi arm bandit problem 
Oracle_MAP<- function(data,
                      init_samples = 30, new_samples = 300, 
                      seed = 0){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = list of J pmfs
  ##  a,b = hyperparameters of the gamma prior on the concentration params
  ##  init_samples = number of starting observations for the MAB
  ##  new_samples = number of sampling step of the MAB
  ##  burnin = length of butnin of each MCMC
  ##  iter = number of iter after burnin of each MCMC
  ##  seed 
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = new_samples, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  species_discovered = rep(0, new_samples) #vector to save the num of discoveries
  
  set.seed(seed)
  J = length(data) #tot number of populations
  
  X = matrix(NA,nrow = J, 
             ncol = init_samples + new_samples) #matrix of observations X[j,i]is X_{j,i}
  
  #sample initial observations
  for(j in 1:J){
    X[j, 1:init_samples] = sample_from_pop(j, data, size = init_samples) 
  }
  
  I = rep(init_samples, J) #initial sample sizes
  
  #normalise data (which is up to norm costant)
  for(j in 1:J){
    data[[j]] = data[[j]] / sum(data[[j]])
  }
  for(newobs in 1:new_samples){ #MAB for cycle
    
    #normalise data (which is up to norm costant)
    
    #estimate prediction prob
    est_prob = NULL
    for(j in 1:J){
      est_prob = c(est_prob, sum(data[[j]][-X[j,!is.na(X[j,])]]))
    }
    print(est_prob)
    where_vec = which(est_prob == max(est_prob))
    where = ifelse(length(where_vec)>1, sample(where_vec,1), where_vec)
    print(est_prob)
    
    #sample a new observation
    x = sample_from_pop(where, data)
    
    #check if species is new
    species_discovered[newobs] = !(x %in% X)
    
    print(c(where,species_discovered[newobs]))
    #add it to the sample
    X[where, I[where]+1] = x
    I[where] = I[where] + 1
    
    setTxtProgressBar(pb, newobs)
  }#for MAB
  
  return(cumsum(species_discovered))
}#function
