#multiarmed bandit for species discovery via mSSP

#functions description: 
#generate_zipf :  returns the pmfs of a number J of zipf distributions 
#                 with overlapping support
#sample_from_pop: returns a sample with replacement from a given pmf
#sample_alpha_DP: returns a sample of the concentration parameter of a DP
#                 from its full conditional in a marginal sampler 
#prob_new_group_DP: ratio between EPPF of DP corresponding to a group 
#                   of equal obs with new value
#prob_group_DP: ratio between EPPF of DP corresponding to a group 
#               of equal obs with old value
#plusDP_MAB:  Multi arm bandit with the +DP, 
#             returns the cumulative number of species discovered
#uniform_MAB: Multi arm bandit with Uniform random choices,
#             returns the cumulative number of species discovered

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
  ## truth = a list of pmfs, truth[[j]][c] = pr(X_{ji} = c)  
  
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
  a = 1
  b = 1
  p = rep(1/3, 3)
  init_samples = 30
  new_samples = 300
  burnin = 10
  iters = 200
  seed = 0
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
        
        r[1, j, which(X[j,]==x)] = 1
        
      } #c
    } #b
  } #a
    
  for(newobs in 1:new_samples){ #MAB for cycle
    #MCMC for latent allocation to restaurants
    
    #vector for prob of new species for each mcmc iteration
    prob_new_species = matrix(0, nrow = tot_iter, ncol = J)
      
    #customer counts in each restaurant
    #l0[j] = num of customers in common restaurant from pop j
    #lj[j] = num of customers in exclusive restraurant from pop j
    l0 = rep(0, J); lj = rep(0, J)
      
    for(iter in 2:(tot_iter+1)){
      #assigned customer to restaurants 
      r[iter,,] = r[iter-1,,]
        for(j in 1:J){ #a
          unique_tab = tabulate(X[j,])
          for(x in which(unique_tab>0)){ #b
            r[iter, j, which(X[j,]==x)] = NA 
            
            if( x %in% X[which(r[iter,,]==0)] ){ #c
              
              r[iter, j, which(X[j,]==x)] = 0  
              
            }else{ #c
              
              p0 = (unique_tab[x])*log(epsilon[j]) +
                prob_new_group_DP(unique_tab[x], alpha[1], 
                                  sum( 1- r[iter, ,], na.rm = TRUE), log = TRUE)

              p1 = (unique_tab[x])*log((1-epsilon[j]))
              temp = X[j,which(r[iter,j,]==1)]
              n_temp = length(temp)
              if(x %in% temp){ #d
                nn_temp = sum(temp == x)
                p1 = p1  + prob_group_DP(unique_tab[x], nn_temp, alpha[j+1], 
                                             n_temp, log = TRUE)
              }else{ #d
                p1 = p1 + prob_new_group_DP(unique_tab[x], alpha[j+1], 
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
      where = sample(which(est_prob == max(est_prob)),1)
      
      #sample a new observation
      x = sample_from_pop(where, data)

      #check if species is new
      species_discovered[newobs] = !(x %in% X)
   
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
        
        r[1, where, which(X[where,]==x)] = 1
        
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
