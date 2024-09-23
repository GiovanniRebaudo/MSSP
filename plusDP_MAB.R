##plusDP_MAB:  Multi arm bandit with the +DP, 
#             returns the cumulative number of species discovered
#             and the matrix of estimated prob of new species

plusDP_MAB<- function(data, a = 6, b = 2,
                      a0 = 1/2, b0 = 2, p = c(0.15,0.15,0.7),
                      init_samples = 30, new_samples = 300, 
                      burnin = 10, iters = 200, seed = 0){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = observations
  ##  a,b = hyperparameters of the gamma prior on the concentration params
  ##  p = hyperparameters of the mixture prior of the mixing proportions
  ##  init_samples = number of starting observations for the MAB
  ##  new_samples = number of sampling step of the MAB
  ##  burnin = length of butnin of each MCMC
  ##  iter = number of iter after burnin of each MCMC
  ##  seed 
  
  if(ncol(data)<(init_samples+new_samples)){
    cat("not enough data provided for", init_samples, "initial samples and", 
        new_samples, "sequential sampling steps")
    cat("\n data should be a matrix with nrow equals to the number of pops")
  }
  
  cat("\n +DP \n")
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = new_samples, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  set.seed(seed)
  J = nrow(data) #tot number of populations
  
  species_discovered = rep(0, new_samples) #vector to save the num of discoveries
  prob_new = matrix(NA, nrow = J, ncol = new_samples) #mat to save probs new
  
  tot_iter = burnin + iters #per each MCMC

  
  X = matrix(NA,nrow = J, 
             ncol = init_samples + new_samples) #matrix of observations X[j,i]is X_{j,i}
  
  #sample initial observations
  for(j in 1:J){
    X[j, 1:init_samples] = data[j, 1:init_samples] 
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
          which_temp = which(X[j,]==x)
          r[iter, j, which_temp] = NA
          
          if( x %in% (X[r[iter,,]==0]) ){ #c
            
            r[iter, j, which_temp] = 0  
            
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
            
            temp = sample(c(0,1), 1,  prob = exp(c(p0, p1)) )
            r[iter, j, which_temp] = temp
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
      alpha[1] = sample_alpha_DP( alpha[1], K, sum(l0), a0, b0) #common restaurant
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
    prob_new[, newobs] = est_prob
    where_vec = which(est_prob == max(est_prob))
    where = ifelse(length(where_vec)>1, sample(where_vec,1), where_vec)
    #print(est_prob)
    
    #sample a new observation
    x = data[where, I[where] + 1] 
    
    #check if species is new
    species_discovered[newobs] = !(x %in% X)
    #print(c(where,species_discovered[newobs]))
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
  
  return(list(discoveries = cumsum(species_discovered), probs = t(prob_new)))
}#function
