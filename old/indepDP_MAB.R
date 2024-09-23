#indepDP_MAB: Multi arm bandit with independent DP, 
#             returns the cumulative number of species discovered
#             and the matrix of estimated prob of new species

indepDP_MAB<- function(data, a = 1, b = 1,
                       init_samples = 30, new_samples = 300, 
                       burnin = 10, iters = 200, seed = 0){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = observations
  ##  a,b = hyperparameters of the gamma prior on the concentration params
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
    prob_new[, newobs] = est_prob
    where_vec = which(est_prob == max(est_prob))
    where = ifelse(length(where_vec)>1, sample(where_vec,1), where_vec)
    #print(est_prob)
    
    #sample a new observation
    x = data[where, I[where] + 1] 
    
    #check if species is new
    species_discovered[newobs] = !(x %in% X)
    if(!(x %in% X)){
      K[where] = K[where]+1
    }
    #print(c(where,species_discovered[newobs]))
    #add it to the sample
    X[where, I[where]+1] = x
    I[where] = I[where] + 1
    
    setTxtProgressBar(pb, newobs)
  }#for MAB
  
  return(list(discoveries = cumsum(species_discovered), probs = prob_new))
}#function
