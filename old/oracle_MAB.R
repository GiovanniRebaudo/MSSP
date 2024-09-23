#Oracle_MAB: Multi arm bandit with oracle knowledge of the distributions
#           returns the cumulative number of species discovered
#             and the matrix of estimated prob of new species

oracle_MAB<- function(data, pmfs,
                      init_samples = 30, new_samples = 300){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = observations
  ##  pmfs = list of J pmfs
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
  
  
  species_discovered = rep(0, new_samples) #vector to save the num of discoveries
  prob_new = matrix(NA, nrow = J, ncol = new_samples) #mat to save probs new
  
  J = nrow(data) #tot number of populations
  
  X = matrix(NA,nrow = J, 
             ncol = init_samples + new_samples) #matrix of observations X[j,i]is X_{j,i}
  
  #sample initial observations
  for(j in 1:J){
    X[j, 1:init_samples] = data[j, 1:init_samples] 
  }
  
  I = rep(init_samples, J) #initial sample sizes
  
  #normalize pmfs (which otherwise are up to norm costant)
  for(j in 1:J){
    pmfs[[j]] = pmfs[[j]] / sum(pmfs[[j]])
  }
  for(newobs in 1:new_samples){ #MAB for cycle
  
    #estimate prediction prob
    est_prob = NULL
    for(j in 1:J){
      est_prob = c(est_prob, sum(pmfs[[j]][-X[!is.na(X)]]))
    }
    prob_new[, newobs] = est_prob
    #print(est_prob)
    where_vec = which(est_prob == max(est_prob))
    where = ifelse(length(where_vec)>1, sample(where_vec,1), where_vec)
    
    #sample a new observation
    x = data[where, I[where] + 1] 
    
    #check if species is new
    species_discovered[newobs] = !(x %in% X)
    
    #print(c(where,species_discovered[newobs]))
    #add it to the sample
    X[where, I[where]+1] = x
    I[where] = I[where] + 1
    
    setTxtProgressBar(pb, newobs)
  }#for MAB
  
  return(list(discoveries = cumsum(species_discovered), probs = prob_new))
}#function

