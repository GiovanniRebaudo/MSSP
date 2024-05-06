#indepPY_MAB:  Multi arm bandit with independent PY, 
#             returns the cumulative number of species discovered
indepPY_MAB<- function(data, a_alpha = 1, b_alpha = 1,
                      a_sigma = 1, b_sigma = 2,
                      init_samples = 30, new_samples = 300, 
                      burnin = 10, seed = 0, 
                      niter_MH = 2000, ada_step = 10,
                      ada_thresh = 0.44,
                      r_ada_input = 0){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = observations
  ##  a_alpha, b_alpha = hyperparameters of the gamma prior on the concent params
  ##  a_sigma, b_sigma = hyperparameters of the beta prior on the discount params
  ##  init_samples = number of starting observations for the MAB
  ##  new_samples = number of sampling step of the MAB
  ##  burnin = length of burnin of each MH
  ##  seed 
  ##  niter_MH = num iter of adaptive metropolis for hyper param
  ##  ada_step, ada_thresh, r_ada_input adaptive metropolis quantitites
  
  
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
  
  X = matrix(NA,nrow = J, 
             ncol = init_samples + new_samples) #matrix of observations X[j,i] is X_{j,i}
  
  #sample initial observations
  for(j in 1:J){
    X[j, 1:init_samples] = data[j, 1:init_samples] 
  }
  
  I = rep(init_samples, J) #initial sample sizes
  
  #initializing the hyperparam for the first run of the MAB
  alpha = rgamma(J, shape = a_alpha, rate = b_alpha)
  sigma = rbeta(J, a_sigma, b_sigma)
  
  # Quantities for adaptive Metropolis for hyper PY
  Prop_sd_logit_sig_j = rep(1, J)
  Move_sigma_j_out    = matrix(nrow=J, ncol=niter_MH)
  Prop_sd_log_alpha_j = rep(1, J)
  Move_alpha_j_out    = matrix(nrow=J, ncol=niter_MH)
  
  # Numerically 0 lowerbound hyperpar
  epsilon_MH = 1e-5
  # Numerically infinite upperbound hyperpar
  Max_val = 1e10
  
  #num of initial unique obs
  nDishes = apply(X, 1, function(x) sum(!is.na(unique(x))))
  
  for(newobs in 1:new_samples){ #MAB for cycle
    
    r_ada = r_ada_input
    
    #MCMC for latent allocation to restaurants
    
    #vector for prob of new species for each mcmc iteration
    prob_new_species = matrix(0, nrow = niter_MH, ncol = J)
    
    #sample hyperparam
    
    for(j in 1:(J)){ #jr
      temp = X[j, !is.na(X[j,])]
      n_temp = I[j]
        
      if(nDishes[j]>1){
        vec_1_to_D_1 = (1:(nDishes[j]-1))
      }else{
        vec_1_to_D_1 = 1
      }
      ell_d_vec = table(temp)
          
      for (iter_MH in 1:niter_MH){
        # Update parameters \alpha
        sigma_old      = sigma[j]
        alpha_old      = alpha[j]
        log_alpha_old  = log(alpha_old)
        # Propose \alpha
        log_alpha_prop = rnorm(1, mean = log_alpha_old, 
                                 sd=Prop_sd_log_alpha_j[j])
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
            
        alpha[j] = alpha_old
            
        # Update parameters \sigma
        log_sigma_old   = log(sigma_old)
        logit_sigma_old = qlogis(sigma_old) # logit function
        # Propose \sigma
        logit_sig_prop = rnorm(1, mean = logit_sigma_old, 
                             sd=Prop_sd_logit_sig_j[j]) 
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
             nDishes[j] * (lgamma(1-sigma_old) - 
                             lgamma(1-sigma_prop)) +
                sum(log(alpha[j] + vec_1_to_D_1 * sigma_prop) - 
                      log(alpha[j] + vec_1_to_D_1 * sigma_old)) +
                sum(lgamma(ell_d_vec - sigma_prop) - 
                      lgamma(ell_d_vec - sigma_old))
              
          # End Likelihood part
              
          move_sigma = (log(runif(1)) < Acc_prob_sigma)
          sigma_old  = ifelse(move_sigma, sigma_prop, sigma_old)
              
        } else {
          move_sigma = FALSE
        }
            
        sigma[j] = sigma_old
        
        #compute the prediction probabilities
        if(j == J){
        prob_new_species[iter_MH,] =  
          (alpha + sigma*nDishes) / (alpha + I)}
        Move_alpha_j_out[j, iter_MH] = move_alpha
        Move_sigma_j_out[j, iter_MH] = move_sigma
        
      }#end iter MH
          
    } #jr
      
    if((iter_MH)%%ada_step == 0){
      r_ada                    = r_ada + ada_step
      ada_delta                = min(0.01, 1/sqrt(iter_MH))
      seq_ada_step             = (r_ada-ada_step):(r_ada)
        
      # (Ada)
      # Update proposal for \sigma_j, j = 1, ..., J
      Accept_sigma_j      = apply(Move_sigma_j_out[,seq_ada_step], 1, mean)
      Dec_which_sigma_j   = Accept_sigma_j < ada_thresh
      Prop_sd_logit_sig_j = ifelse(Dec_which_sigma_j, 
                               exp(log(Prop_sd_logit_sig_j) - ada_delta), 
                               exp(log(Prop_sd_logit_sig_j) + ada_delta))
        
      # Update proposal for \alpha_j, j = 1, ..., J
      Accept_alpha_j      = apply(Move_alpha_j_out[,seq_ada_step], 1, mean)
      Dec_which_alpha_j   = Accept_alpha_j < ada_thresh
      Prop_sd_log_alpha_j = ifelse(Dec_which_alpha_j, 
                                 exp(log(Prop_sd_log_alpha_j) - ada_delta), 
                                 exp(log(Prop_sd_log_alpha_j) + ada_delta))
    }
    # End Update proposal adaptive MH steps
      
    #estimate prediction prob
    est_prob = colMeans(prob_new_species[(burnin+1):niter_MH,])
    prob_new[, newobs] = est_prob
    where_vec = which(est_prob == max(est_prob))
    where = ifelse(length(where_vec)>1, sample(where_vec,1), where_vec)
    #print(est_prob)
      
    #sample a new observation
    x = data[where, I[where] + 1] 
      
    #check if species is new
    species_discovered[newobs] = !(x %in% X)
    if(!(x %in% X)){
      nDishes[where] = nDishes[where] + 1
    }
      
    #print(c(where,species_discovered[newobs]))
    #add it to the sample
    X[where, I[where]+1] = x
    I[where] = I[where] + 1
      
    setTxtProgressBar(pb, newobs)
  }#for MAB
  
  return(list(discoveries = cumsum(species_discovered), probs = prob_new))
}#function