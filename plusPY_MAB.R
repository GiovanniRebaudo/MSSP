#plusPY_MAB:  Multi arm bandit with the +PY, 
#             returns the cumulative number of species discovered
#             and the matrix of estimated prob of new species

plusPY_MAB<- function(data, a_alpha = 1, b_alpha = 1,
                      a_sigma = 1, b_sigma = 2, p = rep(1/3, 3),
                      init_samples = 30, new_samples = 300, 
                      burnin = 100, iters = 400, seed = 0, 
                      niter_MH = 10, ada_step = 5,
                      ada_thresh = 0.44,
                      r_ada_input = 0){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = observations
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
             ncol = init_samples + new_samples) #matrix of observations X[j,i] is X_{j,i}
  
  #sample initial observations
  for(j in 1:J){
    X[j, 1:init_samples] = data[j, 1:init_samples] 
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
        
        if(nDishes>1){
          vec_1_to_D_1 = (1:(nDishes-1))
        }else{
          vec_1_to_D_1 = 1
        }
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
          if(nDishes>1){
            vec_1_to_D_1 = (1:(nDishes-1))
          }else{
            vec_1_to_D_1 = 1
          }
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
  
  return(list(discoveries = cumsum(species_discovered), probs = prob_new))
}#function