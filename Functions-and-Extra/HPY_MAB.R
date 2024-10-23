HPY_MAB <- function(data,
                    a_alpha = 1, b_alpha = 1,
                    init_samples = 30, new_samples = 300, 
                    a_sigma = 1, b_sigma = 2,
                    burnin = 10, iters = 200, seed = 0, 
                    niter_MH = 10, ada_step = 5,
                    ada_thresh = 0.44, r_ada_input = 0){
  ## returns the cumulative number of species discovered
  ## inputs: 
  ##  data = matrix with first batch and future obs
  ##  a_alpha, b_alpha = hyperparameters of the gamma prior on the concent params
  ##  a_sigma, b_sigma = hyperparameters of the beta prior on the discount params
  ##  init_samples = number of starting observations for the MAB
  ##  new_samples = number of sampling step of the MAB
  ##  burnin = length of burnin of each MCMC
  ##  iter = number of iter after burnin of each MCMC
  ##  seed 
  ##  niter_MH = num iter of adaptive metropolis for hyper param
  ##  ada_step, ada_thresh, r_ada_input adaptive metropolis quantitites
  
  
  # BEGIN NESTED FUNCTIONS
  
  HPYP_MCMC_fct = function(
    seed       = 123,
    # seed to be fixed
    Hyperprior = F,
    # learn hyperpar via full Bayes if  Hyperprior==T
    niter_MH   = 5,
    # number of MH iterations for hyperpar update within each steps
    I_j_vec = init_all$I_j_vec,
    Data_vec,
    nGibbsUpdates,
    
    # Numerically 0 lowerbound hyperpar
    epsilon = 1e-5,
    # Numerically infinite upperbound hyperpar
    Max_val = 1e10,
    # Initialized values and quantities named for MCMC: Start here
    theta_vec                 = init_all$theta_vec,
    sigma_vec                 = init_all$sigma_vec,
    theta0                    = init_all$theta0,
    sigma0                    = init_all$sigma0,
    tablesValues              = init_all$tablesValues,
    tableAllocation           = init_all$tableAllocation,
    tableRestaurantAllocation = init_all$tableRestaurantAllocation,
    nPeopleAtTable            = init_all$nPeopleAtTable,
    nTables                   = init_all$nTables,
    maxTableIndex             = init_all$maxTableIndex,
    nTablesInRestaurant       = init_all$nTablesInRestaurant,
    observationDishAllocation = init_all$observationDishAllocation,
    nFreeTables               = init_all$nFreeTables,
    freeTables                = init_all$freeTables,
    # Initialized values and quantities named for MCMC: End here
    shape_theta, 
    # common shape of theta_0, ..., theta_J gamma prior
    rate_theta, 
    # common rate of theta_0, ..., theta_J gamma prior
    a_sigma, 
    # first hyper of sigma_0, ..., sigma_J beta prior
    b_sigma, 
    # Adaptive Metropolis quantities
    ada_step    = 50,
    ada_thresh  = 0.44,
    r_ada       = 0,
    # Output prob new, prob new and last iteration 
    #or additional values for checks
    output   = c("prob_new", "prob and last", "all")
  ){
    
    # Functions to compute probabilities of possible past (for observed dish) tables 
    # Different hyperparameters in different populations
    prob_Table_insample_j = function(model="HPYP"){
      if (model=="HPYP"){
        
        theta_j = theta_vec[indexRestaurant]
        sigma_j = sigma_vec[indexRestaurant]
        # Function to compute prob assignment of past tables
        probNewTable = (nTablesServingCurrentDish - sigma0)/ (nTables + theta0) *
          (theta_j + sigma_j * nTablesInRestaurant[indexRestaurant])
        
        probs = c(nPeopleAtTable[indecesTablesInRestaurant][indecesPossibleTables] 
                  - sigma_j, probNewTable)
      } else if (model=="HGnedin"){
        # TBD
      }
      return(probs)
    }
    
    # compute the predictive prob of new species 
    # conditional on tables and hyperparameters values
    prob_new_species_fct <- function(model="HPYP"){
      if(model=="HPYP"){
        prob_new_species_vec = (theta0+nDishes*sigma0)/(nTables + theta0) *
          (theta_vec + sigma_vec * nTablesInRestaurant)/(theta_vec +I_j_vec)
      }
      return(prob_new_species_vec)
    }
    
    
    set.seed(seed)
    nRest                       = length(I_j_vec)
    nObs                        = sum(I_j_vec)
    dishAllocation              = Data_vec
    nDishes                     = length(unique(Data_vec))
    
    if(Hyperprior){
      # Quantities for adaptive Metropolis quantities
      Prop_sd_logit_sig_j = rep(1, nRest+1)
      Move_sigma_j_out    = matrix(nrow=nRest+1, ncol=nGibbsUpdates)
      Prop_sd_log_theta_j = rep(1, nRest+1)
      Move_theta_j_out    = matrix(nrow=nRest+1, ncol=nGibbsUpdates)
      # We can save less if needed e.g., matrix(nrow=J+1, ncol=ada_step)
    }
    
    # Quantities where to save output
    prob_new_species = matrix(0,nrow = nGibbsUpdates, ncol = nRest)
    
    if(output=="all"){
      tableAllocationAcrossGibbs = matrix(0, nrow = nGibbsUpdates, ncol = nObs)
      theta_vecAcrossGibbs = matrix(0, nrow = nGibbsUpdates, ncol = nRest)
      sigma_vecAcrossGibbs = matrix(0, nrow = nGibbsUpdates, ncol = nRest)
      nTablesInRestaurantAcrossGibbs = matrix(0, nrow=nGibbsUpdates, ncol=nRest)
      theta0AcrossGibbs = double(nGibbsUpdates)
      sigma0AcrossGibbs = double(nGibbsUpdates)
    }
    
    for (iter in 1:nGibbsUpdates) {
      ### ALLOCATE IN-SAMPLE OBSERVATIONS TO TABLES
      #if(iter%%200==0){print(iter)}
      indexCustomerGlobal = 1
      for (indexRestaurant in 1:nRest) {
        
        for (indexCustomerRestaurant in 1:I_j_vec[indexRestaurant]) {
          #indecesTablesInRestaurant = 
          #  (1:maxTableIndex)[tableRestaurantAllocation==indexRestaurant] #bea
          currentTable = tableAllocation[indexCustomerGlobal] # get the current table
          currentDish  = dishAllocation[indexCustomerGlobal] # get the current dish
          nPeopleAtTable[currentTable] = nPeopleAtTable[currentTable] - 1
          #print(nPeopleAtTable[currentTable])
          
          if(nPeopleAtTable[currentTable] == 0) { # free the table
            nFreeTables = nFreeTables +1
            freeTables = c(currentTable,freeTables)
            tableRestaurantAllocation[currentTable] = -1
            nTablesInRestaurant[indexRestaurant] = nTablesInRestaurant[indexRestaurant] - 1
            nTables = nTables - 1
            tablesValues[currentTable] = -1
          }
          indecesTablesInRestaurant = 
              (1:maxTableIndex)[tableRestaurantAllocation==indexRestaurant]
          
          indecesPossibleTables = (tablesValues[indecesTablesInRestaurant] ==
                                     currentDish) #dishAllocation[indexCustomerGlobal])
          
          if(sum(indecesPossibleTables)==0){
            # if no tables in the restaurant is serving the observed dish
            newTableAllocation = -1 # open a new table
          } else {
            # if there are tables in the restaurant serving the observed dish       
            possibleTables = c(indecesTablesInRestaurant[indecesPossibleTables],-1)
            
            nTablesServingCurrentDish = 
              sum(tablesValues == currentDish)#dishAllocation[indexCustomerGlobal])
            
            probs = prob_Table_insample_j(model="HPYP")
            
            newTableAllocation = sample(possibleTables, 1, replace = F, prob = probs)
          }
          
          if(newTableAllocation < 0) {
            nTables = nTables + 1
            if(nFreeTables > 0) { # pick the first free table
              newTableAllocation = freeTables[1]
              freeTables = freeTables[-1]
              nFreeTables = nFreeTables - 1
              nPeopleAtTable[newTableAllocation] = 1
              nTablesInRestaurant[indexRestaurant] = 
                nTablesInRestaurant[indexRestaurant] + 1
              tablesValues[newTableAllocation] = currentDish # dishAllocation[indexCustomerGlobal]
            } else { # create a new table
              nTablesInRestaurant[indexRestaurant] = 
                nTablesInRestaurant[indexRestaurant] + 1
              maxTableIndex = maxTableIndex + 1
              newTableAllocation = maxTableIndex
              nPeopleAtTable = c(nPeopleAtTable,1)
              tablesValues = c(tablesValues, currentDish) #dishAllocation[indexCustomerGlobal])
            }
            # assign the table to the restaurant
            tableRestaurantAllocation[newTableAllocation] = indexRestaurant
          } else{ # the sampled table is already occupied in the restaurant --> 
            # just update the relevant quantities
            nPeopleAtTable[newTableAllocation] = 
              nPeopleAtTable[newTableAllocation] + 1
          }
          
          tableAllocation[indexCustomerGlobal] = newTableAllocation
          
          # If we want to save allocation across Gibbs
          # tableAllocationAcrossGibbs[iter,indexCustomerGlobal] = newTableAllocation
          
          indexCustomerGlobal = indexCustomerGlobal + 1
        }
      }
      
      if(Hyperprior){
        # MH within Gibbs step for hyperparameters
        if(nDishes>1){
          vec_1_to_D_1 = 1:(nDishes-1)
        }
        
        
        ell_d_vec = integer(nDishes)
        for (d in unique(dishAllocation)){
          ell_d_vec[d] = sum(tablesValues == d)
        }
        
        # niter_MH is the number of iteration of the MH within each Gibbs iteration
        # Update parameters \sigma_j \theta_j j = 0, 1, ..., J
        for (iter_MH in 1:niter_MH){
          
          # Update parameters \theta_0
          sigma_old      = sigma0
          theta_old      = theta0
          log_theta_old  = log(theta_old)
          # Propose \theta_0
          log_theta_prop = rnorm(1, mean = log_theta_old, 
                                 sd=Prop_sd_log_theta_j[nRest+1])
          # nRest+1 is position of \theta_0
          theta_prop     = exp(log_theta_prop)
          
          if(epsilon<theta_prop && theta_prop<Max_val){
            # Acc_prob_theta is on the logarithmic scale (consider Jacobian)
            # Prior and Jacobian part
            Acc_prob_theta = shape_theta*(log_theta_prop - log_theta_old)+
              rate_theta*(theta_old-theta_prop)
            
            # Likelihood part
            if(nDishes>1){
              Acc_prob_theta = Acc_prob_theta +
                sum(log(theta_prop + vec_1_to_D_1 * sigma_old) - 
                      log(theta_old  + vec_1_to_D_1 * sigma_old)) +
                lgamma(theta_old + nTables) - lgamma(theta_prop + nTables) +
                lgamma(theta_prop +1)       - lgamma(theta_old +1)
            } else {
              Acc_prob_theta = Acc_prob_theta +
                lgamma(theta_old + nTables) - lgamma(theta_prop + nTables) +
                lgamma(theta_prop +1)       - lgamma(theta_old +1)
            }
            
            
            
            move_theta         = (log(runif(1)) < Acc_prob_theta)
            theta_old          = ifelse(move_theta, theta_prop, theta_old)
          } else {
            move_theta = FALSE
          }
          
          theta0            = theta_old
          
          # Update parameters \sigma_0
          log_sigma_old   = log(sigma_old)
          logit_sigma_old = qlogis(sigma_old) # logit function
          # Propose \sigma_0
          logit_sig_prop = rnorm(1, mean = logit_sigma_old, 
                                 sd=Prop_sd_logit_sig_j[nRest+1]) 
          # nRest+1 is position of \sigma_0
          log_sigma_prop   = plogis(logit_sig_prop, log=T) # inv logistic function
          sigma_prop       = exp(log_sigma_prop) 
          
          if(epsilon<sigma_prop && sigma_prop<1-epsilon){
            # If we propose something numerically out the parameter space 
            # we have to reject
            
            
            # Acc_prob_sigma is on the logarithmic scale (consider Jacobian)
            # Prior and Jacobian part
            Acc_prob_sigma = a_sigma*(log_sigma_prop - log_sigma_old)+
              b_sigma*(log(1-sigma_prop) - log(1-sigma_old))
            
            # Likelihood part (it can be made slightly more effiecient TBD)
            if(nDishes>1){
              Acc_prob_sigma = Acc_prob_sigma + 
                nDishes * (lgamma(1-sigma_old) - 
                             lgamma(1-sigma_prop)) +
                sum(log(theta0 + vec_1_to_D_1 * sigma_prop) - 
                      log(theta0 + vec_1_to_D_1 * sigma_old)) +
                sum(lgamma(ell_d_vec - sigma_prop) - 
                      lgamma(ell_d_vec - sigma_old))
            } else {
              Acc_prob_sigma = Acc_prob_sigma + 
                nDishes * (lgamma(1-sigma_old) - 
                             lgamma(1-sigma_prop)) +
                sum(lgamma(ell_d_vec - sigma_prop) - 
                      lgamma(ell_d_vec - sigma_old))
            }
            
            
            # End Likelihood part
            
            move_sigma = (log(runif(1)) < Acc_prob_sigma)
            
            if(move_sigma){
              sigma_old = sigma_prop
              sigma0    = sigma_old
              # logit_sigma_old = qlogis(sigma_old) # logit function
            }
          } else {
            move_sigma = FALSE
          }
          
          # Save acceptance
          if (iter_MH == niter_MH){
            Move_theta_j_out[nRest+1, iter] = move_theta
            Move_sigma_j_out[nRest+1, iter] = move_sigma
          }
          
          
          for (indexRestaurant in 1:nRest) {
            # Update parameters \theta_j, j = 1, ..., J
            sigma_old      = sigma_vec[indexRestaurant]
            theta_old      = theta_vec[indexRestaurant]
            log_theta_old  = log(theta_old)
            # Propose from adaptive Gaussian on transformed space
            log_theta_prop = rnorm(1, mean = log_theta_old, 
                                   sd=Prop_sd_log_theta_j[indexRestaurant])
            theta_prop     = exp(log_theta_prop)
            
            if(epsilon<theta_prop && theta_prop<Max_val){
              
              # Acc_prob_theta is on the logarithmic scale (consider Jacobian)
              # Prior Gamma and Jacobian part
              Acc_prob_theta = shape_theta*(log_theta_prop - log_theta_old)+
                rate_theta*(theta_old-theta_prop)
              
              # Quantities useful in the log Likelihood part (PYP log EPPF)
              I_j              = I_j_vec[indexRestaurant]
              ell_j            = nTablesInRestaurant[indexRestaurant]
              if(ell_j>1){
                vec_1_to_ell_j_1 = 1:(ell_j-1)
                
                # Likelihood part (PYP log EPPF)
                Acc_prob_theta = Acc_prob_theta +
                  sum(log(theta_prop  + vec_1_to_ell_j_1 *sigma_old) - 
                        log(theta_old + vec_1_to_ell_j_1 *sigma_old)) +
                  lgamma(theta_old + I_j) - lgamma(theta_prop + I_j) +
                  lgamma(theta_prop +1)   - lgamma(theta_old +1)
                # End: Likelihood part (PYP log EPPF)
              } else {
                Acc_prob_theta = Acc_prob_theta +
                  lgamma(theta_old + I_j) - lgamma(theta_prop + I_j) +
                  lgamma(theta_prop +1)   - lgamma(theta_old +1)
              }
              
              move_theta     = (log(runif(1)) < Acc_prob_theta)
              
              if (move_theta){
                theta_old = theta_prop
                theta_vec[indexRestaurant] = theta_old
              }
            } else {
              move_theta = FALSE
            }
            
            
            
            # Update parameters \sigma_j, j = 1, ..., J
            log_sigma_old   = log(sigma_old)
            logit_sigma_old = qlogis(sigma_old) # logit function
            # Propose \sigma_j
            logit_sig_prop = rnorm(1, mean = logit_sigma_old, 
                                   sd=Prop_sd_logit_sig_j[indexRestaurant]) 
            log_sigma_prop   = plogis(logit_sig_prop, log=T) # inv logistic function
            sigma_prop       = exp(log_sigma_prop) 
            
            if(epsilon<sigma_prop && sigma_prop<1-epsilon){
              # If we propose something numerically out the parameter space 
              # we have to reject
              
              
              # Acc_prob_sigma is on the logarithmic scale (consider Jacobian)
              # Prior and Jacobian part
              Acc_prob_sigma = a_sigma*(log_sigma_prop - log_sigma_old)+
                b_sigma*(log(1-sigma_prop) - log(1-sigma_old))
              
              # Likelihood part (it can be made slightly more effiecient TBD)
              indecesTablesInRestaurant = 
                (1:maxTableIndex)[tableRestaurantAllocation==indexRestaurant]
              q_j_vec = nPeopleAtTable[indecesTablesInRestaurant]
              
              if(ell_j>1){
                Acc_prob_sigma = Acc_prob_sigma + 
                  ell_j * (lgamma(1 - sigma_old) - 
                             lgamma(1 - sigma_prop)) +
                  sum(lgamma(q_j_vec  - sigma_prop) - 
                        lgamma(q_j_vec  - sigma_old)) +
                  sum(log(theta_old   + vec_1_to_ell_j_1 * sigma_prop) - 
                        log(theta_old + vec_1_to_ell_j_1 * sigma_old))
                # End Likelihood part
              } else {
                Acc_prob_sigma = Acc_prob_sigma + 
                  ell_j * (lgamma(1 - sigma_old) - 
                             lgamma(1 - sigma_prop)) +
                  sum(lgamma(q_j_vec  - sigma_prop) - 
                        lgamma(q_j_vec  - sigma_old))
                # End Likelihood part
              }
              
              move_sigma                = (log(runif(1)) < Acc_prob_sigma)
              
              if(move_sigma){
                sigma_old       = sigma_prop
                # logit_sigma_old = qlogis(sigma_old) # logit function
                sigma_vec[indexRestaurant] = sigma_old
              }
            } else {
              move_sigma = FALSE
            }
            
            
            # Save acceptance if is the last iteration of MH
            if (iter_MH == niter_MH){
              Move_theta_j_out[indexRestaurant, iter] = move_theta
              Move_sigma_j_out[indexRestaurant, iter] = move_sigma
            }
            
          }
        }  
        # End MH within Gibbs step for hyperparameters
        
        # Update proposal adaptive MH steps
        if(iter%%ada_step == 0){
          r_ada                    = r_ada + ada_step
          ada_delta                = min(0.01, 1/sqrt(iter))
          seq_ada_step             = (r_ada-ada_step):r_ada
          
          # (Ada)
          # Update proposal for \sigma_j, j = 0, 1, ..., J
          Accept_sigma_j      = apply(Move_sigma_j_out[,seq_ada_step], 1, mean)
          Dec_which_sigma_j   = Accept_sigma_j < ada_thresh
          Prop_sd_logit_sig_j = ifelse(Dec_which_sigma_j, 
                                       exp(log(Prop_sd_logit_sig_j) - ada_delta), 
                                       exp(log(Prop_sd_logit_sig_j) + ada_delta))
          
          # Update proposal for \theta_j, j = 0, 1, ..., J
          Accept_theta_j      = apply(Move_theta_j_out[,seq_ada_step], 1, mean)
          Dec_which_theta_j   = Accept_theta_j < ada_thresh
          Prop_sd_log_theta_j = ifelse(Dec_which_theta_j, 
                                       exp(log(Prop_sd_log_theta_j) - ada_delta), 
                                       exp(log(Prop_sd_log_theta_j) + ada_delta))
        }
        # End Update proposal adaptive MH steps
      }
      
      # compute and save the vector of probabilities of new species
      prob_new_species[iter,] = prob_new_species_fct("HPYP")
      
      if(output=="all"){
        tableAllocationAcrossGibbs[iter,] = tableAllocation
        theta_vecAcrossGibbs[iter,] = theta_vec
        sigma_vecAcrossGibbs[iter,] = sigma_vec
        nTablesInRestaurantAcrossGibbs[iter,] = nTablesInRestaurant
        theta0AcrossGibbs[iter] = theta0
        sigma0AcrossGibbs[iter] = sigma0
      }
    }
    
    ## Output
    if(output=="all" && Hyperprior){
      
      return(list(
        tableAllocationAcrossGibbs     = tableAllocationAcrossGibbs,
        theta_vecAcrossGibbs           = theta_vecAcrossGibbs,
        sigma_vecAcrossGibbs           = sigma_vecAcrossGibbs,
        nTablesInRestaurantAcrossGibbs = nTablesInRestaurantAcrossGibbs,
        theta0AcrossGibbs              = theta0AcrossGibbs,
        sigma0AcrossGibbs              = sigma0AcrossGibbs,
        prob_new_species               = prob_new_species,
        Prop_sd_logit_sig_j            = Prop_sd_logit_sig_j,
        Move_sigma_j_out               = Move_sigma_j_out,
        Prop_sd_log_theta_j            = Prop_sd_log_theta_j,
        Move_theta_j_out               = Move_theta_j_out))
      
    } else if (output=="all" && !Hyperprior){
      
      return(list(
        tableAllocationAcrossGibbs     = tableAllocationAcrossGibbs,
        theta_vecAcrossGibbs           = theta_vecAcrossGibbs,
        sigma_vecAcrossGibbs           = sigma_vecAcrossGibbs,
        nTablesInRestaurantAcrossGibbs = nTablesInRestaurantAcrossGibbs,
        theta0AcrossGibbs              = theta0AcrossGibbs,
        sigma0AcrossGibbs              = sigma0AcrossGibbs,
        prob_new_species               = prob_new_species))
    } else if (output=="prob_new"){
      
      return(prob_new_species)
      
    } else if (output=="prob and last"){
      
      return(list(
        theta_vec                 = theta_vec,
        sigma_vec                 = sigma_vec,
        theta0                    = theta0,
        sigma0                    = sigma0,
        tablesValues              = tablesValues,
        tableAllocation           = tableAllocation,
        tableRestaurantAllocation = tableRestaurantAllocation,
        nPeopleAtTable            = nPeopleAtTable,
        nTables                   = nTables,
        maxTableIndex             = maxTableIndex,
        nTablesInRestaurant       = nTablesInRestaurant,
        observationDishAllocation = observationDishAllocation,
        nFreeTables               = nFreeTables,
        freeTables                = freeTables,
        dishAllocation            = dishAllocation,
        I_j_vec                   = I_j_vec,
        prob_new_species          = prob_new_species))
    } else {
      print("no output")
    }
  }
  # Initialize (after a new observation is observed) the values of the MCMC
  # iterations given the last iteration of the previous MCMC
  initSeqHSSP_fct <- function(
    newPop                    = NA,
    newDataPoint              = NA,
    theta_vec                 = out$theta_vec,
    sigma_vec                 = out$sigma_vec,
    theta0                    = out$theta0,
    sigma0                    = out$sigma0,
    tablesValues              = out$tablesValues,
    tableAllocation           = out$tableAllocation,
    tableRestaurantAllocation = out$tableRestaurantAllocation,
    nPeopleAtTable            = out$nPeopleAtTable,
    nTables                   = out$nTables,
    maxTableIndex             = out$maxTableIndex,
    nTablesInRestaurant       = out$nTablesInRestaurant,
    observationDishAllocation = out$observationDishAllocation,
    nFreeTables               = out$nFreeTables,
    freeTables                = out$freeTables,
    dishAllocation            = out$dishAllocation,
    I_j_vec                   = out$I_j_vec
  ) {
    
    indexCustomerGlobal = sum(I_j_vec[1:newPop])+1
    
    labels_1toIj     = 1:(indexCustomerGlobal-1)
    if(newPop==1){
      labels_Ij1toI_j = labels_1toIj
    } else {
      labels_Ij1toI_j = (sum(I_j_vec[1:(newPop-1)])+1):(indexCustomerGlobal-1)
    }
    
    # Add new data point
    Data_vec           = c(dishAllocation[labels_1toIj], newDataPoint,
                           dishAllocation[-labels_1toIj])
    # I_j = I_j + 1
    I_j_vec[newPop]    = I_j_vec[newPop]+1
    
    dishAllocation = Data_vec
    
    if (newDataPoint %in% Data_vec[labels_Ij1toI_j]) {
      # If the dish of new obs is already available in the pop we do not add tables
      # tablesValues              = tablesValues
      # Assign obs to the first table in the restaurant serving the dish
      
      currentTable = (1:maxTableIndex)[tablesValues==newDataPoint & 
                                         tableRestaurantAllocation==newPop][1] #bea
      
      tableAllocation           = c(tableAllocation[labels_1toIj], currentTable,
                                    tableAllocation[-labels_1toIj])
      
      tableRestaurantAllocation = tableRestaurantAllocation
      # allocation of the tables to the restaurant
      
      
      nPeopleAtTable[currentTable] = nPeopleAtTable[currentTable] + 1
      nTables                   = nTables
      maxTableIndex             = maxTableIndex
      nTablesInRestaurant       = nTablesInRestaurant
      nFreeTables               = nFreeTables
      freeTables                = freeTables
      
    } else {
      # If the dish of new obs is not available in the pop 
      
      # Assign obs to a new table
      nTables = nTables + 1
      nTablesInRestaurant[newPop] = 
        nTablesInRestaurant[newPop] + 1
      
      if(nFreeTables > 0) { # pick the first free table
        newTableAllocation = freeTables[1]
        freeTables = freeTables[-1]
        nFreeTables = nFreeTables - 1
        nPeopleAtTable[newTableAllocation] = 1
        tablesValues[newTableAllocation] = newDataPoint
      } else { # create a new table
        maxTableIndex = maxTableIndex + 1
        newTableAllocation = maxTableIndex
        nPeopleAtTable = c(nPeopleAtTable,1)
        tablesValues = c(tablesValues, newDataPoint)
      }
      # assign the table to the restaurant
      tableRestaurantAllocation[newTableAllocation] = newPop
      # assign person to the table
      tableAllocation = c(tableAllocation[labels_1toIj], newTableAllocation, 
                          tableAllocation[-labels_1toIj])
      
      nDishes = length(unique(dishAllocation))
      observationDishAllocation = integer(nDishes)
      for(lab_dish in unique(dishAllocation)){
        observationDishAllocation[lab_dish] = sum(Data_vec==lab_dish)
      }
    }
    
    return(
      list(I_j_vec                   = I_j_vec,
           theta_vec                 = theta_vec,
           sigma_vec                 = sigma_vec,
           theta0                    = theta0,
           sigma0                    = sigma0,
           tablesValues              = tablesValues,
           tableAllocation           = tableAllocation,
           tableRestaurantAllocation = tableRestaurantAllocation,
           nPeopleAtTable            = nPeopleAtTable,
           nTables                   = nTables,
           maxTableIndex             = maxTableIndex,
           nTablesInRestaurant       = nTablesInRestaurant,
           observationDishAllocation = observationDishAllocation,
           dishAllocation            = dishAllocation,
           nFreeTables               = nFreeTables,
           freeTables                = freeTables)
    )
  }
  
  
  ##### INITIALIZATION HSSP
  initHSSP_fct <- function(I_j_vec, 
                           # number of observation in each population
                           Data_vec, 
                           # vector of data ordered by groups i.e.,
                           # X_{1,1}, ..., X_{I_1, 1}, ...., X_{1,J}, ..., X_{I_J, J}
                           model ="HPYP", 
                           # model
                           shape_theta, 
                           # common shape of theta_0, ..., theta_J gamma prior
                           rate_theta, 
                           # common rate of theta_0, ..., theta_J gamma prior
                           a_sigma, 
                           # first hyper of sigma_0, ..., sigma_J beta prior
                           b_sigma, 
                           # second hyper of sigma_0, ..., sigma_J beta prior
                           tablesInit = c("equal", "separate", "manual", "random")
                           # initialization strategy for tables
  ){
    
    # Check if the dishes are labelled in order of arrival
    uniDish = unique(Data_vec)
    if(!all.equal(uniDish,sort(uniDish))){
      print("Error: Data are not in order of arrival")
      stop()
    }
    
    nObs                      = length(Data_vec)
    # total number of observations
    nRest                     = length(I_j_vec)
    # number of populations
    nDishes                   = length(uniDish)
    # number of dishes served in the franchise
    dishAllocation            = Data_vec
    
    observationDishAllocation = integer(nDishes)
    for(lab_dish in uniDish){
      observationDishAllocation[lab_dish] = sum(Data_vec==lab_dish)
    }
    # how many people are eating a certain dish
    
    if(model =="HPYP"){
      ##### INITIALIZATION OF HYPERPARAMETERS WITH THEIR PRIOR MEANS
      theta_vec = rep(shape_theta/rate_theta, nRest)
      sigma_vec = rep(a_sigma/(a_sigma+b_sigma), nRest)
      
      theta0  = shape_theta/rate_theta
      sigma0  = a_sigma/(a_sigma+b_sigma)
      
      ## Check parametrization of gamma and beta in R
      # rSamples =rgamma(1000, shape=shape_theta, rate=rate_theta)
      # mean(rSamples); shape_theta/rate_theta
      # var(rSamples)
      # rSamples =rbeta(1000, shape1=a_sigma, shape2=b_sigma)
      # mean(rSamples);  a_sigma/(a_sigma+b_sigma)
      # var(rSamples)
    }
    
    if (tablesInit == "separate"){
      ##### INITIALIZATION TO ALL DIFFERENT TABLES (and some double notation)
      tableAllocation           = 1:nObs
      # allocation of customers to tables --> 
      # table indexes are global (across the franchise)
      
      tablesValues              = dishAllocation
      # dish served at each table in the franchise
      
      tableRestaurantAllocation = rep(1:nRest, times = I_j_vec)
      # allocation of the table to the restaurant
      
      nPeopleAtTable            = rep(1, nObs)
      # people sitting at each table
      
      nTables                   = nObs 
      # number of occupied tables in the franchise
      
      maxTableIndex             = nObs 
      # max table index (nTables + nFreeTables = maxTableIndex)
      
      nTablesInRestaurant       = I_j_vec
      # contains only the number of occupied tables in each restaurant
      
      ###
      nFreeTables = 0
      
      freeTables = c() # indices of free tables CONSIDER USING A STACK
      
    } else if (tablesInit == "equal") {
      ##### INITIALIZATION TO ALL THE SAME TABLE IF SAME DISH AND POPULATION
      
      K_j_vec_    = K_j_vec_fct(I_j_vec=I_j_vec, Data_vec=Data_vec)
      cum_K_j_vec = K_j_vec_$cum_K_j_vec
      K_j_vec     = K_j_vec_$K_j_vec
      cum_I_j_vec = cumsum(I_j_vec)
      
      tableAllocation           = integer(nObs)
      for(j in 1:nRest){
        lab_ji_vec = 1:I_j_vec[j]
        past_K_j_vec = 0
        if(j!=1){
          lab_ji_vec   = lab_ji_vec + cum_I_j_vec[j-1]
          past_K_j_vec = cum_K_j_vec[j-1]
        }
        tableAllocation[lab_ji_vec] = 
          as.integer(factor(Data_vec[lab_ji_vec], 
                            levels = unique(Data_vec[lab_ji_vec]))) + past_K_j_vec
      }
      # allocation of customers to tables --> 
      # table indexes are global (across the franchise)
      
      tableRestaurantAllocation = rep(1:nRest, times = K_j_vec)
      # allocation of the table to the restaurant
      
      nTables                   = max(tableAllocation)
      # number of occupied tables in the franchise
      
      maxTableIndex             =  nTables
      # max table index (nTables + nFreeTables = maxTableIndex)
      
      tablesValues              = integer(nTables)
      for(tab_lab in 1:nTables){
        tablesValues[tab_lab] = unique(Data_vec[tableAllocation == tab_lab])
      }
      # dish served at each table in the franchise
      
      nPeopleAtTable            = integer(nTables)
      for(tab_lab in 1:nTables){
        nPeopleAtTable[tab_lab] = sum(tableAllocation==tab_lab)
      }
      # people sitting at each table
      
      nTablesInRestaurant       = K_j_vec
      # contains only the number of occupied tables in each restaurant
      
      ###
      nFreeTables = 0
      freeTables = c() # indices of free tables CONSIDER USING A STACK
      
      # in the following nTables will be the total number of tables, 
      # with some of them that might be free, while nTablesInRestaurant 
      # will contain only the number of occupied tables in each restaurant
      
    } else if (tablesInit == "manual") {
      
      print("Error: TBD")
      
    } else if (tablesInit == "random") {
      
      print("Error: TBD")
      
    }
    
    return(list(# initialization values of
      theta_vec                 = theta_vec,
      # theta_1, ..., theta_J
      sigma_vec                 = sigma_vec,
      # sigma_1, ..., sigma_J
      theta0                    = theta0,
      # theta_0
      sigma0                    = sigma0,
      # sigma_0
      # nObs                      = nObs,
      # total number of observations
      # nRest                     = nRest,
      # number of populations
      # nDishes                   = nDishes,             
      # number of dishes served in the franchise
      # dishAllocation            = dishAllocation,
      # observed individual vector of species in order of arrival
      tableAllocation           = tableAllocation,
      # allocation of customers to tables --> 
      # table indexes are global (across the franchise)
      tablesValues              = tablesValues,
      # dish served at each table in the franchise
      tableRestaurantAllocation = tableRestaurantAllocation,
      # allocation of the table to the restaurant
      nPeopleAtTable            = nPeopleAtTable,
      # people sitting at each table
      nTables                   = nTables,
      # number of occupied tables in the franchise
      maxTableIndex             = maxTableIndex,
      # max table index (nTables + nFreeTables = maxTableIndex)
      nTablesInRestaurant       = nTablesInRestaurant,
      # contains only the number of occupied tables in each restaurant
      observationDishAllocation = observationDishAllocation,
      # how many people are eating a certain dish
      nFreeTables = nFreeTables,
      # number of free tables
      freeTables = freeTables 
      # indices of free tables CONSIDER USING A STACK
    )    )
  }
  
  ## END NESTED FUNCTIONS
  
  if(ncol(data)<(init_samples+new_samples)){
    cat("not enough data provided for", init_samples, "initial samples and", 
        new_samples, "sequential sampling steps")
    cat("\n data should be a matrix with nrow equals to the number of pops")
  }
  cat("\n HPY \n")
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = new_samples, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 20,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  set.seed(seed)
  J = nrow(data) #tot number of populations
  tot_iter = burnin + iters #per each MCMC
  
  species_discovered = rep(0, new_samples) 
  #vector to save the num of discoveries
  prob_new_out = matrix(0, nrow = new_samples, ncol = J)
  # vector for prob of new species for each mcmc iteration
  
  
  # Divide e reorganize the data
  J           = nrow(data)
  I_j_vec     = rep(init_samples, J)
  n           = sum(I_j_vec)
  
  # Reorder dish in order of arrival by group
  X_ji_mat    = data[,1:init_samples]
  uniqDish    = c()
  for(j in 1:J){
    uniqDish  = unique(c(uniqDish,unique(as.integer(X_ji_mat[j,]))))
  }
  uniqDishall = unique(c(uniqDish,unique(as.integer(data))))
  
  data = plyr::mapvalues(data, 
                         from = uniqDishall,
                         to   = 1:length(uniqDishall),
                         warn_missing = FALSE)
  
  X_ji_vec = c()
  for (j in 1:J){
    X_ji_vec = c(X_ji_vec, data[j,1:init_samples])
  }
  
  # save more quantities in MCMC for debugging and convergence checks
  species_discovered = logical(new_samples)
  
  
  ####### MAB
  #### Initialization Gibbs
  init_all = initHSSP_fct(I_j_vec     = I_j_vec, 
                          Data_vec    = X_ji_vec,
                          tablesInit  = "separate", #"equal"
                          model       = "HPYP",
                          shape_theta = a_alpha, 
                          rate_theta  = b_alpha, 
                          a_sigma     = a_sigma, 
                          b_sigma     = b_sigma)
  
  # Run a short mcmc with fixed hyper par to better initialize the tables
  init_all = HPYP_MCMC_fct(
    nGibbsUpdates  = 1e3,
    seed           = 123,
    # seed to be fixed
    Hyperprior     = F,
    # learn hyperpar via full Bayes if  Hyperprior==T
    niter_MH       = 1,
    # number of MH iterations for hyperpar update within each steps
    I_j_vec        = I_j_vec,
    Data_vec       = X_ji_vec,
    shape_theta    = a_alpha, 
    rate_theta     = b_alpha, 
    a_sigma        = a_sigma, 
    b_sigma        = b_sigma,
    output         = "prob and last"
  )
  
  # Run MCMC
  init_all = HPYP_MCMC_fct(
    nGibbsUpdates  = 1e3,
    seed           = 123,
    # seed to be fixed
    Hyperprior     = T,
    # learn hyperpar via full Bayes if Hyperprior==T
    niter_MH       = niter_MH,
    # number of MH iterations for hyperpar update within each steps
    I_j_vec        = I_j_vec,
    Data_vec       = X_ji_vec,
    shape_theta    = a_alpha, 
    rate_theta     = b_alpha, 
    a_sigma        = a_sigma, 
    b_sigma        = b_sigma,
    output         = "prob and last"
  )
  
  
  out = HPYP_MCMC_fct(
    nGibbsUpdates  = tot_iter,
    seed           = 123,
    # seed to be fixed
    Hyperprior     = T,
    # learn hyperpar via full Bayes if Hyperprior==T
    niter_MH       = niter_MH,
    # number of MH iterations for hyperpar update within each steps
    I_j_vec        = I_j_vec,
    Data_vec       = X_ji_vec,
    shape_theta    = a_alpha, 
    rate_theta     = b_alpha, 
    a_sigma        = a_sigma, 
    b_sigma        = b_sigma,
    output         = "prob and last"
  )
  
  
  for (iter_new in 1:new_samples){
    # save 
    prob_new_species = out$prob_new_species
    dishAllocation   = out$dishAllocation
    
    #estimate prediction prob
    est_prob = colMeans(prob_new_species[(burnin+1):tot_iter,])
    prob_new_out[iter_new, ] = est_prob
    
    # Choose optimal arm
    where_vec = which(est_prob == max(est_prob))
    newj = ifelse(length(where_vec)>1, sample(where_vec,1), where_vec)
    
    # Pick new obs
    newObs = data[newj, out$I_j_vec[newj]+1]
    
    # Check if a new species is discovered
    species_new = !(newObs %in% dishAllocation)
    species_discovered[iter_new] = species_new
    
    if (species_new){
      # Relabel the dishes if new dish
      newObsLab = max(dishAllocation)+1
      
      if(newObs!=newObsLab){
        temp = max(c(dishAllocation,max(data)))+1
        data = plyr::mapvalues(data, 
                               from = c(newObsLab, newObs),
                               to   = c(temp, newObsLab),
                               warn_missing = FALSE)
        newObs = newObsLab
      }
    }
    
    init_all = initSeqHSSP_fct(newPop = newj,
                               newDataPoint = newObs)
    
    if(iter_new<new_samples){
      # Run MCMC
      out = HPYP_MCMC_fct(
        nGibbsUpdates  = tot_iter,
        seed           = 123,
        # seed to be fixed
        Hyperprior     = T,
        # learn hyperpar via full Bayes if Hyperprior==T
        niter_MH       = niter_MH,
        # number of MH iterations for hyperpar update within each steps
        Data_vec       = init_all$dishAllocation,
        ada_step       = ada_step,
        ada_thresh     = ada_thresh,
        r_ada          = r_ada_input,
        shape_theta    = a_alpha, 
        rate_theta     = b_alpha, 
        a_sigma        = a_sigma, 
        b_sigma        = b_sigma,
        output         = "prob and last"
      )
    }
    setTxtProgressBar(pb, iter_new)
  }#for MAB
  
  return(list(discoveries = cumsum(species_discovered), probs = prob_new_out))
}
