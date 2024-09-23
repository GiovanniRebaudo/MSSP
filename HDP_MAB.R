HDP_MAB <- function(data,
                    a_alpha = 1, b_alpha = 1/2,
                    a_alpha0 = 1, b_alpha0 = 1/3,
                    init_samples = 30, new_samples = 300, 
                    burnin = 10, iters = 200, seed = 0, 
                    niter_MH = 10, ada_step = 10,
                    ada_thresh = 0.44, r_ada_input = 0){
  ## returns the cumulative number of species discovered
  ## inputs: 
  ##  data = matrix with first batch and future obs
  ##  a_alpha, b_alpha = hyperparameters of the gamma prior on the concent params
  ##  init_samples = number of starting observations for the MAB
  ##  new_samples = number of sampling step of the MAB
  ##  burnin = length of burnin of each MCMC
  ##  iter = number of iter after burnin of each MCMC
  ##  seed 
  ##  niter_MH = num iter of adaptive metropolis for hyper param
  ##  ada_step, ada_thresh, r_ada_input adaptive metropolis quantitites
  
  
  # BEGIN NESTED FUNCTIONS
  
  HDP_MCMC_fct = function(
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
    theta0                    = init_all$theta0,
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
    prob_Table_insample_j = function(model="HDP"){
      if (model=="HDP"){
        
        theta_j = theta_vec[indexRestaurant]
        # Function to compute prob assignment of past tables
        probNewTable = (nTablesServingCurrentDish)/ (nTables + theta0)*theta_j
        
        probs = c(nPeopleAtTable[indecesTablesInRestaurant][indecesPossibleTables],
                  probNewTable)
      } else if (model=="HGnedin"){
        # TBD
      }
      return(probs)
    }
    
    # compute the predictive prob of new species 
    # conditional on tables and hyperparameters values
    prob_new_species_fct <- function(model="HDP"){
      if(model=="HDP"){
        prob_new_species_vec = theta0/(nTables + theta0) *
          theta_vec/(theta_vec + I_j_vec)
      }
      return(prob_new_species_vec)
    }
    
    
    set.seed(seed)
    nRest                       = length(I_j_vec)
    nObs                        = sum(I_j_vec)
    dishAllocation              = Data_vec
    nDishes                     = length(unique(Data_vec))
    
    
    # Quantities where to save output
    prob_new_species = matrix(0,nrow = nGibbsUpdates, ncol = nRest)
    
    if(output=="all"){
      tableAllocationAcrossGibbs = matrix(0, nrow = nGibbsUpdates, ncol = nObs)
      theta_vecAcrossGibbs = matrix(0, nrow = nGibbsUpdates, ncol = nRest)
      nTablesInRestaurantAcrossGibbs = matrix(0, nrow=nGibbsUpdates, ncol=nRest)
      theta0AcrossGibbs = double(nGibbsUpdates)
    }
    
    for (iter in 1:nGibbsUpdates) {
      ### ALLOCATE IN-SAMPLE OBSERVATIONS TO TABLES
      #if(iter%%200==0){print(iter)}
      indexCustomerGlobal = 1
      for (indexRestaurant in 1:nRest) {
        
        for (indexCustomerRestaurant in 1:I_j_vec[indexRestaurant]) {
          indecesTablesInRestaurant = 
            (1:maxTableIndex)[tableRestaurantAllocation==indexRestaurant]
          currentTable = tableAllocation[indexCustomerGlobal] # get the current table
          currentDish  = dishAllocation[indexCustomerGlobal] # get the current dish
          nPeopleAtTable[currentTable] = nPeopleAtTable[currentTable] - 1
          
          if(nPeopleAtTable[currentTable] == 0) { # free the table
            nFreeTables = nFreeTables +1
            freeTables = c(currentTable,freeTables)
            tableRestaurantAllocation[currentTable] = -1
            nTablesInRestaurant[indexRestaurant] = nTablesInRestaurant[indexRestaurant] - 1
            nTables = nTables - 1
            tablesValues[currentTable] = -1
          }
          
          indecesPossibleTables = (tablesValues[indecesTablesInRestaurant] ==
                                     dishAllocation[indexCustomerGlobal])
          
          if(sum(indecesPossibleTables)==0){
            # if no tables in the restaurant is serving the observed dish
            newTableAllocation = -1 # open a new table
          } else {
            # if there are tables in the restaurant serving the observed dish       
            possibleTables = c(indecesTablesInRestaurant[indecesPossibleTables],-1)
            
            nTablesServingCurrentDish = 
              sum(tablesValues == dishAllocation[indexCustomerGlobal])
            
            probs = prob_Table_insample_j(model="HDP")
            
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
              tablesValues[newTableAllocation] = dishAllocation[indexCustomerGlobal]
            } else { # create a new table
              nTablesInRestaurant[indexRestaurant] = 
                nTablesInRestaurant[indexRestaurant] + 1
              maxTableIndex = maxTableIndex + 1
              newTableAllocation = maxTableIndex
              nPeopleAtTable = c(nPeopleAtTable,1)
              tablesValues = c(tablesValues,dishAllocation[indexCustomerGlobal])
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
        theta0 = sample_alpha_DP(theta0, nDishes, nTables, a_alpha0 , b_alpha0)
        for (indexRestaurant in 1:nRest) {
          theta_vec[indexRestaurant] = 
            sample_alpha_DP(theta_vec[indexRestaurant], 
                            nTablesInRestaurant[indexRestaurant], 
                            I_j_vec[indexRestaurant], a_alpha , b_alpha)
        } 
      }
      # compute and save the vector of probabilities of new species
      prob_new_species[iter,] = prob_new_species_fct("HDP")
    }
    
    ## Output
    if(output=="all" && Hyperprior){
      
      return(list(
        tableAllocationAcrossGibbs     = tableAllocationAcrossGibbs,
        theta_vecAcrossGibbs           = theta_vecAcrossGibbs,
        nTablesInRestaurantAcrossGibbs = nTablesInRestaurantAcrossGibbs,
        theta0AcrossGibbs              = theta0AcrossGibbs,
        prob_new_species               = prob_new_species,
        Prop_sd_log_theta_j            = Prop_sd_log_theta_j,
        Move_theta_j_out               = Move_theta_j_out))
      
    } else if (output=="all" && !Hyperprior){
      
      return(list(
        tableAllocationAcrossGibbs     = tableAllocationAcrossGibbs,
        theta_vecAcrossGibbs           = theta_vecAcrossGibbs,
        nTablesInRestaurantAcrossGibbs = nTablesInRestaurantAcrossGibbs,
        theta0AcrossGibbs              = theta0AcrossGibbs,
        prob_new_species               = prob_new_species))
    }
    
    else if (output=="prob_new"){
      
      return(prob_new_species)
      
    } else if (output=="prob and last"){
      
      return(list(
        theta_vec                 = theta_vec,
        theta0                    = theta0,
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
    theta0                    = out$theta0,
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
           theta0                    = theta0,
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
                           model ="HDP", 
                           # model
                           shape_theta, 
                           # common shape of theta_0, ..., theta_J gamma prior
                           rate_theta, 
                           # common rate of theta_0, ..., theta_J gamma prior
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
    
    if(model =="HDP"){
      ##### INITIALIZATION OF HYPERPARAMETERS WITH THEIR PRIOR MEANS
      theta_vec = rep(shape_theta/rate_theta, nRest)
      
      theta0  = a_alpha0/b_alpha0
      
      ## Check parametrization of gamma and beta in R
      # rSamples =rgamma(1000, shape=shape_theta, rate=rate_theta)
      # mean(rSamples); shape_theta/rate_theta
      # var(rSamples)
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
      theta0                    = theta0,
      # theta_0
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
  cat("\n HDP \n")
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = new_samples, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
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
                          model       = "HDP",
                          shape_theta = a_alpha, 
                          rate_theta  = b_alpha)
  
  # Run a short mcmc with fixed hyper par to better initialize the tables
  init_all = HDP_MCMC_fct(
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
    output         = "prob and last"
  )
  
  # Run MCMC
  init_all = HDP_MCMC_fct(
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
    output         = "prob and last"
  )
  
  
  out = HDP_MCMC_fct(
    nGibbsUpdates  = tot_iter,
    seed           = 123,
    # seed to be fixed
    Hyperprior     = T,
    # learn hyperpar via full Bayes if Hyperprior==T
    niter_MH       = niter_MH,
    # number of MH iterations for hyperpar update within each steps
    I_j_vec        = I_j_vec,
    Data_vec       = X_ji_vec,
    ada_step       = ada_step,
    ada_thresh     = ada_thresh,
    r_ada          = r_ada_input,
    shape_theta    = a_alpha, 
    rate_theta     = b_alpha, 
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
      out = HDP_MCMC_fct(
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
        output         = "prob and last"
      )
    }
    setTxtProgressBar(pb, iter_new)
  }#for MAB
  
  return(list(discoveries = cumsum(species_discovered), probs = prob_new_out))
}
