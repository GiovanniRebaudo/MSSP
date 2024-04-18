# Functions used in the codes accompanying 
# "Multivariate Species Sampling Process (MSSP)"

# Functions to compute the urn schemes
urn_DP_norm <- function(freq_minus, theta_PY){
  return(c(freq_minus,theta_PY)/(sum(freq_minus)+theta_PY))
}

urn_PY <- function(freq_minus, theta_PY, sigma_PY){
  #  K<-length(freq_minus)
  return(c(freq_minus-sigma_PY,theta_PY+K*sigma_PY))
}

urn_Dir <- function(freq_minus, beta_DM, K_max){
  #  K<-length(freq_minus)
  return(c(freq_minus+beta_DM,beta_DM*(K_max-K)*(K_max>K)))
}

urn_GN_norm <- function(freq_minus, gamma_GN){
  #  K_S<-length(unique(freq_minus))
  tot_freq_minus = sum(freq_minus)
  unorm_out      = c((freq_minus+1)*(tot_freq_minus-K_S+gamma_GN),K_S^2-K_S*gamma_GN)
  return(unorm_out/(tot_freq_minus^2+gamma_GN*tot_freq_minus))
}


# Functions to compute probabilities of possible past (for observed dish) tables 
prob_Table_insample = function(model="HPYP"){
  if (model=="HPYP"){
    # Function to compute prob assignment of past tables
    probNewTable = (nTablesServingCurrentDish - sigma0)/ (nTables + theta0) *
      (theta + sigma * nTablesInRestaurant[indexRestaurant])
    
    probs = c(nPeopleAtTable[indecesTablesInRestaurant][indecesPossibleTables] 
              - sigma, probNewTable)
  } else if (model=="HGnedin"){
    # TBD
  }
  return(probs)
}

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


### Compute the numbers of different species within pop (K_j_vec) and 
### their cumulative (cum_K_j_vec)
K_j_vec_fct <- function(
    I_j_vec,
    # number of observation in each population
    Data_vec
    # vector of data ordered by groups i.e.,
    # X_{1,1}, ..., X_{I_1, 1}, ...., X_{1,J}, ..., X_{I_J, J}
                        ){
  J = length(I_j_vec)
  K_j_vec = integer(J)
  K_j_vec = unique(Data_vec[1:I_j_vec[1]]) 
  for(j in 1:J){
    lab_ji_vec = 1:I_j_vec[j]
    if(j!=1){lab_ji_vec = lab_ji_vec+cum_I_j_vec[j-1]}
    K_j_vec[j] = length(unique(Data_vec[lab_ji_vec]))
  }
  cum_K_j_vec = cumsum(K_j_vec)
  return(list(K_j_vec=K_j_vec, cum_K_j_vec = cumsum(K_j_vec)))
}


### Compute the empirical pEPPF unnormalized
emp_pEPPF_un_fct <- function(
    I_j_vec,
    # number of observation in each population
    Data_vec
    # vector of data ordered by groups i.e.,
    # X_{1,1}, ..., X_{I_1, 1}, ...., X_{1,J}, ..., X_{I_J, J}
){
  J            = length(I_j_vec)
  Xstar_d_vec  = unique(Data_vec)
  D            = length(Xstar_d_vec)
  emp_pEPPF_un = matrix(0, nrow=D, ncol=J)
  
  for(j in 1:J){
    lab_ji_vec = 1:I_j_vec[j]
    if(j!=1){lab_ji_vec = lab_ji_vec+cum_I_j_vec[j-1]}
    for (d in Xstar_d_vec){
      emp_pEPPF_un[d,j] = sum(Data_vec[lab_ji_vec]==d)
    }
  }
  return(emp_pEPPF_un)
}



##### INITIALIZATION HSSP
initHSSP_fct <- function(I_j_vec, 
                     # number of observation in each population
                     Data_vec, 
                     # vector of data ordered by groups i.e.,
                     # X_{1,1}, ..., X_{I_1, 1}, ...., X_{1,J}, ..., X_{I_J, J}
                     model ="HPYP", 
                     # model
                     shape_theta = shape_theta, 
                     # common shape of theta_0, ..., theta_J gamma prior
                     rate_theta = rate_theta, 
                     # common rate of theta_0, ..., theta_J gamma prior
                     a_sigma = a_sigma, 
                     # first hyper of sigma_0, ..., sigma_J beta prior
                     b_sigma = b_sigma, 
                     # second hyper of sigma_0, ..., sigma_J beta prior
                     tablesInit = c("equal", "separate", "manual", "random")
                     # initialization strategy for tables
                     ){
  
  # Check if the dishes are labelled in order of arrival
  uniDish = unique(Data_vec)
  if(!uniDish==sort(uniDish)){
    print("Error: Data are not in order of arrival")
    stop()
  }
  
  nObs                      = n = length(Data_vec)
  nRest                     = J = length(I_j_vec)
  # number of populations
  nDishes                   = D  = length(unique(Data_vec))
  # number of dishes served in the franchise
  dishAllocation            = Data_vec
  
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
    
    tablesValues              = dishAllocation
    # dish served at each table in the franchise
    
    tableRestaurantAllocation = rep(1:nRest, times = I_j_vec)
    # allocation of the table to the restaurant
    
    nPeopleAtTable            = rep(1,n)
    # people sitting at each table
    
    nTables                   = n 
    # number of occupied tables in the franchise
    
    maxTableIndex             = n 
    # max table index (nTables + nFreeTables = maxTableIndex)
    
    nTablesInRestaurant       = I_j_vec
    # contains only the number of occupied tables in each restaurant
    
    observationDishAllocation = Data_vec
    
    # how many people are eating a certain dish
    
    ###
    nFreeTables = 0
    
    freeTables = c() # indices of free tables CONSIDER USING A STACK
    
  } else if (tablesInit == "equal") {
    
    ##### INITIALIZATION TO ALL THE SAME TABLE IF SAME DISH AND POPULATION
    tableAllocation           = integer(n)
    for(j in 1:nRest){
      lab_ji_vec = 1:I_j_vec[j]
      past_K_j_vec = 0
      if(j!=1){
        lab_ji_vec   = lab_ji_vec+cum_I_j_vec[j-1]
        past_K_j_vec = cum_K_j_vec[j-1]
      }
      tableAllocation[lab_ji_vec] = Data_vec[lab_ji_vec] + past_K_j_vec
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
    
    observationDishAllocation = integer(nDishes)
    for(lab_dish in 1:nDishes){
      observationDishAllocation[lab_dish] = sum(Data_vec==lab_dish)
    }
    
    # how many people are eating a certain dish
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
 
  return( # initialization values of
    theta_vec      = theta_vec,
    # theta_1, ..., theta_J
    sigma_vec      = sigma_vec,
    # sigma_1, ..., sigma_J
    theta0         = theta0,
    # theta_0
    sigma0         = sigma0,
    # sigma_0
    nObs           = nObs
    # total number of observations
    nRest          = nRest
    # number of populations
    nDishes        = nDishes              
    # number of dishes served in the franchise
    dishAllocation = dishAllocation          
    
         )
}