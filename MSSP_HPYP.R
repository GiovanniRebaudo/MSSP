# Codes accompanying "Entropy Regularization in Probabilistic Clustering"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.14
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(salso)   
library(ggplot2) 
Fixed_Shared_Hyper = FALSE
Fixed_Diff_Hyper   = FALSE
Random_Diff_Hyper  = FALSE
# Load functions
source("MSSP_fcts.R")
Save_Plot = FALSE

set.seed(123)

# Multivariate species simulations/truth

J           = 3 # Number of populations

I_j_vec     = rep(100,J) 
cum_I_j_vec = cumsum(I_j_vec)
# I_j_vec = (I_1, ...,I_J) vector of sample size in different population
n           = sum(I_j_vec) # tot number of observations

# X_ji_vec = integer(n) 
# Values of all observations (X_{1,1},...,X_{1,I_1},...,X_{J,1},...,X_{J,I_J})
X_ji_vec = c(rep(1,I_j_vec[1]), # X_1i_vec = (X_{1,1},...,X_{1,I_1})
             rep(c(rep(2,15),rep(3,20),rep(1,15)),2), # X_2i_vec
             rep(1:10,10)) # X_3i_vec

Xstar_d_vec = unique(X_ji_vec) # observed dishes (dishes=species)
D           = length(Xstar_d_vec) # overall number of dishes

# Check 
if (n != length(X_ji_vec)){print("error: n != length(X_ji_vec)"); stop()}
if (sum(Xstar_d_vec != 1:D)){
  print("error: labels of dish are not ordered"); stop()}
if (J != length(I_j_vec)){print("error: J != length(I_j_vec)"); stop()}

# Truth of in sample species
table(X_ji_vec) # overall species freq. (n_{.,d})_{d=1}^{D}

# Numbers of different species within pop
K_j_vec = integer(J)
K_j_vec = unique(X_ji_vec[1:I_j_vec[1]]) 
for(j in 1:J){
  lab_ji_vec = 1:I_j_vec[j]
  if(j!=1){lab_ji_vec = lab_ji_vec+cum_I_j_vec[j-1]}
  K_j_vec[j] = length(unique(X_ji_vec[lab_ji_vec]))
}
K_j_vec
cum_K_j_vec = cumsum(K_j_vec)

# Empirical pEPPF unnormalized
emp_pEPPF_un = matrix(0, nrow=D, ncol=J)
for(j in 1:J){
  lab_ji_vec = 1:I_j_vec[j]
  if(j!=1){lab_ji_vec = lab_ji_vec+cum_I_j_vec[j-1]}
  for (d in Xstar_d_vec){
    emp_pEPPF_un[d,j] = sum(X_ji_vec[lab_ji_vec]==d)
  }
}

emp_pEPPF_un
emp_pEPPF_un/n

# Numerically 0
epsilon = 1e-14
# Numerically infinite
M       = 1e100

# MCMC quantities
nGibbsUpdates             = 2e4



set.seed(123)
##### INITIALIZATION
nRest                     = J
nObs                      = n
nDishes                   = D  # number of dishes served in the franchise
dishAllocation            = X_ji_vec

#####
if(FALSE){
  ##### INITIALIZATION TO ALL DIFFERENT TABLES (and some double notation)
  tableAllocation           = 1:nObs
  
  tablesValues              = dishAllocation
  # dish served at each table in the franchise
  
  tableRestaurantAllocation = rep(1:J, times = I_j_vec)
  # allocation of the table to the restaurant its fixed in hssp 
  # (not in nested variations)
  
  nPeopleAtTable            = rep(1,n)
  # people sitting at each table
  
  nTables                   = n 
  # number of occupied tables in the franchise
  
  maxTableIndex             = n 
  # max table index (nTables + nFreeTables = maxTableIndex)
  
  nTablesInRestaurant       = I_j_vec
  # contains only the number of occupied tables in each restaurant
  
  observationDishAllocation = X_ji_vec
  
  # how many people are eating a certain dish
  
  ###
  nFreeTables = 0
  
  freeTables = c() # indices of free tables CONSIDER USING A STACK
  
  tableAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates, ncol = nObs) 
  # tables to which customers are allocated across Gibbs
  
} else if (FALSE){
  ##### INITIALIZATION TO ALL THE SAME TABLE IF SAME DISH AND POPULATION
  
  tableAllocation           = integer(n)
  for(j in 1:J){
    lab_ji_vec = 1:I_j_vec[j]
    past_K_j_vec = 0
    if(j!=1){
      lab_ji_vec   = lab_ji_vec+cum_I_j_vec[j-1]
      past_K_j_vec = cum_K_j_vec[j-1]
    }
    tableAllocation[lab_ji_vec] = X_ji_vec[lab_ji_vec] + past_K_j_vec
  }
  # allocation of customers to tables --> 
  # table indexes are global (across the franchise)
  
  nTables                   = max(tableAllocation)
  # number of occupied tables in the franchise
  
  maxTableIndex             =  n 
  # max table index (nTables + nFreeTables = maxTableIndex)
  
  tablesValues              = dishAllocation
  # dish served at each table in the franchise
  
  
  tableRestaurantAllocation = rep(1:J, times = I_j_vec)
  # allocation of the table to the restaurant it is fixed across mcmc in hssp 
  # (not fixed in nested variations)
  
  nPeopleAtTable            = rep(1,n)
  # people sitting at each table

  nTablesInRestaurant       = K_j_vec
  # contains only the number of occupied tables in each restaurant
  
  observationDishAllocation = X_ji_vec
  # how many people are eating a certain dish
  ###
  nFreeTables = 0
  freeTables = c() # indices of free tables CONSIDER USING A STACK
  tableAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObs) 
  # tables to which customers are allocated across Gibbs
}

# in the following nTables will be the total number of tables, 
# with some of them that might be free, while nTablesInRestaurant 
# will contain only the number of occupied tables in each restaurant


# RANDOM HYPER-PARAMETERS DIFFERENT ACROSS POPULATIONS 
# hyperparameters of Gamma distribution of \theta_j j = 0, 1, ..., J
shape_theta = 1
rate_theta  = 1
# Hyperparameters of Beta distribution of \sigma_j j = 0, 1, ..., J
a_sigma = 1
b_sigma = 1
# Adaptive Metropolis quantities
ada_step   = 50
ada_thresh = 0.44
r_ada      = 0
niter_MH   = 5 # Number of MH iterations used to update hyperparameters
# in each Gibbs iteration

# Quantities for adaptive Metropolis quantities
Prop_sd_logit_sig_j    = rep(1,J+1)
Move_sigma_j_out       = matrix(nrow=J+1, ncol=nGibbsUpdates)
Prop_sd_log_theta_j    = rep(1,J+1)
Move_theta_j_out       = matrix(nrow=J+1, ncol=nGibbsUpdates)

##### INITIALIZATION OF HYPERPARAMETERS WITH THEIR PRIOR MEANS
theta_vec = rep(shape_theta/rate_theta, nRest)
sigma_vec = rep(a_sigma/(a_sigma+b_sigma), nRest)

theta0  = shape_theta/rate_theta
sigma0  = a_sigma/(a_sigma+b_sigma)

### Gibbs Sampler (past tables) 
# RANDOM HYPER-PARAMETERS DIFFERENT ACROSS POPULATIONS 
for (r in 1:nGibbsUpdates) {
  ### ALLOCATE IN-SAMPLE OBSERVATIONS TO TABLES
  if(r%%200==0){print(r)}
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
      
      tableAllocationAcrossGibbs[r,indexCustomerGlobal] = newTableAllocation
      
      indexCustomerGlobal = indexCustomerGlobal + 1
    }
  }
  
  
  # MH within Gibbs step for hyperparameters
  vec_1_to_D   = 1:nDishes
  vec_1_to_D_1 = 1:(nDishes-1)
  
  ell_d_vec = integer(nDishes)
  for (d in unique(dishAllocation)){
    ell_d_vec[d] <- sum(tablesValues == d)
  }

    
  # niter_MH is the number of iteration of the MH within each Gibbs iteration
  # Update parameters \sigma_j \theta_j j = 0, 1, ..., J
  for (iter_MH in niter_MH){
    
    # Update parameters \theta_0
    sigma_old      = sigma0
    theta_old      = theta0
    log_theta_old  = log(theta_old)
    # Propose \theta_0
    log_theta_prop = rnorm(1, mean = log_theta_old, 
                           sd=Prop_sd_log_theta_j[nRest+1])
    # nRest+1 is position of \theta_0
    theta_prop     = exp(log_theta_prop)
    
    # Acc_prob_theta is on the logarithmic scale (consider Jacobian)
    # Prior and Jacobian part
    Acc_prob_theta = shape_theta*(log_theta_prop - log_theta_old)+
      rate_theta*(theta_old-theta_prop)
    
    # Likelihood part
    
    Acc_prob_theta = Acc_prob_theta +
      sum(log(theta_prop + vec_1_to_D *sigma_old) - 
            log(theta_old + vec_1_to_D *sigma_old)) +
      lgamma(theta_old + nTables) - lgamma(theta_old +1) + 
      lgamma(theta_prop +1) - lgamma(theta_old + nTables)
    
    
    move_theta         = (log(runif(1)) < Acc_prob_theta)
    
    theta_old          = ifelse(move_theta, theta_prop, theta_old)
    theta0            = theta_old
    

    
    # Update parameters \sigma_0
    log_sigma_old   = log(sigma_old)
    log_1_sigma_old = log(1-sigma_old)
    logit_sigma_old = log_sigma_old-log_1_sigma_old # logit function
    # Propose \sigma_0
    logit_sig_prop = rnorm(1, mean = logit_sigma_old, 
                           sd=Prop_sd_logit_sig_j[nRest+1]) 
    # nRest+1 is position of \sigma_0
    sigma_prop       = 1/(1+ exp(-logit_sig_prop)) # logistic function
    log_sigma_prop   = log(sigma_prop)
    log_1_sigma_prop = log(1-sigma_prop)
    
    if(epsilon<sigma_prop && sigma_prop<1-epsilon){
      # If we propose something numerically out the parameter space 
      # we have to reject
      log_sigma_prop   = log(sigma_prop)
      log_1_sigma_prop = log(1-sigma_prop)
    
    
    # Acc_prob_theta is on the logarithmic scale (consider Jacobian)
    # Prior and Jacobian part
    Acc_prob_sigma = a_sigma*(log_sigma_prop - log_sigma_old)+
      b_sigma*(log_1_sigma_prop - log_1_sigma_old)
    
    # Likelihood part (it can be made slightly more effiecient TBD)
    Acc_prob_sigma = Acc_prob_sigma + 
      nDishes *(lgamma(1-sigma_old) - lgamma(1-sigma_prop))+
      sum(lgamma(ell_d_vec - sigma_prop) + lgamma(ell_d_vec - sigma_old))+
      sum(log(theta0 + vec_1_to_D_1 * sigma_prop) - 
            log(theta0 + vec_1_to_D_1 * sigma_old))
    # End Likelihood part
    
    move_sigma = (log(runif(1)) < Acc_prob_sigma)
    
    if(move_sigma){
      sigma_old       = sigma_prop
      log_1_sigma_old = log(1-sigma_old)
      logit_sigma_old = log_sigma_old-log_1_sigma_old # logit function
    }
    } else {
      move_sigma = FALSE
    }
    
    sigma0    = sigma_old
    
    # Save acceptance
    if (iter_MH == niter_MH){
      Move_theta_j_out[nRest+1, r] = move_theta
      Move_sigma_j_out[nRest+1, r] = move_sigma
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
      
      # Acc_prob_theta is on the logarithmic scale (consider Jacobian)
      # Prior Gamma and Jacobian part
      Acc_prob_theta = shape_theta*(log_theta_prop - log_theta_old)+
        rate_theta*(theta_old-theta_prop)
      
      # Quantities useful in the log Likelihood part (PYP log EPPF)
      ell_j            = nTablesInRestaurant[indexRestaurant]
      vec_1_to_ell_j   = 1:ell_j
      vec_1_to_ell_j_1 = 1:(ell_j-1)
      I_j              = I_j_vec[j]
      
      # Likelihood part (PYP log EPPF)
      Acc_prob_theta = Acc_prob_theta +
        sum(log(theta_prop + vec_1_to_ell_j *sigma_old) - 
              log(theta_old + vec_1_to_ell_j *sigma_old)) +
        lgamma(theta_old + I_j) - lgamma(theta_old +1) + 
        lgamma(theta_prop +1) - lgamma(theta_old + I_j)
      # End: Likelihood part (PYP log EPPF)
      
      move_theta     = (log(runif(1)) < Acc_prob_theta)
      
      theta_old = ifelse(move_theta, theta_prop, theta_old)
      
      theta_vec[indexRestaurant] = theta_old
      

      # Update parameters \sigma_j, j = 1, ..., J
      log_sigma_old   = log(sigma_old)
      log_1_sigma_old = log(1-sigma_old)
      logit_sigma_old = log_sigma_old-log_1_sigma_old # logit function
      # Propose \sigma_j
      logit_sig_prop = rnorm(1, mean = logit_sigma_old, 
                             sd=Prop_sd_logit_sig_j[indexRestaurant]) 
      sigma_prop       = 1/(1+ exp(-logit_sig_prop)) # logistic function
      
      if(epsilon<sigma_prop && sigma_prop<1-epsilon){
        # If we propose something numerically out the parameter space 
        # we have to reject
        log_sigma_prop   = log(sigma_prop)
        log_1_sigma_prop = log(1-sigma_prop)
        
        
        # Acc_prob_theta is on the logarithmic scale (consider Jacobian)
        # Prior and Jacobian part
        Acc_prob_sigma = a_sigma*(log_sigma_prop - log_sigma_old)+
          b_sigma*(log_1_sigma_prop - log_1_sigma_old)
        
        # Likelihood part (it can be made slightly more effiecient TBD)
        indecesTablesInRestaurant = 
          (1:maxTableIndex)[tableRestaurantAllocation==indexRestaurant]
        q_j_vec = nPeopleAtTable[indecesTablesInRestaurant]
        
        Acc_prob_sigma = Acc_prob_sigma + 
          ell_j *(lgamma(1 - sigma_old) - lgamma(1 - sigma_prop))+
          sum(lgamma(q_j_vec - sigma_prop) + lgamma(q_j_vec - sigma_old))+
          sum(log(theta_old + vec_1_to_ell_j_1 * sigma_prop) - 
                log(theta_old + vec_1_to_ell_j_1 * sigma_old))
        # End Likelihood part
        
        move_sigma                = (log(runif(1)) < Acc_prob_sigma)
        
        if(move_sigma){
          sigma_old       = sigma_prop
          log_1_sigma_old = log(1-sigma_old)
          logit_sigma_old = log_sigma_old-log_1_sigma_old # logit function
        }
      } else {
        move_sigma = FALSE
      }
      
      sigma_vec[indexRestaurant] = sigma_old
    

      # Save acceptance if is the last iteration of MH
      if (iter_MH == niter_MH){
        Move_theta_j_out[indexRestaurant, r] = move_theta
        Move_sigma_j_out[indexRestaurant, r] = move_sigma
      }
      
    }
    
  }
  # End MH within Gibbs step for hyperparameters
  
  # Update proposal adaptive MH steps
  if(r%%ada_step == 0){
    r_ada                    = r_ada + ada_step
    ada_delta                = min(0.01, 1/sqrt(r))
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

prob_new_species_vec = (theta0+nDishes*sigma0)/(nTables + theta0) *
  (theta_vec + sigma_vec * nTablesInRestaurant)/(theta_vec +I_j_vec)



# Optimal arm
which.max(prob_new_species_vec)








