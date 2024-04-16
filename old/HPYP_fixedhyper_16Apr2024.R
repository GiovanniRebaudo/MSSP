# HPYP with fixed parameters

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
nDishes                   = D
dishAllocation            = X_ji_vec

#####
if(FALSE){
  ##### INITIALIZATION TO ALL DIFFERENT TABLES (and some double notation)
  tableAllocation           = 1:nObs
  
  tablesValues              = dishAllocation
  # dish served at each table in the franchise
  
  tableRestaurantAllocation = rep(1:J, times = I_j_vec)
  
  nPeopleAtTable            = rep(1,n)
  
  nTables                   = n 
  # number of occupied tables in the franchise
  
  maxTableIndex             = n 
  # max table index (nTables + nFreeTables=maxTableIndex)
  
  nTablesInRestaurant       = I_j_vec
  
  observationDishAllocation = X_ji_vec
  
  dishesCounts              = as.vector(table(observationDishAllocation)) 
  # how many people are eating a certain dish
  ###
  nFreeTables = 0
  
  freeTables = c() # indices of free tables CONSIDER USING A STACK
  
  tableAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObs) 
  # tables to which customers are allocated across Gibbs

if(Fixed_Shared_Hyper){
  # Hyper-parameter settings
  alpha  = 0.5
  gamma  = 10 # 10
  sigma0 = 0.5
  theta0 = 10 # 10
  sigma  = 0.5
  theta  = 10 # 10
  
  
  
  ### Gibbs Sampler (past tables) 
  # FIXED HYPER-PARAMETERS EQUAL ACROSS POPULATIONS
  for (r in 1:nGibbsUpdates) {
    ### ALLOCATE IN-SAMPLE OBSERVATIONS TO TABLES
    if(r%%200==0){print(r)}
    indexCustomerGlobal = 1
    for (indexRestaurant in 1:nRest) {
      
      for (indexCustomerRestaurant in 1:I_j_vec[indexRestaurant]) {
        indecesTablesInRestaurant = 
          (1:maxTableIndex)[tableRestaurantAllocation==indexRestaurant]
        currentTable = tableAllocation[indexCustomerGlobal] # get the current table
        currentDish = dishAllocation[indexCustomerGlobal] # get the current dish
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
          
          probs = prob_Table_insample(model="HPYP")
          
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
            tablesValues = c(tablesValues, dishAllocation[indexCustomerGlobal])
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
  }
  
  ### Conditional predictive probabilities (1-step ahead)
  
  # Compute conditional on past tables and obs the probabilities 
  # of 1-step ahead future observations and tables
  
  # TBD
  
  
  # Function to compute prob of discovering a new species in each population
  # FIXED HYPER-PARAMETERS EQUAL ACROSS POPULATIONS
  prob_new_species_vec = (theta0+nDishes*sigma0)/(nTables + theta0) *
    (theta + sigma * nTablesInRestaurant)/(theta +I_j_vec)
  ##
}

if(Fixed_Shared_Hyper){
  # FIXED HYPER-PARAMETERS DIFFERENT ACROSS POPULATIONS 
  sigma_vec = rep(0.5,J)
  theta_vec = rep(10, J)
  
  ### Gibbs Sampler (past tables) 
  # FIXED HYPER-PARAMETERS DIFFERENT ACROSS POPULATIONS
  for (r in 1:nGibbsUpdates) {
    ### ALLOCATE IN-SAMPLE OBSERVATIONS TO TABLES
    if(r%%200==0){print(r)}
    indexCustomerGlobal = 1
    for (indexRestaurant in 1:nRest) {
      
      for (indexCustomerRestaurant in 1:I_j_vec[indexRestaurant]) {
        indecesTablesInRestaurant = 
          (1:maxTableIndex)[tableRestaurantAllocation==indexRestaurant]
        currentTable = tableAllocation[indexCustomerGlobal] # get the current table
        currentDish = dishAllocation[indexCustomerGlobal] # get the current dish
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
  }
  
  prob_new_species_vec = (theta0+nDishes*sigma0)/(nTables + theta0) *
    (theta_vec + sigma_vec * nTablesInRestaurant)/(theta_vec +I_j_vec)
}