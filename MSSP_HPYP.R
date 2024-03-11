# Codes accompanying "Entropy Regularization in Probabilistic Clustering"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.14
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(salso)   
library(ggplot2) 

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

             
# Hyper-parameter settings
alpha  = 0.5
gamma  = 10 # 10
sigma0 = 0.5
theta0 = 10 # 10
sigma  = 0.5
theta  = 10 # 10

# MCMC quantities
niter  = 2e2

set.seed(123)
####### CRF (check prior genero dal modello)
sigma0 = 0.5
theta0 = 10
sigma  = 0.5
theta  = 10

nGibbsUpdates             = 2e4


##### INITIALIZATION TO ALL DIFFERENT TABLES (and some double notation)
nRest                     = J
nObs                      = n
tableAllocation           = 1:nObs
dishAllocation            = X_ji_vec
tablesValues              = dishAllocation
tableRestaurantAllocation = rep(1:J, times = I_j_vec)
nPeopleAtTable            = rep(1,n)
maxTableIndex             = n # largest table index
nTables                   = n # number of non-empty tables
nTablesInRestaurant       = I_j_vec
observationDishAllocation = X_ji_vec
dishesCounts              = as.vector(table(observationDishAllocation)) 
# how many people are eating a certain dish
nDishes                   = D
##### 


# Initialization
nFreeTables = 0
freeTables = c() # indices of free tables CONSIDER USING A STACK
tableAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObs) 
# tables to which customers are allocated across Gibbs

# in the following nTables will be the total number of tables, 
# with some of them that might be free, while nTablesInRestaurant 
# will contain only the number of occupied tables in each restaurant


### Gibbs Sampler (past tables)
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


### Conditional predictive probabilities





