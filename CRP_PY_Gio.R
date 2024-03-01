n=1000
sigma = 0.5
theta = 10

CRP_PY = function(n,sigma,theta){
  tableAllocation = integer(length = n)
  nTables = 1
  tableCounts = 1
  
  tableAllocation[1] = 1
  
  for (i in 2:n) {
    tableAllocation[i] = sample(1:(nTables+1),1,replace = F,prob = c(tableCounts-sigma,theta+nTables*sigma))
    
    if(tableAllocation[i] == nTables+1) {
      nTables = nTables + 1
      tableCounts = c(tableCounts,1)
    } else {tableCounts[tableAllocation[i]] = tableCounts[tableAllocation[i]] + 1}
  }
  
  return(list(tableAllocation = tableAllocation, tableCounts = tableCounts, nTables = nTables))
}

nn = seq(from = 10000, to = 50000, by = 1000)
distinctValues = double(length = length(nn))

count = 1
for (n in nn) {
  print(count)
  crp_py = CRP_PY(n,sigma,theta)
  distinctValues[count] = length(table(crp_py$tableAllocation))
  count = count+1
}

plot(sqrt(nn), distinctValues)

####### USEFUL FUNCTIONS
CRF_PY_sample_from_prior = function(n,sigma0,theta0,sigma,theta) {
  # n contains n_1,...,n_J
  nObs = sum(n)
  nRest = length(n)
  tableAllocation = integer(length = nObs) # allocation of customers to tables --> table indexes are global (across the franchise)
  dishAllocation = integer(length = nObs) # allocation of customers to dishes --> dish indexes are global (across the franchise)
  nTables = 0 # number of tables in the whole franchise
  tablesValues = c() # dish served at each table in the frachise
  tableRestaurantAllocation = c() # allocation of the table to the restaurant
  nTablesInRestaurant = integer(length = nRest)
  nPeopleAtTable = c() # people sitting at each table
  nDishes = 0 # number of dishes served in the franchise
  dishesCounts = c() # count of the number of tables serving dish 1,2,... in the franchise
  indexCustomerGlobal = 1
  for (indexRestaurant in 1:nRest) {
    ####### initialize the restaurant
    firstTableInRestaurant = nTables + 1 # first table of the restaurant
    tableAllocation[indexCustomerGlobal] = firstTableInRestaurant
    nTables = nTables + 1
    nPeopleAtTable = c(nPeopleAtTable,1)
    
    if(indexRestaurant == 1){
      # initialize the whole franchise
      tablesValues = 1
      dishAllocation[indexCustomerGlobal] = 1
      dishesCounts = 1
      nDishes = nDishes + 1
    } else{
      # sample the dish served at the first table
      dish = sample(1:(nDishes+1),1,replace = F,prob = c(dishesCounts-sigma0,theta0+nDishes*sigma0))
      if(dish == nDishes+1) {
        nDishes = nDishes + 1
        dishesCounts = c(dishesCounts,1)
      } else {dishesCounts[dish] = dishesCounts[dish] + 1}
      tablesValues = c(tablesValues,dish)
      dishAllocation[indexCustomerGlobal] = dish
    }
    tableRestaurantAllocation[nTables] = indexRestaurant # allocate the table to the restaurant
    indexCustomerGlobal = indexCustomerGlobal + 1
    
    tableCountsInRestaurant = 1 # count of people at the tables in the restaurant
    nTablesInRestaurant[indexRestaurant] = 1 # number of tables in the restaurant
    ####### allocate remaining customers
    for (indexCustomerRestaurant in 2:n[indexRestaurant]) {
      # allocate customers to tables
      tableAllocationInRestaurant = sample(1:(nTablesInRestaurant[indexRestaurant]+1),1,replace = F,prob = c(tableCountsInRestaurant-sigma,theta+nTablesInRestaurant[indexRestaurant]*sigma))
      tableAllocationGlobal = firstTableInRestaurant - 1 + tableAllocationInRestaurant
      tableAllocation[indexCustomerGlobal] = tableAllocationGlobal
      
      if(tableAllocationInRestaurant == nTablesInRestaurant[indexRestaurant]+1) {
        nTablesInRestaurant[indexRestaurant] = nTablesInRestaurant[indexRestaurant] + 1
        nTables = nTables + 1
        tableRestaurantAllocation[nTables] = indexRestaurant # allocate the table to the restaurant
        tableCountsInRestaurant = c(tableCountsInRestaurant,1)
        nPeopleAtTable = c(nPeopleAtTable,1)
        # sample the dish served at the new table
        dish = sample(1:(nDishes+1),1,replace = F,prob = c(dishesCounts-sigma0,theta0+nDishes*sigma0))
        if(dish == nDishes+1) {
          nDishes = nDishes + 1
          dishesCounts = c(dishesCounts,1)
        } else {dishesCounts[dish] = dishesCounts[dish] + 1}
        tablesValues = c(tablesValues,dish)
        dishAllocation[indexCustomerGlobal] = dish
        
      } else {
        tableCountsInRestaurant[tableAllocationInRestaurant] = tableCountsInRestaurant[tableAllocationInRestaurant] + 1
        nPeopleAtTable[tableAllocationGlobal] = nPeopleAtTable[tableAllocationGlobal] + 1
      }
      
      dishAllocation[indexCustomerGlobal] = tablesValues[tableAllocationGlobal]
      
      indexCustomerGlobal = indexCustomerGlobal + 1
    }
  }
  
  return(list(dishAllocation = dishAllocation,
              tableAllocation = tableAllocation,
              tablesValues = tablesValues,
              tableRestaurantAllocation = tableRestaurantAllocation,
              nTables = nTables,
              nTablesInRestaurant = nTablesInRestaurant,
              nPeopleAtTable = nPeopleAtTable,
              dishesCounts = dishesCounts,
              nDishes = nDishes))
}

####### CRF (check prior genero dal modello)
n=c(500, 300, 100)
sigma0 = 0.5
theta0 = 10
sigma = 0.5
theta = 10

# scheme to perform gibbs sampling from CRF with no given observations
set.seed(1234)
crf_py = CRF_PY_sample_from_prior(n,sigma0,theta0,sigma,theta)
crf_py

nGibbsUpdates = 2000
dishAllocation = crf_py$dishAllocation
tableAllocation = crf_py$tableAllocation
tablesValues = crf_py$tablesValues
tableRestaurantAllocation = crf_py$tableRestaurantAllocation
nPeopleAtTable = crf_py$nPeopleAtTable
dishesCounts = crf_py$dishesCounts
nTables = crf_py$nTables
nTablesInRestaurant = crf_py$nTablesInRestaurant
nDishes = crf_py$nDishes
nObs = sum(n)
nRest = length(n)

##### INITIALIZATION TO ALL DIFFERENT TABLES
tableAllocation = 1:nObs
tablesValues = dishAllocation
tableRestaurantAllocation = rep(c(1,2,3),times = n)
nPeopleAtTable = rep(1,nObs)
maxTableIndex = nObs # largest table index
nTables = nObs # number of non-empy tables
nTablesInRestaurant = n
##### 

nFreeTables = 0
nFreeDishes = 0
freeTables = c() # indeces of free tables CONSIDER USING A STACK
freeDishes = c() # indeces of the free dishes CONSIDER USING A STACK

# in the following nTables will be the total number of tables, with some of them that might be free,
# while nTablesInRestaurant willcontain only the number of occupied tables in each restaurant

### SSM TABLE ALLOCATION IN-SAMPLE OBSERVATIONS AND SEQUENTIAL PREDICTIVE SAMPLING FOR OUT-OF-SAMPLE
tableAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObs) # tables to which customers are allocated across gibbs

mOutOfSample = c(100,100,100)
nObsOutOfSample = sum(mOutOfSample)
predictiveTableAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObsOutOfSample) # tables to which out-of-sample obs seat in across gibbs
predictiveDishAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObsOutOfSample) # dish eaten by out-of-sample obs in across gibbs

for (r in 1:nGibbsUpdates) {
  ### ALLOCATE IN-SAMPLE OBSERVATIONS TO TABLES
  if(r%%200==0){print(r)}
  indexCustomerGlobal = 1
  for (indexRestaurant in 1:nRest) {
    
    for (indexCustomerRestaurant in 1:n[indexRestaurant]) {
      indecesTablesInRestaurant = (1:maxTableIndex)[tableRestaurantAllocation==indexRestaurant]
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
      
      indecesPossibleTables = (tablesValues[indecesTablesInRestaurant] == dishAllocation[indexCustomerGlobal])
      possibleTables = c(indecesTablesInRestaurant[indecesPossibleTables],-1)
      
      probNewTable = (theta + sigma*nTablesInRestaurant[indexRestaurant])/(nTables + theta0)
      nTablesServingCurrentDish = sum(tablesValues == dishAllocation[indexCustomerGlobal])
      if(nTablesServingCurrentDish > 0) {
        probNewTable = probNewTable*(nTablesServingCurrentDish - sigma0)
      }
      
      probs = c(nPeopleAtTable[indecesTablesInRestaurant][indecesPossibleTables] - sigma, probNewTable)
      
      newTableAllocation = sample(possibleTables,1,replace = F, prob = probs)
      
      if(newTableAllocation < 0) {
        nTables = nTables + 1
        if(nFreeTables > 0) { # pick the first free table
          newTableAllocation = freeTables[1]
          freeTables = freeTables[-1]
          nFreeTables = nFreeTables - 1
          nPeopleAtTable[newTableAllocation] = 1
          nTablesInRestaurant[indexRestaurant] = nTablesInRestaurant[indexRestaurant] + 1
          tablesValues[newTableAllocation] = dishAllocation[indexCustomerGlobal]
        } else { # create a new table
          nTablesInRestaurant[indexRestaurant] = nTablesInRestaurant[indexRestaurant] + 1
          maxTableIndex = maxTableIndex + 1
          newTableAllocation = maxTableIndex
          nPeopleAtTable = c(nPeopleAtTable,1)
          tablesValues = c(tablesValues,dishAllocation[indexCustomerGlobal])
        }
        # assign the table to the restaurant
        tableRestaurantAllocation[newTableAllocation] = indexRestaurant
      } else{ # the sampled table is already occupied in the restaurant --> just update the relevant quantities
        nPeopleAtTable[newTableAllocation] = nPeopleAtTable[newTableAllocation] + 1
      }
      
      tableAllocation[indexCustomerGlobal] = newTableAllocation
      
      tableAllocationAcrossGibbs[r,indexCustomerGlobal] = newTableAllocation
      
      indexCustomerGlobal = indexCustomerGlobal + 1
    }
  }
  
  ### SAMPLE TABLE AND DISH FOR OUT-OF-SAMPLE OBSERVATIONS
  # first we redifine the quantities we need that will be modified in the predictive step but should not be touched in the in-sample part
  nTablesPredictive = nTables
  maxTableIndexPredictive = maxTableIndex
  nTablesInRestaurantPredictive = nTablesInRestaurant
  nFreeTablesPredictive = nFreeTables
  freeTablesPredictive = freeTables
  nPeopleAtTablePredictive = nPeopleAtTable
  tablesValuesPredictive = tablesValues
  tableRestaurantAllocationPredictive = tableRestaurantAllocation
  nDishesPredictive = nDishes
  dishesCountsPredictive = dishesCounts
  
  indexCustomerPredictive = 1
  for (indexRestaurant in 1:nRest) {
    
    for (indexCustomerRestaurant in 1:mOutOfSample[indexRestaurant]) {
      # sample the table
      indecesTablesInRestaurant = (1:maxTableIndexPredictive)[tableRestaurantAllocationPredictive==indexRestaurant]
      
      probNewTable = theta + nTablesInRestaurantPredictive[indexRestaurant]*sigma
      probs = c(nPeopleAtTablePredictive[indecesTablesInRestaurant] - sigma, probNewTable)
      
      newTableAllocation = sample(c(indecesTablesInRestaurant,-1),1,replace = F, prob = probs)
      
      if(newTableAllocation < 0) { # sample the dish
        nTablesPredictive = nTablesPredictive + 1
        nTablesInRestaurantPredictive[indexRestaurant] = nTablesInRestaurantPredictive[indexRestaurant] + 1
        if(nFreeTablesPredictive > 0) { # pick the first free table
          newTableAllocation = freeTablesPredictive[1]
          freeTablesPredictive = freeTablesPredictive[-1]
          nFreeTablesPredictive = nFreeTablesPredictive - 1
          nPeopleAtTablePredictive[newTableAllocation] = 1
        } else { # create a new table
          maxTableIndexPredictive = maxTableIndexPredictive + 1
          newTableAllocation = maxTableIndexPredictive
          nPeopleAtTablePredictive = c(nPeopleAtTablePredictive,1)
          tablesValuesPredictive = c(tablesValuesPredictive,-1) # the dish of the table will be sampled in the next step
        }
        
        # now sample the dish and update tablesValues (we cannot have free dishes here, as we allocate sequentially, but never reallocate)
        dish = sample(1:(nDishesPredictive+1),1,replace = F,prob = c(dishesCountsPredictive-sigma0,theta0+nDishesPredictive*sigma0))
        if(dish > nDishesPredictive) {
          nDishesPredictive = nDishesPredictive + 1
          dishesCountsPredictive = c(dishesCountsPredictive,1)
        } else {dishesCountsPredictive[dish] = dishesCountsPredictive[dish] + 1}
        
        # update the table value, i.e. the dish served at the table
        tablesValuesPredictive[newTableAllocation] = dish
        
      } else{ # the sampled table is already occupied in the restaurant --> just update the relevant quantities
        nPeopleAtTablePredictive[newTableAllocation] = nPeopleAtTablePredictive[newTableAllocation] + 1
      }
      
      predictiveTableAllocationAcrossGibbs[r,indexCustomerPredictive] = newTableAllocation
      predictiveDishAllocationAcrossGibbs[r,indexCustomerPredictive] = tablesValuesPredictive[newTableAllocation]

      indexCustomerPredictive = indexCustomerPredictive + 1
    }
  }
}

tableAllocationAcrossGibbs


predictiveTableAllocationAcrossGibbs
predictiveDishAllocationAcrossGibbs

burnIn = 400

# distinct values in predictive samples
# nDistValues = matrix(nrow = nGibbsUpdates - burnIn,ncol = nRest)
# cutpointsIndeces = c(0,cumsum(mOutOfSample))
# par(mfrow=c(nRest,1))
# for(i in 1:nRest){
#   leftIndex = cutpointsIndeces[i]+1
#   rightIndex = cutpointsIndeces[i+1]
#   for (r in 1:(nGibbsUpdates-burnIn)) {
#     nDistValues[r,i] = length(table(predictiveDishAllocationAcrossGibbs[r+burnIn,leftIndex:rightIndex]))
#   }
#   hist(nDistValues[,i])
# }

# number of new species in predictive samples
nNewSpecies = matrix(nrow = nGibbsUpdates - burnIn,ncol = nRest)
cutpointsIndeces = c(0,cumsum(mOutOfSample))
par(mfrow=c(nRest,1))
for(i in 1:nRest){
  leftIndex = cutpointsIndeces[i]+1
  rightIndex = cutpointsIndeces[i+1]
  for (r in 1:(nGibbsUpdates-burnIn)) {
    nNewSpecies[r,i] = length(setdiff(predictiveDishAllocationAcrossGibbs[r+burnIn,leftIndex:rightIndex],dishAllocation[leftIndex:rightIndex]))
  }
  hist(nNewSpecies[,i])
}






##############################################
# NHPYP
##############################################

nByGroup = c(300, 200, 200, 100,100) # number of obs in each group
firstIndividuals = 1 + cumsum(c(0,nByGroup[-length(nByGroup)]))
lastIndividuals = cumsum(nByGroup)
alpha = 0.5
gamma = 10
sigma0 = 0.5
theta0 = 10
sigma = 0.5
theta = 10

# simulate observations
set.seed(123)
# groupRestaurantAllocation = c(1,1,2,2,3) # allocation of groups by restaurant
# observationRestaurantAllocation = rep(groupRestaurantAllocation,times = n) # allocation of obs to restaurant
n = c(nByGroup[1]+nByGroup[2],nByGroup[3]+nByGroup[4],nByGroup[5])
crf_py = CRF_PY_sample_from_prior(n,sigma0,theta0,sigma,theta)
crf_py

# define quatities for gibbs sampling
# INITIALIZATION TO ALL DIFFERENT TABLES AND DIFFERENT RESTAURANTS
nGibbsUpdates = 200
nObs = sum(nByGroup)
nGroups = length(nByGroup)
nRest = nGroups
maxRestIndex = nGroups
nTables = nObs # number of occupied tables in the franchise
maxTableIndex = nObs # max table index (nTables+nFreeTables=maxTableIndex)
#nTablesInGroup = nByGroup
# nTablesInRestaurant = nByGroup
observationGroupAllocation = rep(1:nGroups, times = nByGroup)
observationRestaurantAllocation = observationGroupAllocation
observationDishAllocation = crf_py$dishAllocation
observationTableAllocation = 1:nObs
groupRestaurantAllocation = 1:nRest
nGroupsInRestaurant = rep(1,nRest)
tablesValues = observationDishAllocation
tableRestaurantAllocation = observationRestaurantAllocation
nPeopleAtTable = rep(1,nObs)

dishesCounts = as.vector(table(observationDishAllocation)) # how many people are eating a certain dish
nDishes = length(dishesCounts)


nFreeRestaurants = 0
nFreeTables = 0
nFreeDishes = 0
freeRestaurants = c() # indeces of free restaurants CONSIDER USING A STACK
freeTables = c() # indeces of free tables CONSIDER USING A STACK
freeDishes = c() # indeces of the free dishes CONSIDER USING A STACK

# in the following nTables will be the total number of tables, with some of them that might be free,
# while nTablesInRestaurant will contain only the number of occupied tables in each restaurant

observationTableAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObs) # tables to which customers are allocated across gibbs
groupRestaurantAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nGroups) # restaurants to which the groups are allocated

mOutOfSample = c(100,100,100,100,100)
nObsOutOfSample = sum(mOutOfSample)
predictiveTableAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObsOutOfSample) # tables to which out-of-sample obs seat in across gibbs
predictiveObservationDishAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObsOutOfSample) # dish eaten by out-of-sample obs in across gibbs


for (r in 1:nGibbsUpdates) {
  ### ALLOCATE IN-SAMPLE OBSERVATIONS TO TABLES
  if(r%%200==0){print(r)}
  indexCustomerGlobal = 1
  for (indexGroup in 1:nGroups) {
    indexRestaurant = groupRestaurantAllocation[indexGroup]
    
    ##### compute quantities without group J
    # observationTableAllocation_noJ = observationTableAllocation
    # observationTableAllocation_noJ[firstIndividuals[indexGroup]:lastIndividuals[indexGroup]] = -1
    observationRestaurantAllocation_noJ = observationRestaurantAllocation
    observationRestaurantAllocation_noJ[firstIndividuals[indexGroup]:lastIndividuals[indexGroup]] = -1
    # observationDishAllocation_noJ = observationDishAllocation
    # observationDishAllocation_noJ[firstIndividuals[indexGroup]:lastIndividuals[indexGroup]] = -1
    # nPeopleAtTable_noJ = nPeopleAtTable - as.integer(table(factor(observationTableAllocation[observationGroupAllocation==indexGroup], levels = 1:maxTableIndex))) # as.integer(table(factor(observationTableAllocation_noJ, levels = 1:maxTableIndex)))
    # nTables_noJ = sum(nPeopleAtTable_noJ>0)
    # 
    # freeTables_noJ = (1:maxTableIndex)[nPeopleAtTable_noJ == 0]
    # nFreeTables_noJ = length(freeTables_noJ)
    # 
    # tablesValues_noJ = tablesValues
    # tablesValues_noJ[freeTables_noJ] = -1
    # 
    # tableRestaurantAllocation_noJ = tableRestaurantAllocation
    # tableRestaurantAllocation_noJ[freeTables_noJ] = -1
    # 
    # nTablesInRestaurant_noJ = sum(tableRestaurantAllocation_noJ==indexRestaurant)
    nPeopleInRestaurant_noJ = sum(observationRestaurantAllocation_noJ==indexRestaurant)
    # 
    dishesCountsInGroup = as.vector(table(factor(observationDishAllocation[observationGroupAllocation==indexGroup], levels = 1:nDishes))) # number of people eating dish d in group J
    dishesCounts_noJ = dishesCounts - dishesCountsInGroup
    # 
    # nDishes_noJ = nDishes - sum(dishesCounts_noJ==0)
    # 
    # # nTablesServingCurrentDish = sum(tablesValues_noJ == observationDishAllocation[indexCustomerGlobal])
    ##### end computation of quantities without group J

    ##### initialize quantities for computation of acceptance probability for the "current" state
    nPeopleInRestaurantCS = nPeopleInRestaurant_noJ
    dishesCountsCS = dishesCounts_noJ
    ##### end initialization quantities for acceptance prob
    
    
    for (indexCustomerGroup in 1:nByGroup[indexGroup]) {
      logP_CS = 0 # log of P[X_ji=x|____] to be used in MH
      indecesTablesInRestaurant = (1:maxTableIndex)[tableRestaurantAllocation==indexRestaurant]
      currentTable = observationTableAllocation[indexCustomerGlobal] # get the current table
      currentDish = observationDishAllocation[indexCustomerGlobal] # get the current dish
      nPeopleAtTable[currentTable] = nPeopleAtTable[currentTable] - 1
      
      if(nPeopleAtTable[currentTable] == 0) { # free the table
        nFreeTables = nFreeTables +1
        freeTables = c(currentTable,freeTables)
        tableRestaurantAllocation[currentTable] = -1
        nTables = nTables - 1
        tablesValues[currentTable] = -1
      }
      
      nTablesInRestaurant = sum(tableRestaurantAllocation==indexRestaurant)
      
      indecesPossibleTables = (tablesValues[indecesTablesInRestaurant] == observationDishAllocation[indexCustomerGlobal])
      possibleTables = c(indecesTablesInRestaurant[indecesPossibleTables],-1)
      
      probNewTable = (theta + sigma*nTablesInRestaurant)/(nTables + theta0)
      nTablesServingCurrentDish = sum(tablesValues == observationDishAllocation[indexCustomerGlobal])
      if(nTablesServingCurrentDish > 0) {
        probNewTable = probNewTable*(nTablesServingCurrentDish - sigma0)
      }
      
      probs = c(nPeopleAtTable[indecesTablesInRestaurant][indecesPossibleTables] - sigma, probNewTable)
      
      newTableAllocation = sample(possibleTables,1,replace = F, prob = probs)
      
      if(newTableAllocation < 0) {
        nTables = nTables + 1
        if(nFreeTables > 0) { # pick the first free table
          newTableAllocation = freeTables[1]
          freeTables = freeTables[-1]
          nFreeTables = nFreeTables - 1
          nPeopleAtTable[newTableAllocation] = 1
          tablesValues[newTableAllocation] = observationDishAllocation[indexCustomerGlobal]
          tableRestaurantAllocation[newTableAllocation] = indexRestaurant # assign table to restaurant
        } else { # create a new table
          maxTableIndex = maxTableIndex + 1
          newTableAllocation = maxTableIndex
          nPeopleAtTable = c(nPeopleAtTable,1)
          tablesValues = c(tablesValues,observationDishAllocation[indexCustomerGlobal])
          tableRestaurantAllocation = c(tableRestaurantAllocation,indexRestaurant) # assign table to restaurant
        }
        
      } else{ # the sampled table is already occupied in the restaurant --> just update the relevant quantities
        nPeopleAtTable[newTableAllocation] = nPeopleAtTable[newTableAllocation] + 1
      }
      
      observationTableAllocation[indexCustomerGlobal] = newTableAllocation
      
      ##### COMPUTE THE CONTRIBUTION PART OF THE LOG ACCEPTANCE PROBABILITY
      individualsToRemove = c((firstIndividuals[indexGroup]-1+indexCustomerGroup):lastIndividuals[indexGroup])
      nPeopleAtTableCS = nPeopleAtTable - as.integer(table(factor(observationTableAllocation[individualsToRemove], levels = 1:maxTableIndex))) # as.integer(table(factor(observationTableAllocation_noJ, levels = 1:maxTableIndex)))
      nTablesCS = sum(nPeopleAtTableCS>0)
      
      # freeTablesCS = (1:maxTableIndex)[nPeopleAtTableCS == 0]
      # tableRestaurantAllocationCS = tableRestaurantAllocation
      # tableRestaurantAllocationCS[freeTablesCS] = -1
      # nTablesInRestaurantCS = sum(tableRestaurantAllocationCS==indexRestaurant)
      
      nTablesInRestaurantCS = length(table(observationTableAllocation[-individualsToRemove][observationRestaurantAllocation[-individualsToRemove]==indexRestaurant]))
      
      
      #tablesValuesCS = tablesValues
      #tablesValuesCS[freeTablesCS] = -1
      
      # compute log of P[X_ji=x|____]
      if(dishesCountsCS[currentDish] == 0){ # X_ji = "new"
        nDishesCS = nDishesCS + 1
        logP_CS = logP_CS + log(theta + nTablesInRestaurantCS*sigma) - 
                  log(theta + nPeopleInRestaurantCS) + 
                  log(theta0 + nDishesCS*sigma0) -
                  log(theta0 + nTablesCS)
      } else{
        # compute number of tables serving the current dish in the whole franchise, excluding individualsToRemove
        nTablesServingCurrentDishCS = length(table(observationTableAllocation[-individualsToRemove][observationDishAllocation[-individualsToRemove] == currentDish]))
        
        # get the tables in current restaurant serving the current dish (excluding future observations)
        freeTablesCS = (1:maxTableIndex)[nPeopleAtTableCS == 0] # which(nPeopleAtTableCS == 0)
        tableRestaurantAllocationCS = tableRestaurantAllocation
        tableRestaurantAllocationCS[freeTablesCS] = -1
        tablesInRestaurantServingCurrentDish = (1:maxTableIndex)[((tableRestaurantAllocationCS == indexRestaurant)&(tablesValues == currentDish))]
        
        #nTablesInRestaurantServingCurrentDishCS = length(tablesInRestaurantServingCurrentDish)
        
        logP_CS = logP_CS - log(theta + nPeopleInRestaurantCS) +
                  log(sum(nPeopleAtTableCS[tablesInRestaurantServingCurrentDish] - sigma) + 
                        (theta + nTablesInRestaurantCS*sigma)*(nTablesServingCurrentDishCS - sigma0)/(theta0 + nTablesCS))
      }
      
      ### update quantities of "current state"
      nPeopleInRestaurantCS = nPeopleInRestaurantCS + 1
      dishesCountsCS[currentDish] = dishesCountsCS[currentDish] + 1
      
      indexCustomerGlobal = indexCustomerGlobal + 1
    }
    
    #### MH STEP
    # we begin the MH proposal sampling
    nRestMH = nRest
    maxRestIndexMH = maxRestIndex
    groupRestaurantAllocationMH = groupRestaurantAllocation
    nGroupsInRestaurantMH = nGroupsInRestaurant
    nGroupsInRestaurantMH[indexRestaurant] = nGroupsInRestaurantMH[indexRestaurant] - 1
    nFreeRestaurantsMH = nFreeRestaurants
    freeRestaurantsMH = freeRestaurants
    
    if(nGroupsInRestaurantMH[indexRestaurant] == 0) {
      nRestMH = nRestMH - 1
      nFreeRestaurantsMH = nFreeRestaurantsMH + 1
      freeRestaurantsMH = c(indexRestaurant,freeRestaurantsMH)
    }
    
    #### PROPOSE A NEW RESTAURANT AND TABLE ALLOCATION
    #### sample restaurant
    nonEmptyRest = (1:maxRestIndexMH)[nGroupsInRestaurantMH>0] # restaurants with at least one group assigned
    possibleRestaurants = c(nonEmptyRest,-1)
    probs =c(nGroupsInRestaurantMH[nonEmptyRest] - alpha, gamma + alpha*nRestMH)
    
    newRestaurantAllocation = sample(possibleRestaurants,1,replace = F, prob = probs)
    
    if(newRestaurantAllocation < 0) {
      nRestMH = nRestMH + 1
      if(nFreeRestaurantsMH > 0) { # pick the first free restaurant
        newRestaurantAllocation = freeRestaurantsMH[1]
        freeRestaurantsMH = freeRestaurantsMH[-1]
        nFreeRestaurantsMH = nFreeRestaurantsMH - 1
        nGroupsInRestaurantMH[newRestaurantAllocation] = 1
        groupRestaurantAllocationMH[indexGroup] = newRestaurantAllocation # assign group to restaurant
      } else { # create a new restaurant
        maxRestIndexMH = maxRestIndexMH + 1
        newRestaurantAllocation = maxRestIndexMH
        nGroupsInRestaurantMH = c(nGroupsInRestaurantMH,1)
        groupRestaurantAllocationMH = c(groupRestaurantAllocationMH,newRestaurantAllocation) # assign group to restaurant
      }
      
    } else{ # the sampled restaurants contains already some groups --> just update the relevant quantities
      nGroupsInRestaurantMH[newRestaurantAllocation] = nGroupsInRestaurantMH[newRestaurantAllocation] + 1
    }
    
    #### sample table allocation
    
    # initialize MH quantities excluding current group j
    maxTableIndexMH = maxTableIndex
    
    observationTableAllocationMH = observationTableAllocation
    observationTableAllocationMH[firstIndividuals[indexGroup]:lastIndividuals[indexGroup]] = -1
    observationRestaurantAllocationMH = observationRestaurantAllocation_noJ # computed above
    observationDishAllocationMH = observationDishAllocation
    observationDishAllocationMH[firstIndividuals[indexGroup]:lastIndividuals[indexGroup]] = -1
    nPeopleAtTableMH = nPeopleAtTable - as.integer(table(factor(observationTableAllocation[observationGroupAllocation==indexGroup], levels = 1:maxTableIndex))) # as.integer(table(factor(observationTableAllocation_noJ, levels = 1:maxTableIndex)))
    nTablesMH = sum(nPeopleAtTableMH>0)
    
    freeTablesMH = (1:maxTableIndexMH)[nPeopleAtTableMH == 0]
    nFreeTablesMH = length(freeTablesMH)
    
    tablesValuesMH = tablesValues
    tablesValuesMH[freeTablesMH] = -1
    
    tableRestaurantAllocationMH = tableRestaurantAllocation
    tableRestaurantAllocationMH[freeTables_noJ] = -1
    
    nPeopleInRestaurantMH = nPeopleInRestaurant_noJ # computed above
    
    dishesCountsMH = dishesCounts_noJ
    
    nDishesMH = nDishes - sum(dishesCountsMH==0)
    
    # # nTablesServingCurrentDish = sum(tablesValues_noJ == observationDishAllocation[indexCustomerGlobal])
    ##### end computation of quantities without group J
    
    
    # #### UPDATE TABLE CONFIGURATION IN THE CURRENT GROUP (analogous to the sampling from the full conditionals)
    for (indexCustomerGroup in 1:nByGroup[indexGroup]) { # this loop should be joined with the loop for the full conditionals!!!!!
      logP_MH = 0 # log of P[X_ji=x|____ S'] to be used in MH
      indecesTablesInRestaurant = (1:maxTableIndexMH)[tableRestaurantAllocationMH==indexRestaurant]
      currentTable = observationTableAllocation[indexCustomerGlobal] # get the current table
      currentDish = observationDishAllocation[indexCustomerGlobal] # get the current dish
      nPeopleAtTableMH[currentTable] = nPeopleAtTableMH[currentTable] - 1
      
      if(nPeopleAtTableMH[currentTable] == 0) { # free the table
        nFreeTablesMH = nFreeTablesMH +1
        freeTablesMH = c(currentTable,freeTablesMH)
        tableRestaurantAllocationMH[currentTable] = -1
        nTablesMH = nTablesMH - 1
        tablesValuesMH[currentTable] = -1
      }
      
      nTablesInRestaurant = sum(tableRestaurantAllocationMH==indexRestaurant)
      
      indecesPossibleTables = (tablesValuesMH[indecesTablesInRestaurant] == observationDishAllocationMH[indexCustomerGlobal])
      possibleTables = c(indecesTablesInRestaurant[indecesPossibleTables],-1)
      
      probNewTable = (theta + sigma*nTablesInRestaurant)/(nTables + theta0)
      nTablesServingCurrentDish = sum(tablesValues == observationDishAllocation[indexCustomerGlobal])
      if(nTablesServingCurrentDish > 0) {
        probNewTable = probNewTable*(nTablesServingCurrentDish - sigma0)
      }
      
      probs = c(nPeopleAtTable[indecesTablesInRestaurant][indecesPossibleTables] - sigma, probNewTable)
      
      newTableAllocation = sample(possibleTables,1,replace = F, prob = probs)
      
      if(newTableAllocation < 0) {
        nTablesMH = nTablesMH + 1
        if(nFreeTablesMH > 0) { # pick the first free table
          newTableAllocation = freeTablesMH[1]
          freeTablesMH = freeTablesMH[-1]
          nFreeTablesMH = nFreeTablesMH - 1
          nPeopleAtTableMH[newTableAllocation] = 1
          tablesValuesMH[newTableAllocation] = observationDishAllocationMH[indexCustomerGlobal]
          tableRestaurantAllocationMH[newTableAllocation] = indexRestaurant # assign table to restaurant
        } else { # create a new table
          maxTableIndexMH = maxTableIndexMH + 1
          newTableAllocationMH = maxTableIndexMH
          nPeopleAtTableMH = c(nPeopleAtTableMH,1)
          tablesValuesMH = c(tablesValuesMH,observationDishAllocationMH[indexCustomerGlobal])
          tableRestaurantAllocationMH = c(tableRestaurantAllocationMH,indexRestaurant) # assign table to restaurant
        }
        
      } else{ # the sampled table is already occupied in the restaurant --> just update the relevant quantities
        nPeopleAtTableMH[newTableAllocation] = nPeopleAtTableMH[newTableAllocation] + 1
      }
      
      observationTableAllocationMH[indexCustomerGlobal] = newTableAllocation
      
      # compute log of P[X_ji=x|____S']
      if(dishesCountsMH[currentDish] == 0){ # X_ji = "new"
        logP_MH = logP_MH + log(theta + nTablesInRestaurant*sigma) - # nTablesInRestaurant refers to the MH configuration and is defined above
          log(theta + nPeopleInRestaurantMH) + 
          log(theta0 + nDishesMH*sigma0) -
          log(theta0 + nTablesMH)
      } else{
        # compute number of tables serving the current dish in the whole franchise, excluding future individuals
        nTablesServingCurrentDishMH = length(table(observationTableAllocationMH[observationDishAllocationMH == currentDish]))
        
        # get the tables in current restaurant serving the current dish (excluding future observations)
        tablesInRestaurantServingCurrentDish = (1:maxTableIndexMH)[((tableRestaurantAllocationMH == indexRestaurant)&(tablesValuesMH == currentDish))]
        
        logP_MH = logP_MH - log(theta + nPeopleInRestaurantMH) +
          log(sum(nPeopleAtTableMH[tablesInRestaurantServingCurrentDish] - sigma) + 
                (theta + nTablesInRestaurant*sigma)*(nTablesServingCurrentDishMH - sigma0)/(theta0 + nTablesMH))
      }
      
      ### update quantities of MH
      nPeopleInRestaurantMH = nPeopleInRestaurantMH + 1
      dishesCountsMH[currentDish] = dishesCountsMH[currentDish] + 1
      observationDishAllocationMH[indexCustomerGlobal] = currentDish
    }
    pAccept = min(exp(logP_MH - logP_CS),1)
    accept = runif(1) < pAccept
    print(pAccept)
    if(accept == T){
      nRest = nRestMH
      maxRestIndex = maxRestIndexMH
      groupRestaurantAllocation[indexGroup] = groupRestaurantAllocationMH[indexGroup]
      nGroupsInRestaurant = nGroupsInRestaurantMH
      nFreeRestaurants = nFreeRestaurantsMH
      freeRestaurants = freeRestaurantsMH
      
      maxTableIndex = maxTableIndexMH
      observationTableAllocation = observationTableAllocationMH
      observationRestaurantAllocation = observationRestaurantAllocationMH
      nPeopleAtTable = nPeopleAtTableMH
      nTables = nTablesMH
      freeTables = freeTablesMH
      nFreeTables = nFreeTablesMH
      tablesValues = tablesValuesMH
      tableRestaurantAllocation = tableRestaurantAllocationMH
    }
    
    ##### END OF METROPOLIS HASTINGS STEP
    
    
    ### EDIT THIS PART TO UPDATE THE TABLES ALLOCATIONS FOR THE WHOLE GROUP
    groupRestaurantAllocationAcrossGibbs[r, indexGroup] = groupRestaurantAllocation[indexGroup]
    observationTableAllocationAcrossGibbs[r,firstIndividuals[indexGroup]:lastIndividuals[indexGroup]] = observationTableAllocation[firstIndividuals[indexGroup]:lastIndividuals[indexGroup]]
  }
  
  ### SAMPLE TABLE AND DISH FOR OUT-OF-SAMPLE OBSERVATIONS
  # THIS HAS BEEN MOVED OUT OF THE LOOP FOR THE MOMENT
}









####################### PREDICTIVE PART
# first we redifine the quantities we need that will be modified in the predictive step but should not be touched in the in-sample part
nTablesPredictive = nTables
maxTableIndexPredictive = maxTableIndex
nTablesInRestaurantPredictive = nTablesInRestaurant
nFreeTablesPredictive = nFreeTables
freeTablesPredictive = freeTables
nPeopleAtTablePredictive = nPeopleAtTable
tablesValuesPredictive = tablesValues
tableRestaurantAllocationPredictive = tableRestaurantAllocation
nDishesPredictive = nDishes
dishesCountsPredictive = dishesCounts

indexCustomerPredictive = 1
for (indexRestaurant in 1:nRest) {
  
  for (indexCustomerRestaurant in 1:mOutOfSample[indexRestaurant]) {
    # sample the table
    indecesTablesInRestaurant = (1:maxTableIndexPredictive)[tableRestaurantAllocationPredictive==indexRestaurant]
    
    probNewTable = theta + nTablesInRestaurantPredictive[indexRestaurant]*sigma
    probs = c(nPeopleAtTablePredictive[indecesTablesInRestaurant] - sigma, probNewTable)
    
    newTableAllocation = sample(c(indecesTablesInRestaurant,-1),1,replace = F, prob = probs)
    
    if(newTableAllocation < 0) { # sample the dish
      nTablesPredictive = nTablesPredictive + 1
      nTablesInRestaurantPredictive[indexRestaurant] = nTablesInRestaurantPredictive[indexRestaurant] + 1
      if(nFreeTablesPredictive > 0) { # pick the first free table
        newTableAllocation = freeTablesPredictive[1]
        freeTablesPredictive = freeTablesPredictive[-1]
        nFreeTablesPredictive = nFreeTablesPredictive - 1
        nPeopleAtTablePredictive[newTableAllocation] = 1
      } else { # create a new table
        maxTableIndexPredictive = maxTableIndexPredictive + 1
        newTableAllocation = nTablesPredictive
        nPeopleAtTablePredictive = c(nPeopleAtTablePredictive,1)
        tablesValuesPredictive = c(tablesValuesPredictive,-1) # the dish of the table will be sampled in the next step
      }
      
      # now sample the dish and update tablesValues (we cannot have free dishes here, as we allocate sequentially, but never reallocate)
      dish = sample(1:(nDishesPredictive+1),1,replace = F,prob = c(dishesCountsPredictive-sigma0,theta0+nDishesPredictive*sigma0))
      if(dish > nDishesPredictive) {
        nDishesPredictive = nDishesPredictive + 1
        dishesCountsPredictive = c(dishesCountsPredictive,1)
      } else {dishesCountsPredictive[dish] = dishesCountsPredictive[dish] + 1}
      
      # update the table value, i.e. the dish served at the table
      tablesValuesPredictive[newTableAllocation] = dish
      
    } else{ # the sampled table is already occupied in the restaurant --> just update the relevant quantities
      nPeopleAtTablePredictive[newTableAllocation] = nPeopleAtTablePredictive[newTableAllocation] + 1
    }
    
    predictiveTableAllocationAcrossGibbs[r,indexCustomerPredictive] = newTableAllocation
    predictiveObservationDishAllocationAcrossGibbs[r,indexCustomerPredictive] = tablesValuesPredictive[newTableAllocation]
    
    indexCustomerPredictive = indexCustomerPredictive + 1
  }
}

####################### END PREDICTIVE PART




















##### TRIAL PART WITH RESAMPLING FROM PRIOR DISTRIBUTION
observationTableAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObs) # tables to which customers are allocated across gibbs
observationDishAllocationAcrossGibbs = matrix(0,nrow = nGibbsUpdates,ncol = nObs) # dishes to which customers are allocated across gibbs

for (r in 1:nGibbsUpdates) {
  indexCustomerGlobal = 1
  for (indexRestaurant in 1:nRest) {
    #nPeopleAtTableInRestaurant = nPeopleAtTable[indecesTablesInRestaurant]
    #tablesValuesInRestaurant = tablesValues[indecesTablesInRestaurant]
    # nTablesInRestaurant = length(indecesTablesInRestaurant)
    for (indexCustomerRestaurant in 1:n[indexRestaurant]) {
      indecesTablesInRestaurant = (1:nTables)[tableRestaurantAllocation==indexRestaurant]
      currentTable = tableAllocation[indexCustomerGlobal] # get the current table
      nPeopleAtTable[currentTable] = nPeopleAtTable[currentTable] - 1
      if(nPeopleAtTable[currentTable] == 0) { # free the table, and possibly the dish
        nFreeTables = nFreeTables +1
        freeTables = c(currentTable,freeTables)
        tableRestaurantAllocation[currentTable] = -1
        currentDish = tablesValues[currentTable] # get the current dish
        tablesValues[currentTable] = -1 # probably unnecessary
        dishesCounts[currentDish] = dishesCounts[currentDish] - 1 # update dishes counts
        if(dishesCounts[currentDish] == 0) { # free the dish
          nFreeDishes = nFreeDishes +1
          freeDishes = c(currentDish,freeDishes)
        }
      }
      
      # indexCurrentTableInRestaurant = min(which(nPeopleAtTable[indecesTablesInRestaurant]<0)) # find the index of the table in the restaurant
      indecesNonEmptyTablesInRestaurant = indecesTablesInRestaurant[nPeopleAtTable[indecesTablesInRestaurant]>0]
      possibleTables = c(indecesNonEmptyTablesInRestaurant,-1)
      probs = c(nPeopleAtTable[indecesNonEmptyTablesInRestaurant] - sigma,theta + length(indecesNonEmptyTablesInRestaurant)*sigma)
      newTableAllocation = sample(possibleTables,1,replace = F, prob = probs)
      
      if(newTableAllocation < 0) { # allocate to new table and sample the dish
        # allocate to new table
        if(nFreeTables > 0) { # pick the first free table
          newTableAllocation = freeTables[1]
          freeTables = freeTables[-1]
          nFreeTables = nFreeTables - 1
          nPeopleAtTable[newTableAllocation] = 1
        } else { # create a new table
          nTables = nTables + 1
          newTableAllocation = nTables
          nPeopleAtTable = c(nPeopleAtTable,1)
          tablesValues = c(tablesValues,-1) # the dish of the table will be sampled in the next step
        }
        # assign the table to the restaurant
        tableRestaurantAllocation[newTableAllocation] = indexRestaurant
        
        # sample the dish
        indecesNonEmptyDishes = (1:nDishes)[dishesCounts>0]
        possibleDishes = c((1:nDishes)[indecesNonEmptyDishes],-1) # the last dish will be either the current dish or a new one, depending if it was a unique value or not
        probs = c(dishesCounts[indecesNonEmptyDishes]-sigma0,theta0+length(indecesNonEmptyDishes)*sigma0)
        newDishAllocation = sample(possibleDishes,1,replace = F, prob = probs)
        if(newDishAllocation < 0) { # sample the new dish
          if(nFreeDishes > 0){ # pick the first free dish
            newDishAllocation = freeDishes[1]
            freeDishes = freeDishes[-1]
            nFreeDishes = nFreeDishes - 1
            dishesCounts[newDishAllocation] = 1
          } else{
            nDishes = nDishes + 1
            newDishAllocation = nDishes
            dishesCounts = c(dishesCounts,1)
          }
        } else{ # the sampled dish is already served in the franchise --> just update the relevant quantities
          dishesCounts[newDishAllocation] = dishesCounts[newDishAllocation] + 1
        }
        # update the table value, i.e. the dish served at the table
        tablesValues[newTableAllocation] = newDishAllocation
      } else{ # the sampled table is already occupied in the restaurant --> just update the relevant quantities
        nPeopleAtTable[newTableAllocation] = nPeopleAtTable[newTableAllocation] + 1
      }
      tableAllocation[indexCustomerGlobal] = newTableAllocation
      
      tableAllocationAcrossGibbs[r,indexCustomerGlobal] = newTableAllocation
      dishAllocationAcrossGibbs[r,indexCustomerGlobal] = tablesValues[newTableAllocation]
      
      indexCustomerGlobal = indexCustomerGlobal + 1
    }
  }
}


