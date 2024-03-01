rm(list = ls())
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
      nDishes = nDishes + 1
      dishesCounts = 1
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
    tableRestaurantAllocation = c(tableRestaurantAllocation,indexRestaurant) # allocate the table to the restaurant
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
        tableRestaurantAllocation = c(tableRestaurantAllocation,indexRestaurant) # allocate the table to the restaurant
        tableCountsInRestaurant = c(tableCountsInRestaurant,1)
        nPeopleAtTable = c(nPeopleAtTable,1)
        # sample the dish served at the new table
        dish = sample(1:(nDishes+1),1,replace = F,prob = c(dishesCounts-sigma0,theta0+nDishes*sigma0))
        if(dish == nDishes+1) {
          nDishes = nDishes + 1
          dishesCounts = c(dishesCounts,1)
        } else {dishesCounts[dish] = dishesCounts[dish] + 1}
        tablesValues = c(tablesValues,dish)
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


nByGroup = c(50, 50, 50, 50, 50)*2 # c(300, 300, 300, 300, 400) # number of obs in each group
firstIndividuals = 1 + cumsum(c(0,nByGroup[-length(nByGroup)]))
lastIndividuals = cumsum(nByGroup)
alpha = 0.5
gamma = 10 # 10
sigma0 = 0.5
theta0 = 10 # 10
sigma = 0.5
theta = 10 # 10

# simulate observations
set.seed(123)
# groupRestaurantAllocation = c(1,1,2,2,3) # allocation of groups by restaurant
# observationRestaurantAllocation = rep(groupRestaurantAllocation,times = n) # allocation of obs to restaurant
n = c(nByGroup[1]+nByGroup[2],nByGroup[3]+nByGroup[4],nByGroup[5])
crf_py = CRF_PY_sample_from_prior(n,sigma0,theta0,sigma,theta)
crf_py

# define quatities for gibbs sampling
# INITIALIZATION TO ALL DIFFERENT TABLES AND DIFFERENT RESTAURANTS
nGibbsUpdates = 1000 # 300
nObs = sum(nByGroup)
nGroups = length(nByGroup)
nTables = nObs # number of occupied tables in the franchise
maxTableIndex = nObs # max table index (nTables+nFreeTables=maxTableIndex)
#nTablesInGroup = nByGroup
# nTablesInRestaurant = nByGroup
observationGroupAllocation = rep(1:nGroups, times = nByGroup)
observationDishAllocation = crf_py$dishAllocation
observationTableAllocation = 1:nObs
groupRestaurantAllocation = 1:nGroups
observationRestaurantAllocation = groupRestaurantAllocation[observationGroupAllocation]
nGroupsInRestaurant = as.integer(table(groupRestaurantAllocation)) # rep(1,nRest)
nRest = length(table(groupRestaurantAllocation))
maxRestIndex = nRest
tablesValues = observationDishAllocation
tableRestaurantAllocation = observationRestaurantAllocation
nPeopleAtTable = rep(1,nObs)


# trial part with manually-generated observations
# observationDishAllocation[firstIndividuals[1]:lastIndividuals[1]] = c(rep(1,nByGroup[1]/2),rep(3,nByGroup[1]/2))
# observationDishAllocation[firstIndividuals[2]:lastIndividuals[2]] = c(rep(2,nByGroup[2]/2),rep(3,nByGroup[2]/2))
# observationDishAllocation[firstIndividuals[3]:lastIndividuals[3]] = c(rep(1,nByGroup[3]/2),rep(2,nByGroup[3]/2))
# observationDishAllocation[firstIndividuals[4]:lastIndividuals[4]] = c(rep(1,nByGroup[4]/2),rep(3,nByGroup[4]/2))
# observationDishAllocation[firstIndividuals[5]:lastIndividuals[5]] = c(rep(2,nByGroup[2]/2),rep(4,nByGroup[2]/2)) # rep(2,nByGroup[5])



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


############### PRELIMINARY PLOT
# hist(x = crf_py$dishAllocation, breaks = as.integer(names(table(crf_py$dishAllocation))))
# 
# par(mfrow = c(3,2))
# for (i in 1:nGroups) {
#   hist(x = crf_py$dishAllocation[observationGroupAllocation==i], breaks = as.integer(names(table(crf_py$dishAllocation))))
# }
# 
# par(mfrow = c(1,1))


# to delete
r=1
indexCustomerGlobal=1
indexGroup=1

### trial

nByGroup = c(50, 50, 50, 50, 100)*2

nGibbsUpdates = 1000 # 300
nObs = sum(nByGroup)
nGroups = length(nByGroup)
nTables = nObs # number of occupied tables in the franchise
maxTableIndex = nObs # max table index (nTables+nFreeTables=maxTableIndex)
#nTablesInGroup = nByGroup
# nTablesInRestaurant = nByGroup
observationGroupAllocation = rep(1:nGroups, times = nByGroup)

observationDishAllocation          = integer(nObs)
observationDishAllocation[1:200]   = rep(1,100);
observationDishAllocation[201:400] = rep(c(rep(2,30),rep(3,40),rep(1,30)),2)
observationDishAllocation[401:600] = rep(1:100,2)

# observationDishAllocation = crf_py$dishAllocation


observationTableAllocation = 1:nObs
groupRestaurantAllocation  = 1:nGroups
observationRestaurantAllocation = groupRestaurantAllocation[observationGroupAllocation]
nGroupsInRestaurant = as.integer(table(groupRestaurantAllocation)) # rep(1,nRest)
nRest = length(table(groupRestaurantAllocation))
maxRestIndex = nRest
tablesValues = observationDishAllocation
tableRestaurantAllocation = observationRestaurantAllocation
nPeopleAtTable = rep(1,nObs)


dishesCounts = as.vector(table(observationDishAllocation)) # how many people are eating a certain dish
nDishes = length(dishesCounts)



set.seed(123)
### end trial

for (r in 1:nGibbsUpdates) {
  ### ALLOCATE IN-SAMPLE OBSERVATIONS TO TABLES
  if(r%%20==0){print(r)}
  indexCustomerGlobal = 1
  for (indexGroup in 1:nGroups) {
    indexRestaurant = groupRestaurantAllocation[indexGroup]
    
    ##### compute quantities without group J
    observationRestaurantAllocation_noJ = observationRestaurantAllocation
    individualsInCurrentGroup = firstIndividuals[indexGroup]:lastIndividuals[indexGroup]
    observationRestaurantAllocation_noJ[individualsInCurrentGroup] = -1

    # nPeopleInRestaurant_noJ = sum(observationRestaurantAllocation_noJ==indexRestaurant)

    dishesCountsInGroup = as.vector(table(factor(observationDishAllocation[individualsInCurrentGroup],
                                                 levels = 1:nDishes))) # number of people eating dish d in group J
    dishesCounts_noJ = dishesCounts - dishesCountsInGroup
    nDishes_noJ = nDishes - sum(dishesCounts_noJ==0)
    ##### end computation of quantities without group J
    
    ##### initialize quantities for computation of acceptance probability for the "current" state
    # nPeopleInRestaurantCS = nPeopleInRestaurant_noJ
    dishesCountsCS = dishesCounts_noJ
    nDishesCS = nDishes_noJ
    ##### end initialization quantities for acceptance prob
    
    # define global index to be used in the MH part to start from the right global index
    indexCustomerGlobalMH = indexCustomerGlobal
    
    logP_CS = 0 # log of P[X_ji=x|____S] to be used in MH
    for (indexCustomerGroup in 1:nByGroup[indexGroup]) {
      currentTable = observationTableAllocation[indexCustomerGlobal] # get the current table
      currentDish = observationDishAllocation[indexCustomerGlobal] # get the current dish (this definition could be avoided)
      nPeopleAtTable[currentTable] = nPeopleAtTable[currentTable] - 1
      
      if(nPeopleAtTable[currentTable] == 0) { # free the table
        nFreeTables = nFreeTables + 1
        freeTables = c(currentTable,freeTables)
        tableRestaurantAllocation[currentTable] = -1
        nTables = nTables - 1
        tablesValues[currentTable] = -1
      }
      
      nTablesInRestaurant = sum(tableRestaurantAllocation==indexRestaurant)
      
      indecesTablesInRestaurant = (1:maxTableIndex)[tableRestaurantAllocation==indexRestaurant] # indeces of tables in the restaurant
      indecesPossibleTables = (tablesValues[indecesTablesInRestaurant] == currentDish) # tables in the restaurant serving the current dish
      possibleTables = c(indecesTablesInRestaurant[indecesPossibleTables],-1)
      
      probNewTable = (theta + sigma*nTablesInRestaurant)/(nTables + theta0)
      nTablesServingCurrentDish = sum(tablesValues == currentDish)
      if(nTablesServingCurrentDish > 0) {
        probNewTable = probNewTable*(nTablesServingCurrentDish - sigma0)
      }
      
      probs = c(nPeopleAtTable[indecesTablesInRestaurant][indecesPossibleTables] - sigma, probNewTable)
      
      newTableAllocation = sample(possibleTables,1,replace = F, prob = probs)
      
      if(newTableAllocation < 0) {
        ##### COMPUTE THE CONTRIBUTION PART OF THE LOG ACCEPTANCE PROBABILITY
        individualsToRemove = c((firstIndividuals[indexGroup]-1+indexCustomerGroup):lastIndividuals[indexGroup])
        nPeopleAtTableCS = as.integer(table(factor(observationTableAllocation[-individualsToRemove], levels = 1:maxTableIndex))) # nPeopleAtTable - as.integer(table(factor(observationTableAllocation[individualsToRemove], levels = 1:maxTableIndex))) # as.integer(table(factor(observationTableAllocation_noJ, levels = 1:maxTableIndex))) # this is wrong since it removes twice the current observation from the count of people at table 

        nTablesCS = sum(nPeopleAtTableCS>0)
        
        # compute nTablesInRestaurantCS
        # freeTablesCS = (1:maxTableIndex)[nPeopleAtTableCS == 0]
        # tableRestaurantAllocationCS = tableRestaurantAllocation
        # tableRestaurantAllocationCS[freeTablesCS] = -1
        # nTablesInRestaurantCS = sum(tableRestaurantAllocationCS==indexRestaurant)
        
        # nTablesInRestaurantCS = 
        #   length(table(observationTableAllocation[-individualsToRemove][observationRestaurantAllocation[-individualsToRemove]==indexRestaurant]))
        # 
        
        # compute log of P[X_ji=x|____]
        if(dishesCountsCS[currentDish] == 0){ # X_ji = "new"
          logP_CS = logP_CS +
            log(theta0 + nDishesCS*sigma0) -
            log(theta0 + nTablesCS)
          nDishesCS = nDishesCS + 1
        } else{
          # compute number of tables serving the current dish in the whole franchise, excluding individualsToRemove
          # nTablesServingCurrentDishCS = length(table(observationTableAllocation[-individualsToRemove][observationDishAllocation[-individualsToRemove] == currentDish])) 
          # alternative computation:
          freeTablesCS = (1:maxTableIndex)[nPeopleAtTableCS == 0] # which(nPeopleAtTableCS == 0)
          tablesValuesCS = tablesValues
          tablesValuesCS[freeTablesCS] = -1
          nTablesServingCurrentDishCS = sum(tablesValuesCS==currentDish)
          
          # get the tables in current restaurant serving the current dish (excluding future observations)
          # freeTablesCS = (1:maxTableIndex)[nPeopleAtTableCS == 0] # which(nPeopleAtTableCS == 0)
          # tableRestaurantAllocationCS = tableRestaurantAllocation
          # tableRestaurantAllocationCS[freeTablesCS] = -1
          # tablesInRestaurantServingCurrentDishCS = (1:maxTableIndex)[((tableRestaurantAllocationCS == indexRestaurant)&(tablesValues == currentDish))]
          
          # tablesValuesCS = tablesValues
          # tablesValuesCS[freeTablesCS] = -1
          
          #nTablesInRestaurantServingCurrentDishCS = length(tablesInRestaurantServingCurrentDish)
          
          logP_CS = logP_CS +
            log(nTablesServingCurrentDishCS - sigma0) -
            log(theta0 + nTablesCS)
          
        }
        
        ### update quantities of "current state"
        # nPeopleInRestaurantCS = nPeopleInRestaurantCS + 1
        dishesCountsCS[currentDish] = dishesCountsCS[currentDish] + 1
        
        ##### UPDATE GIBBS SAMPLING QUANTITIES
        nTables = nTables + 1
        if(nFreeTables > 0) { # pick the first free table
          newTableAllocation = freeTables[1]
          freeTables = freeTables[-1]
          nFreeTables = nFreeTables - 1
          nPeopleAtTable[newTableAllocation] = 1
          tablesValues[newTableAllocation] = currentDish
          tableRestaurantAllocation[newTableAllocation] = indexRestaurant # assign table to restaurant
        } else { # create a new table
          maxTableIndex = maxTableIndex + 1
          newTableAllocation = maxTableIndex
          nPeopleAtTable = c(nPeopleAtTable,1)
          tablesValues = c(tablesValues,currentDish)
          tableRestaurantAllocation = c(tableRestaurantAllocation,indexRestaurant) # assign table to restaurant
        }
        
      } else{ # the sampled table is already occupied in the restaurant --> just update the relevant quantities
        nPeopleAtTable[newTableAllocation] = nPeopleAtTable[newTableAllocation] + 1
      }
      
      observationTableAllocation[indexCustomerGlobal] = newTableAllocation
      
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
        groupRestaurantAllocationMH[indexGroup] = newRestaurantAllocation # assign group to restaurant
      }
      
    } else{ # the sampled restaurants contains already some groups --> just update the relevant quantities
      groupRestaurantAllocationMH[indexGroup] = newRestaurantAllocation # assign group to restaurant
      nGroupsInRestaurantMH[newRestaurantAllocation] = nGroupsInRestaurantMH[newRestaurantAllocation] + 1
    }
    
    #### sample table allocation
    
    # initialize MH quantities excluding current group j
    maxTableIndexMH = maxTableIndex
    
    observationTableAllocationMH = observationTableAllocation
    observationTableAllocationMH[individualsInCurrentGroup] = -1
    observationRestaurantAllocationMH = observationRestaurantAllocation_noJ # computed above
    observationDishAllocationMH = observationDishAllocation
    observationDishAllocationMH[individualsInCurrentGroup] = -1
    nPeopleAtTableMH = as.integer(table(factor(observationTableAllocation[-individualsInCurrentGroup], levels = 1:maxTableIndex))) # nPeopleAtTable - as.integer(table(factor(observationTableAllocation[individualsInCurrentGroup], levels = 1:maxTableIndex))) # as.integer(table(factor(observationTableAllocation_noJ, levels = 1:maxTableIndex)))
    nTablesMH = sum(nPeopleAtTableMH>0)
    
    freeTablesMH = (1:maxTableIndexMH)[nPeopleAtTableMH == 0]
    nFreeTablesMH = length(freeTablesMH)
    
    tablesValuesMH = tablesValues
    tablesValuesMH[freeTablesMH] = -1
    
    tableRestaurantAllocationMH = tableRestaurantAllocation
    tableRestaurantAllocationMH[freeTablesMH] = -1
    
    # nPeopleInRestaurantMH = nPeopleInRestaurant_noJ # computed above
    
    dishesCountsMH = dishesCounts_noJ
    
    nDishesMH = nDishes - sum(dishesCountsMH==0)
    
    # # nTablesServingCurrentDish = sum(tablesValues_noJ == observationDishAllocation[indexCustomerGlobal])
    ##### end computation of quantities without group J
    
    
    # #### UPDATE TABLE CONFIGURATION IN THE CURRENT GROUP (analogous to the sampling from the full conditionals)
    logP_MH = 0 # log of P[X_ji=x|____ S'] to be used in MH
    for (indexCustomerGroup in 1:nByGroup[indexGroup]) { # this loop should be joined with the loop for the full conditionals!!!!!
      indecesTablesInRestaurant = (1:maxTableIndexMH)[tableRestaurantAllocationMH==newRestaurantAllocation]
      currentDish = observationDishAllocation[indexCustomerGlobalMH] # get the current dish
      observationRestaurantAllocationMH[indexCustomerGlobalMH] = newRestaurantAllocation

      nTablesInRestaurant = sum(tableRestaurantAllocationMH==newRestaurantAllocation)
      
      indecesPossibleTables = (tablesValuesMH[indecesTablesInRestaurant] == observationDishAllocation[indexCustomerGlobalMH])
      possibleTables = c(indecesTablesInRestaurant[indecesPossibleTables],-1)
      
      probNewTable = (theta + sigma*nTablesInRestaurant)/(nTablesMH + theta0)
      nTablesServingCurrentDishMH = sum(tablesValuesMH == observationDishAllocation[indexCustomerGlobalMH])
      if(nTablesServingCurrentDishMH > 0) {
        probNewTable = probNewTable*(nTablesServingCurrentDishMH - sigma0)
      }
      
      probs = c(nPeopleAtTableMH[indecesTablesInRestaurant][indecesPossibleTables] - sigma, probNewTable)
      
      newTableAllocation = sample(possibleTables,1,replace = F, prob = probs)
      
      if(newTableAllocation < 0) {
        # compute log of P[X_ji=x|____S']
        if(dishesCountsMH[currentDish] == 0){ # X_ji = "new"
          logP_MH = logP_MH +
            log(theta0 + nDishesMH*sigma0) -
            log(theta0 + nTablesMH)
        } else{
          # compute number of tables serving the current dish in the whole franchise, excluding future individuals
          nTablesServingCurrentDishMH = length(table(observationTableAllocationMH[observationDishAllocationMH == currentDish]))
          
          # get the tables in current restaurant serving the current dish (excluding future observations)
          # tablesInRestaurantServingCurrentDish = (1:maxTableIndexMH)[((tableRestaurantAllocationMH == newRestaurantAllocation)&(tablesValuesMH == currentDish))]
          
          logP_MH = logP_MH +
            log(nTablesServingCurrentDishMH - sigma0) -
            log(theta0 + nTablesMH)
        }
        
        ### update quantities of MH
        # nPeopleInRestaurantMH = nPeopleInRestaurantMH + 1
        dishesCountsMH[currentDish] = dishesCountsMH[currentDish] + 1
        observationDishAllocationMH[indexCustomerGlobalMH] = currentDish
        
        nTablesMH = nTablesMH + 1
        if(nFreeTablesMH > 0) { # pick the first free table
          newTableAllocation = freeTablesMH[1]
          freeTablesMH = freeTablesMH[-1]
          nFreeTablesMH = nFreeTablesMH - 1
          nPeopleAtTableMH[newTableAllocation] = 1
          tablesValuesMH[newTableAllocation] = observationDishAllocationMH[indexCustomerGlobalMH]
          tableRestaurantAllocationMH[newTableAllocation] = newRestaurantAllocation # assign table to restaurant
        } else { # create a new table
          maxTableIndexMH = maxTableIndexMH + 1
          newTableAllocationMH = maxTableIndexMH
          nPeopleAtTableMH = c(nPeopleAtTableMH,1)
          tablesValuesMH = c(tablesValuesMH,observationDishAllocationMH[indexCustomerGlobalMH])
          tableRestaurantAllocationMH = c(tableRestaurantAllocationMH,newRestaurantAllocation) # assign table to restaurant
        }
        
      } else{ # the sampled table is already occupied in the restaurant --> just update the relevant quantities
        nPeopleAtTableMH[newTableAllocation] = nPeopleAtTableMH[newTableAllocation] + 1
      }
      
      observationTableAllocationMH[indexCustomerGlobalMH] = newTableAllocation
      
      # update indexCustomerGlobalMH 
      indexCustomerGlobalMH = indexCustomerGlobalMH + 1
    }
    pAccept = min(exp(logP_MH - logP_CS),1)
    accept = runif(1) < pAccept
    # print(pAccept)
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
    observationTableAllocationAcrossGibbs[r,firstIndividuals[indexGroup]:lastIndividuals[indexGroup]] = 
      observationTableAllocation[firstIndividuals[indexGroup]:lastIndividuals[indexGroup]]
  }
  
  ### SAMPLE TABLE AND DISH FOR OUT-OF-SAMPLE OBSERVATIONS
  # THIS HAS BEEN MOVED OUT OF THE LOOP FOR THE MOMENT
}


print(groupRestaurantAllocationAcrossGibbs[(nGibbsUpdates-99):nGibbsUpdates,])

summary(groupRestaurantAllocationAcrossGibbs)

library(salso)

salso(groupRestaurantAllocationAcrossGibbs[-200,])

# groupRestaurantAllocationAcrossGibbs[201:300,]

