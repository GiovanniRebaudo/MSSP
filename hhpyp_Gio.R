nByGroup = c(300, 200, 300, 100, 500) # number of obs in each group
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
nGibbsUpdates = 400
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
  if(r%%20==0){print(r)}
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
    
    # define global index to be used in the MH part to start from the right global index
    indexCustomerGlobalMH = indexCustomerGlobal
    
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
      
      #GIO
      nDishesCS = length(dishesCountsCS>0)
      
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
      currentDish = observationDishAllocation[indexCustomerGlobalMH] # get the current dish

      nTablesInRestaurant = sum(tableRestaurantAllocationMH==indexRestaurant)
      
      indecesPossibleTables = (tablesValuesMH[indecesTablesInRestaurant] == observationDishAllocation[indexCustomerGlobalMH])
      possibleTables = c(indecesTablesInRestaurant[indecesPossibleTables],-1)
      
      probNewTable = (theta + sigma*nTablesInRestaurant)/(nTables + theta0)
      nTablesServingCurrentDish = sum(tablesValues == observationDishAllocation[indexCustomerGlobalMH])
      if(nTablesServingCurrentDish > 0) {
        probNewTable = probNewTable*(nTablesServingCurrentDish - sigma0)
      }
      
      probs = c(nPeopleAtTableMH[indecesTablesInRestaurant][indecesPossibleTables] - sigma, probNewTable)
      
      newTableAllocation = sample(possibleTables,1,replace = F, prob = probs)
      
      if(newTableAllocation < 0) {
        nTablesMH = nTablesMH + 1
        if(nFreeTablesMH > 0) { # pick the first free table
          newTableAllocation = freeTablesMH[1]
          freeTablesMH = freeTablesMH[-1]
          nFreeTablesMH = nFreeTablesMH - 1
          nPeopleAtTableMH[newTableAllocation] = 1
          tablesValuesMH[newTableAllocation] = observationDishAllocationMH[indexCustomerGlobalMH]
          tableRestaurantAllocationMH[newTableAllocation] = indexRestaurant # assign table to restaurant
        } else { # create a new table
          maxTableIndexMH = maxTableIndexMH + 1
          newTableAllocationMH = maxTableIndexMH
          nPeopleAtTableMH = c(nPeopleAtTableMH,1)
          tablesValuesMH = c(tablesValuesMH,observationDishAllocationMH[indexCustomerGlobalMH])
          tableRestaurantAllocationMH = c(tableRestaurantAllocationMH,indexRestaurant) # assign table to restaurant
        }
        
      } else{ # the sampled table is already occupied in the restaurant --> just update the relevant quantities
        nPeopleAtTableMH[newTableAllocation] = nPeopleAtTableMH[newTableAllocation] + 1
      }
      
      observationTableAllocationMH[indexCustomerGlobalMH] = newTableAllocation
      
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
      observationDishAllocationMH[indexCustomerGlobalMH] = currentDish
      
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
    observationTableAllocationAcrossGibbs[r,firstIndividuals[indexGroup]:lastIndividuals[indexGroup]] = observationTableAllocation[firstIndividuals[indexGroup]:lastIndividuals[indexGroup]]
  }
  
  ### SAMPLE TABLE AND DISH FOR OUT-OF-SAMPLE OBSERVATIONS
  # THIS HAS BEEN MOVED OUT OF THE LOOP FOR THE MOMENT
}


########### Commenti Gio #######
# 1) Bisogna lanciare la funzione "CRF_PY_sample_from_prior" da file CRP_PY.R
# 2) errore: "nDishesCS" not found.
# 3) errore: "freeTables_noJ" not found.


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
