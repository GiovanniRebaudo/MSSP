# Functions used in the codes accompanying 
# "Multivariate Species Sampling Process (MSSP)"

# Functions to compute the urn schemes
urn_DP_norm <- function(freq_minus,theta_PY){
  return(c(freq_minus,theta_PY)/(sum(freq_minus)+theta_PY))
}

urn_PY <- function(freq_minus,theta_PY,sigma_PY){
  #  K<-length(freq_minus)
  return(c(freq_minus-sigma_PY,theta_PY+K*sigma_PY))
}

urn_Dir <- function(freq_minus,beta_DM,K_max){
  #  K<-length(freq_minus)
  return(c(freq_minus+beta_DM,beta_DM*(K_max-K)*(K_max>K)))
}

urn_GN_norm <- function(freq_minus,gamma_GN){
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