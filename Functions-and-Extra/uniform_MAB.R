#uniform_MAB: Multi arm bandit with Uniform random choices,
#             returns the cumulative number of species discovered

uniform_MAB <- function(data,
                        init_samples = 30, new_samples = 300, 
                        seed = 0){
  ## returns the cumulative number of species discovered
  ##inputs: 
  ##  data = observations
  ##  init_samples = number of starting observations for the MAB
  ##  new_samples = number of sampling step of the MAB
  ##  seed 
  
  cat("\n Uniform \n")
  if(ncol(data)<(init_samples+new_samples)){
    cat("not enough data provided for", init_samples, "initial samples and", 
        new_samples, "sequential sampling steps")
    cat("\n data should be a matrix with nrow equals to the number of pops")
  }
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = new_samples, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  set.seed(seed)
  J = nrow(data) #tot number of populations
  
  species_discovered = rep(0, new_samples) #vector to save the num of discoveries
  
  X = matrix(NA,nrow = J, 
             ncol = init_samples + new_samples) #matrix of observations X[j,i]is X_{j,i}
  
  #sample initial observations
  for(j in 1:J){
    X[j, 1:init_samples] = data[j, 1:init_samples] 
  }
  
  I = rep(init_samples, J) #initial sample sizes
  
  for(newobs in 1:new_samples){ #MAB for cycle
    
    #estimate prediction prob
    where = sample(J,1)
    
    #sample a new observation
    x = data[where, I[where] + 1] 
    
    #check if species is new
    species_discovered[newobs] = !(x %in% X)
    
    X[where, I[where]+1] = x
    I[where] = I[where] + 1
    
    
    setTxtProgressBar(pb, newobs)
  }#for MAB

  return(list(discoveries = cumsum(species_discovered)))
}#function
