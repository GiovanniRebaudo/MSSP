# Codes accompanying "Entropy Regularization in Probabilistic Clustering"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.14
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(salso)   
library(ggplot2) 
library(plyr)
# Load functions
source("utils.R")
source("MSSP_fcts.R")
Save_Plot = FALSE

set.seed(123)

# ##### Multivariate species simulations/truth
# J           = 3 # Number of populations
# I_j_vec     = rep(100,J) 
# cum_I_j_vec = cumsum(I_j_vec)
# # I_j_vec = (I_1, ...,I_J) vector of sample size in different population
# n           = sum(I_j_vec) # tot number of observations
# 
# # X_ji_vec = integer(n) 
# # Values of all observations (X_{1,1},...,X_{1,I_1},...,X_{J,1},...,X_{J,I_J})
# X_ji_vec = c(rep(1,I_j_vec[1]), # X_1i_vec = (X_{1,1},...,X_{1,I_1})
#              rep(c(rep(2,15),rep(3,20),rep(1,15)),2), # X_2i_vec
#              rep(1:10,10)) # X_3i_vec
# 
# 
# Xstar_d_vec = unique(X_ji_vec) # observed dishes (dishes=species)
# D           = length(Xstar_d_vec) # overall number of dishes
# # Check 
# if (n != length(X_ji_vec)){print("error: n != length(X_ji_vec)"); stop()}
# if (sum(Xstar_d_vec != 1:D)){
#   print("error: labels of dish are not ordered"); stop()}
# if (J != length(I_j_vec)){print("error: J != length(I_j_vec)"); stop()}


# generate true pmf
truth_par = c(rep(1.3, 2), rep(2, 6))
truth     = generate_zipf(param = truth_par, 
                     tot_species = 3000, j_species = 2500, seed = 0)
# how many new sample? 
init_samples = 30
new_samples  = 300

# sample data
J           = length(truth_par) 
I_j_vec     = rep(init_samples, J)
n           = sum(I_j_vec)
data = c()
dataNewLs = list()
for (j in 1:J){
  dataj =  sample_from_pop(j, truth, size = init_samples+new_samples)
  data = c(data, dataj[1:init_samples])
  dataNewLs[[j]] = dataj[(init_samples+1):new_samples]
}

# Rerorder species labels in order of arrival by group
X_ji_vec = as.integer(factor(data, levels = unique(data)))
for (j in 1:J){
  dataNewLs[[j]] = plyr::mapvalues(dataNewLs[[j]], 
                                   from = unique(data), 
                                   to   = unique(X_ji_vec))
}



### Preliminaries and data summaries
# Truth of in sample species
table(X_ji_vec) # overall species freq. (n_{.,d})_{d=1}^{D}

# Numbers of different species within pop
K_j_vec_fct(I_j_vec=I_j_vec, Data_vec=X_ji_vec)

emp_pEPPF_un =emp_pEPPF_un_fct(I_j_vec=I_j_vec, Data_vec=X_ji_vec)

emp_pEPPF_un
round(emp_pEPPF_un/n,2) # empirical pEPPF


shape_theta    = 1
rate_theta     = 1
a_sigma        = 1 
b_sigma        = 1
niter_MH       = 5

# save more quantities in MCMC for debugging and convergence checks
species_discovered = logical(new_samples)

if(F){
#### Initialization Gibbs
init_all = initHSSP_fct(I_j_vec     = I_j_vec, 
                        Data_vec    = X_ji_vec,
                        tablesInit  = "equal", # "separate"
                        model       = "HPYP",
                        shape_theta = shape_theta, 
                        rate_theta  = rate_theta, 
                        a_sigma     = a_sigma, 
                        b_sigma     = b_sigma)


output = "all"# c("prob_new", "prob and last", "all")

# Run a short mcmc with fixed hyper par to better initialize the tables
init_all = HPYP_MCMC_fct(
  nGibbsUpdates  = 2e3,
  seed           = 123,
  # seed to be fixed
  Hyperprior     = F,
  # learn hyperpar via full Bayes if  Hyperprior==T
  niter_MH       = 1,
  # number of MH iterations for hyperpar update within each steps
  I_j_vec        = I_j_vec,
  Data_vec       = X_ji_vec,
  shape_theta    = shape_theta, 
  rate_theta     = rate_theta, 
  a_sigma        = a_sigma, 
  b_sigma        = b_sigma,
  output         = "prob and last"
)

Hyperprior = T

# Run MCMC
out = HPYP_MCMC_fct(
  nGibbsUpdates  = 1e4,
  seed           = 123,
  # seed to be fixed
  Hyperprior     = Hyperprior,
  # learn hyperpar via full Bayes if Hyperprior==T
  niter_MH       = niter_MH,
  # number of MH iterations for hyperpar update within each steps
  I_j_vec        = I_j_vec,
  Data_vec       = X_ji_vec,
  shape_theta    = shape_theta, 
  rate_theta     = rate_theta, 
  a_sigma        = a_sigma, 
  b_sigma        = b_sigma,
  output         = output
)



# Additional checks and debugging
if(output=="prob_new"){
  output_prob = out
} else {
  output_prob = out$prob_new_species
}


nGibbsUpdates = nrow(output_prob)
iter_considered = 1:nGibbsUpdates # iter_considered = burnin:nGibbsUpdates

# Check predictive probabilities
ggplot(data = data.frame(cbind(iter_considered, output_prob)), 
       aes(x = iter_considered)) + 
  geom_line(aes(y = V2), col=1) + 
  geom_line(aes(y = V3), col=2) +
  geom_line(aes(y = V4), col=3) +
  labs(x="iter", y = "prob new") 


burnin        = min(nGibbsUpdates/2,500)
# Choose the optimal arm
which.max(colMeans(output_prob[burnin:nGibbsUpdates,]))
colMeans(output_prob[burnin:nGibbsUpdates,])

if(output=="all"){
  P = ggplot(data = data.frame(cbind(iter_considered, output_prob)), 
         aes(x = iter_considered)) + 
    geom_line(aes(y = V2), col=1) + 
    geom_line(aes(y = V3), col=2) +
    geom_line(aes(y = V4), col=3) +
    labs(x="iter", y = "prob new") 
  
  print(P)

  
  P = ggplot(data = 
           data.frame(cbind(iter_considered, 
                            out$theta_vecAcrossGibbs[iter_considered,])), 
         aes(x = iter_considered)) + 
    geom_line(aes(y = V2), col=1) + 
    geom_line(aes(y = V3), col=2) +
    geom_line(aes(y = V4), col=3) +
    labs(x="iter", y = "thetaj") 
  
  print(P)
  
  P = ggplot(data = 
           data.frame(cbind(iter_considered, 
                            out$sigma_vecAcrossGibbs[iter_considered,])), 
         aes(x = iter_considered)) + 
    geom_line(aes(y = V2), col=1) + 
    geom_line(aes(y = V3), col=2) +
    geom_line(aes(y = V4), col=3) +
    labs(x="iter", y = "sigmaj") 
  
  print(P)
  
  P = ggplot(data = 
           data.frame(cbind(iter_considered, 
                            out$theta0AcrossGibbs[iter_considered])), 
         aes(x = iter_considered)) + 
    geom_line(aes(y = V2), col=1) + 
    labs(x="iter", y = "theta0") 
  
  print(P)
  
  P = ggplot(data = 
           data.frame(cbind(iter_considered, 
                            out$sigma0AcrossGibbs[iter_considered])), 
         aes(x = iter_considered)) + 
    geom_line(aes(y = V2), col=1) + 
    labs(x="iter", y = "sigma0") 
  
  print(P)
  
  P = ggplot(data = 
           data.frame(cbind(iter_considered, 
                      out$nTablesInRestaurantAcrossGibbs[iter_considered,])), 
         aes(x = iter_considered)) + 
    geom_line(aes(y = V2), col=1) + 
    geom_line(aes(y = V3), col=2) +
    geom_line(aes(y = V4), col=3) +
    labs(x="iter", y = "ntablesj") 
  
  print(P)
}


if(output=="all" && Hyperprior){
  Move_sigma_j_out = out$Move_sigma_j_out
  print(rowMeans(Move_sigma_j_out))
  Move_theta_j_out = out$Move_theta_j_out
  print(rowMeans(Move_theta_j_out))
  # out$Prop_sd_log_theta_j
  # out$Prop_sd_logit_sig_j
}

}


####### MAB
#### Initialization Gibbs
init_all = initHSSP_fct(I_j_vec     = I_j_vec, 
                        Data_vec    = X_ji_vec,
                        tablesInit  = "equal", # "separate"
                        model       = "HPYP",
                        shape_theta = shape_theta, 
                        rate_theta  = rate_theta, 
                        a_sigma     = a_sigma, 
                        b_sigma     = b_sigma)

# Run a short mcmc with fixed hyper par to better initialize the tables
init_all = HPYP_MCMC_fct(
  nGibbsUpdates  = 1e3,
  seed           = 123,
  # seed to be fixed
  Hyperprior     = F,
  # learn hyperpar via full Bayes if  Hyperprior==T
  niter_MH       = 1,
  # number of MH iterations for hyperpar update within each steps
  I_j_vec        = I_j_vec,
  Data_vec       = X_ji_vec,
  shape_theta    = shape_theta, 
  rate_theta     = rate_theta, 
  a_sigma        = a_sigma, 
  b_sigma        = b_sigma,
  output         = "prob and last"
)

# Run MCMC
out = HPYP_MCMC_fct(
  nGibbsUpdates  = 4e3,
  seed           = 123,
  # seed to be fixed
  Hyperprior     = T,
  # learn hyperpar via full Bayes if Hyperprior==T
  niter_MH       = niter_MH,
  # number of MH iterations for hyperpar update within each steps
  I_j_vec        = I_j_vec,
  Data_vec       = X_ji_vec,
  shape_theta    = shape_theta, 
  rate_theta     = rate_theta, 
  a_sigma        = a_sigma, 
  b_sigma        = b_sigma,
  output         = "prob and last"
)



for (iter_new in 1:new_samples){
  # save 
  prob_new_species = out$prob_new_species
  dishAllocation   = out$dishAllocation
  
  #
  nGibbsUpd = nrow(prob_new_species)
  burnin    = min(1000, nGibbsUpd/2)
  # Choose optimal arm
  newj = which.max(colMeans(prob_new_species[burnin:nGibbsUpd,]))
  # Pick new obs
  newObs = dataNewLs[[newj]][iter_new]
  
  # Check if a new species is discovered
  species_discovered[iter_new] = !(newObs %in% dishAllocation)
  
  init_all = initSeqHSSP_fct(newPop = newj,
                             newDataPoint = newObs)
  
  if(iter_new<new_samples){
    # Run MCMC
    out = HPYP_MCMC_fct(
      nGibbsUpdates  = 4e2,
      seed           = 123,
      # seed to be fixed
      Hyperprior     = T,
      # learn hyperpar via full Bayes if Hyperprior==T
      niter_MH       = 3,
      # number of MH iterations for hyperpar update within each steps
      I_j_vec        = init_all$I_j_vec,
      Data_vec       = init_all$dishAllocation,
      shape_theta    = shape_theta, 
      rate_theta     = rate_theta, 
      a_sigma        = a_sigma, 
      b_sigma        = b_sigma,
      output         = "prob and last"
    )
  }
}





