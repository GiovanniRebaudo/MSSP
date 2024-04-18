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

##### Multivariate species simulations/truth
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
niter_MH       = 1


#### Initialization Gibbs
init_all = initHSSP_fct(I_j_vec     = I_j_vec, 
                        Data_vec    = X_ji_vec,
                        tablesInit  = "equal", # "separate"
                        model       = "HPYP",
                        shape_theta = shape_theta, 
                        rate_theta  = rate_theta, 
                        a_sigma     = a_sigma, 
                        b_sigma     = b_sigma)




output = HPYP_MCMC_fct(
  nGibbsUpdates  = 1e4,
  seed           = 123,
  # seed to be fixed
  Hyperprior     = F,
  # learn hyperpar via full Bayes if  Hyperprior==T
  niter_MH       = niter_MH,
  # number of MH iterations for hyperpar update within each steps
  I_j_vec        = I_j_vec,
  Data_vec       = X_ji_vec,
  shape_theta    = shape_theta, 
  rate_theta     = rate_theta, 
  a_sigma        = a_sigma, 
  b_sigma        = b_sigma
)


nGibbsUpdates = nrow(output)
burnin        = min(nGibbsUpdates/2,500)

# Check predictive probabilities
ggplot(data = data.frame(cbind(1:nGibbsUpdates,output)), aes(x = X1)) + 
  geom_line(aes(y = X2), col=1) + 
  geom_line(aes(y = X3), col=2) +
  geom_line(aes(y = X4), col=3) +
  labs(x="iter", y = "") 

# Choose the optimal arm
which.max(colMeans(output[burnin:nGibbsUpdates,]))

colMeans(output[burnin:nGibbsUpdates,])






