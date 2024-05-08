# Multiarmed bandit for species discovery via mSSP
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
source("mSSPmab.R")
source("MAB_HPYP_fct.R")


###############which true pmf? 
J = 8
ordered = FALSE

if(!ordered){
  #generate true pmf
  pmfs = generate_zipf(param = c(rep(1.3, 4), rep(2, 4)), 
                       tot_species = 3000, j_species = 2500, seed = 0)
}else{
  pmfs = generate_zipf_reorder(param = c(rep(1.3, 4), rep(2, 4)), 
                               tot_species = 3000, j_species = 2500, seed = 0)
}


###############how many initial and new samples? 
init_samples = 30 #per each pop
new_samples = 300

############### how many replicas?
tot_replica = 10

###############initialize for more replicas
results_HPY      = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_HDP      = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_random   = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_oracle   = matrix(NA, nrow = new_samples, ncol = tot_replica)
est_prob_new_HPY = vector("list", tot_replica)
est_prob_new_HDP = vector("list", tot_replica)


############### Gibbs samplers
for(replica in 1:tot_replica){
  cat("\nReplica", replica, "out of", tot_replica, "\n")
  
  ############### Sample observations for fair comparison of methods
  X = sample_from_pop_all(truth = pmfs, size = init_samples + new_samples,
                          seed = replica, verbose = FALSE)
  
  if(F){
    #solve MAB decisions via HPY
    results_HPY_temp = HPY_MAB(data = X,
                               a_alpha = 1, b_alpha = 1,
                               init_samples = init_samples,
                               new_samples = new_samples,
                               a_sigma = 1, b_sigma = 2,
                               burnin = 100, iters = 200, seed = 0,
                               niter_MH = 10, ada_step = 10,
                               ada_thresh = 0.44, r_ada_input = 0)
    results_HPY[,replica] = results_HPY_temp$discoveries
    est_prob_new_HPY[[replica]] = results_HPY_temp$probs
  }
  
  #solve MAB decisions via HDP
  results_HDP_temp = HDP_MAB(data = X,
                             a_alpha = 1, b_alpha = 1,
                             init_samples = init_samples,
                             new_samples = new_samples,
                             burnin = 100, iters = 200, seed = 0,
                             niter_MH = 10, ada_step = 10,
                             ada_thresh = 0.44, r_ada_input = 0)
  results_HDP[,replica] = results_HDP_temp$discoveries
  est_prob_new_HDP[[replica]] = results_HDP_temp$probs
  
  #solve MAB decision via uniform
  results_random_temp = uniform_MAB(data = X, new_samples = new_samples, 
                                    seed = 0)
  results_random[,replica] = results_random_temp$discoveries
  
  #solve MAB decision via oracle
  results_oracle_temp = oracle_MAB(data = X, pmfs = pmfs)
  results_oracle[,replica] = results_oracle_temp$discoveries
}

if(F){
  result_HPY_mean      = rowMeans(results_HPY)
} else {
  load("./Data-and-Results/result_HPY_mean.RData")
}

result_HDP_mean      = rowMeans(results_HDP)
results_random_mean  = rowMeans(results_random)
results_oracle_mean  = rowMeans(results_oracle)

if(F){
  save(result_HPY_mean,     file="./Data-and-Results/result_HPY_mean.RData")
  save(result_HDP_mean,     file="./Data-and-Results/result_HDP_mean.RData")
  save(results_random_mean, file="./Data-and-Results/results_random_mean.RData")
  save(results_oracle_mean, file="./Data-and-Results/results_oracle_mean.RData")
}


################################################################################
# Plot results 

# prepare data matrix
names = c("HPY", "HDP", "Uniform", "Oracle")
num_model_to_compare = length(names)
model = c()
for(mm in 1:num_model_to_compare){
  model = c(model, rep(names[mm], new_samples))
}
data_plot <- data.frame(
  time = rep(1:new_samples, num_model_to_compare),
  model = model,
  value = c(result_HPY_mean, result_HDP_mean,
            results_random_mean, results_oracle_mean))

if(!ordered){
  # Plotting
  ggplot(data_plot, aes(x = time, y = value, color = as.factor(model)) )+
    geom_line(size=1.2) +
    theme_minimal() +  # Use minimal theme for polished look
    labs(x = "Additional Samples", y = "Discoveries") +  # Set axis labels
    scale_color_brewer(palette = "Dark2") +  # Choose color palette
    theme(
      legend.position = "right",  # Position legend
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)  # Center plot title
    ) +
    ggtitle("Simulated data (unordered Zipf) - results")  # Set plot title
}else{
  # Plotting
  ggplot(data_plot, aes(x = time, y = value, color = as.factor(model)) )+
    geom_line(size=1.2) +
    theme_minimal() +  # Use minimal theme for polished look
    labs(x = "Additional Samples", y = "Discoveries") +  # Set axis labels
    scale_color_brewer(palette = "Dark2") +  # Choose color palette
    theme(
      legend.position = "right",  # Position legend
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)  # Center plot title
    ) +
    ggtitle("Simulated data (ordered Zipf) - results")  # Set plot title
}

# Average number of species discovered 
sum(diff(result_HPY_mean))     / new_samples
sum(diff(result_HDP_mean))     / new_samples
sum(diff(results_random_mean)) / new_samples
sum(diff(results_oracle_mean)) / new_samples
