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
source("utils.R")
source("MAB_HPYP_fct.R")


###############which true pmf? 
J = 8
ordered = TRUE

if(!ordered){
  #generate true pmf
  pmfs = generate_zipf(param = c(rep(1.3, 2), rep(2, 6)), 
                         tot_species = 3000, j_species = 2500, seed = 0)
}else{
  pmfs = generate_zipf_reorder(param = c(rep(1.3, 2), rep(2, 6)), 
                     tot_species = 3000, j_species = 2500, seed = 0)
}


###############how many initial and new samples? 
init_samples = 30 #per each pop
new_samples = 300

############### how many replicas?
tot_replica = 10

###############initialize for more replicas
results_HPY = matrix(NA, nrow=new_samples, ncol= tot_replica)
est_prob_new_HPY = vector("list", tot_replica)


############### Gibbs samplers
for(seed in 1:tot_replica){
  cat("\nReplica", seed, "out of", tot_replica, "\n")
  
  ############### Sample observations for fair comparison of methods
  X = sample_from_pop_all(truth = pmfs, size = init_samples + new_samples,
                          seed = seed, verbose = FALSE)
  
  #solve MAB decisions via HPY
  results_HPY_temp = HPY_MAB(data = X,
                             a_alpha = 1, b_alpha = 1,
                             init_samples = init_samples, 
                             new_samples = new_samples, 
                             a_sigma = 1, b_sigma = 2,
                             burnin = 10, iters = 30, seed = 0, 
                             niter_MH = 10, ada_step = 10,
                             ada_thresh = 0.44, r_ada_input = 0)
  

  results_HPY[,seed] = results_HPY_temp$discoveries
  est_prob_new_HPY[[seed]] = results_HPY_temp$probs
}

result_avg_HPY = rowMeans(results_HPY)

################################################################################
# Plot results 

# prepare data matrix
num_model_to_compare = 1
names = c("HPY")
model = c()
for(mm in 1:num_model_to_compare){
  model = c(model, rep(names[mm], new_samples))
}
data_plot <- data.frame(
  time = rep(1:new_samples, num_model_to_compare),
  model = model,
  value = c(result_avg_HPY))

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
sum(diff(result_avg_HPY)) / new_samples
