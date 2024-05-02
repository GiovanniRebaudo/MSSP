#multiarmed bandit for species discovery via mSSP
rm(list = ls())
library(rstudioapi) # version 0.15.0
library(ggplot2) # version 3.5.0 

#set working directory to Source file directory
#code to set the working directory to the current folder from RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("mSSPmab.R")

#which truth? 
ordered = FALSE

#how many new sample? 
new_samples = 300

if(!ordered){
  #generate true pmf
  data = generate_zipf(param = c(rep(1.3, 2), rep(2, 6)), 
                         tot_species = 3000, j_species = 2500, seed = 0)
}else{
  data = generate_zipf_reorder(param = c(rep(1.3, 2), rep(2, 6)), 
                     tot_species = 3000, j_species = 2500, seed = 0)
}


#solve MAB decisions via plusDP
results_plusDP = plusDP_MAB(data, new_samples = new_samples, seed = 0)

#solve MAB decisions via plusPY
results_plusPY = plusPY_MAB(data, new_samples = new_samples, seed = 0)

results_indepDP = indepDP_MAB(data, new_samples = new_samples, seed = 0)

#solve MAB decisions via plusMD (not available)
#results_plusMD = plusMD_MAB(data, new_samples = new_samples, seed = 0)

#solve MAB decision via uniform
results_random = uniform_MAB(data, new_samples = new_samples, seed = 0)

#solve MAB decision via oracle
results_oracle = oracle_MAB(data, seed = 0)

#plot results 

# prepare data matrix
num_model_to_compare = 5
names = c("+DP", "+PY", "Ind DP", "Uniform", "Oracle")
model = c()
for(mm in 1:num_model_to_compare){
  model = c(model, rep(names[mm], new_samples))
}
data_plot <- data.frame(
  time = rep(1:new_samples, num_model_to_compare),
  model = model,
  value = c(results_plusDP, results_plusPY, results_indepDP, results_random, 
            results_oracle = Oracle_MAP(data)))

if(!ordered){
  # Plotting
  ggplot(data_plot, aes(x = time, y = value, color = as.factor(model)) )+
    geom_line() +
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
    geom_line() +
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

#average number of species discovered 
sum(diff(results_plusDP)) / new_samples
sum(diff(results_plusPY)) / new_samples
sum(diff(results_indepDP)) / new_samples
sum(diff(results_random)) / new_samples
sum(diff(results_oracle)) / new_samples
results_oracle