#multiarmed bandit for species discovery via mSSP
rm(list = ls())
library(rstudioapi) # version 0.15.0
library(ggplot2) # version 3.5.0 

#set working directory to Source file directory
#code to set the working directory to the current folder from RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("mSSPmab.R")

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

################"Empirical Bayes prob of tie across"
#ptie matrix
ptie = matrix(NA, nrow = J, ncol = J)
row = matrix(rep(1:J,J), nrow = J)
col = t(row)

for(j in 1:J){
  for(jj in j:J){
    ptie[j,jj] = sum( ( pmfs[[j]] / sum( pmfs[[j]] ) ) * 
                        ( pmfs[[jj]] / sum( pmfs[[jj]] ) ) )
  }
}

temp = t(ptie)
ptie[row>col] = temp[row>col]

#plot prob tie matrix
x = paste0("Pop", seq(1,J))
y = paste0("Pop", seq(1,J))
data = expand.grid(X=x, Y=y)
data$ptie = as.vector(ptie)

# Heatmap 
ggplot(data, aes(X, Y, fill= ptie)) + 
  geom_tile()+
  geom_text(aes(label = format(ptie, scientific = TRUE, digits = 3) ),
            color = "white")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank()) + 
  guides(fill=guide_legend(title="Prob. tie"))

###############how many initial and new samples? 
init_samples = 30 #per each pop
new_samples = 300

###############how many replicas?
tot_replica = 10

###############initialize for more replicas
results_plusDP = 0 
results_plusPY = 0
results_indepDP = 0
results_indepPY = 0
results_random  = 0
results_oracle  = matrix(NA, nrow = new_samples, ncol = tot_replica)
est_prob_new_plusDP = vector("list", tot_replica)
est_prob_new_plusPY = vector("list", tot_replica)
est_prob_new_indepDP = vector("list", tot_replica)
est_prob_new_indepPY = vector("list", tot_replica)
est_prob_new_oracle = vector("list", tot_replica)

###############gibbs samplers
replica = 0
for(seed in 1:tot_replica){
  
  replica = replica + 1
  
  cat("\nReplica", replica, "out of", tot_replica, "\n")
  
  ###############sample observations for fair comparison of methods
  X = sample_from_pop_all(truth = pmfs, size = init_samples + new_samples,
                          seed = seed, verbose = FALSE)
  
  #solve MAB decisions via plusDP
  results_plusDP_temp = plusDP_MAB(data = X, new_samples = new_samples, 
                                   seed = 0)
  results_plusDP[,seed] = results_plusDP_temp$discoveries
  est_prob_new_plusDP[[seed]] = results_plusDP_temp$probs
  
  #solve MAB decisions via plusPY
  results_plusPY_temp = plusPY_MAB(data = X, new_samples = new_samples, 
                                   seed = 0)
  results_plusPY[,seed] = results_plusPY_temp$discoveries
  est_prob_new_plusPY[[seed]] = results_plusPY_temp$probs
  
  #solve MAB decisions via indepDP
  results_indepDP_temp = indepDP_MAB(data = X, new_samples = new_samples, 
                                     seed = 0)
  results_indepDP[,seed] = results_indepDP_temp$discoveries
  est_prob_new_indepDP[[seed]] = results_indepDP_temp$probs
  
  #solve MAB decisions via indepPY (not available)
  results_indepPY_temp = indepPY_MAB(data = X, new_samples = new_samples, seed = 0)
  results_indepPY[,seed] = results_indepPY_temp$discoveries
  est_prob_new_indepPY[[seed]] = results_indepPY_temp$probs
  
  #solve MAB decisions via plusMD (not available)
  #results_plusMD = plusMD_MAB(data, new_samples = new_samples, seed = 0)
  
  #solve MAB decision via uniform
  results_random_temp = uniform_MAB(data = X, new_samples = new_samples, 
                                    seed = 0)
  results_random = results_random + results_random_temp$discoveries

  results_oracle_temp = oracle_MAB(data = X, pmfs = pmfs)
  results_oracle[,seed] = results_oracle_temp$discoveries
  est_prob_new_oracle[[seed]] = results_oracle_temp$probs
}

results_plusDP = rowMeans( results_plusDP )
results_plusPY = rowMeans( results_plusPY )
results_indepDP = rowMeans( results_indepDP )
results_indepPY = rowMeans( results_indepPY )
results_random  = rowMeans( results_random )
results_oracle  = rowMeans( results_oracle )


################################################################################
#plot results 

# prepare data matrix
num_model_to_compare = 6
names = c("Uniform", "Ind DP", "Ind PY", "+DP", "+PY", "Oracle")
model = c()
for(mm in 1:num_model_to_compare){
  model = c(model, rep(names[mm], new_samples))
}
data_plot <- data.frame(
  time = rep(1:new_samples, num_model_to_compare),
  model = model,
  value = c(results_random, results_indepDP, results_indepPY,
            results_plusDP, results_plusPY,
            results_oracle))

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

#average number of species discovered 
sum(diff(results_plusDP)) / new_samples
sum(diff(results_plusPY)) / new_samples
sum(diff(results_indepDP)) / new_samples
sum(diff(results_random)) / new_samples
sum(diff(results_oracle)) / new_samples
