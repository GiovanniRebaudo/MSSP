# Multiarmed bandit for species discovery via mSSP - simulation study
rm(list = ls())
library(rstudioapi) # version 0.15.0
library(ggplot2) # version 3.5.0 
library(readxl) # version 1.4.3 

#set working directory to Source file directory
#code to set the working directory to the current folder from RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("mSSPmab.R")

############### true pmf
J = 8
pmfs = generate_zipf(param = c(rep(1.3, 4), rep(2, 4)), 
                     tot_species = 3000, j_species = 2500, seed = 0)

################"Plot prob of tie"
# Compute ptie matrix
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

# Plot prob tie matrix
x = paste0("Group", seq(1,J))
y = paste0("Group", seq(1,J))
data = expand.grid(X=x, Y=y)
data$ptie = as.vector(ptie)

# Heatmap 
ggplot(data, aes(X, Y, fill= ptie)) + 
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1)+
  scale_fill_gradient2(low = "#FFFFCC",
                       high = "#075AFF") +
  geom_text(aes(label = format(ptie, scientific = TRUE, digits = 1) ),
            color = "black")+
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank()) + 
  guides(fill=guide_legend(title="Prob. tie")) 







############### How many initial and new samples? 
init_samples = 30 # in each pop
new_samples = 300

############### How many replicas?
seed_replicas = seq(1,20)
tot_replica = length(seed_replicas)
n_workers = 4L # Set to 1L for the same calculation run sequentially.

############### Initialize matrices and list to save for more replicas
results_plusDP = matrix(NA, nrow = new_samples, ncol = tot_replica) 
results_plusPY = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_indepDP = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_indepPY = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_random  = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_oracle  = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_HPY  = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_HDP  = matrix(NA, nrow = new_samples, ncol = tot_replica)

est_prob_new_plusDP = vector("list", tot_replica)
est_prob_new_plusPY = vector("list", tot_replica)
est_prob_new_indepDP = vector("list", tot_replica)
est_prob_new_indepPY = vector("list", tot_replica)
est_prob_new_oracle = vector("list", tot_replica)
est_prob_new_HDP = vector("list", tot_replica)
est_prob_new_HPY = vector("list", tot_replica)
true_prob_new = vector("list", tot_replica)
############### Gibbs samplers
run_simulation_replica = function(seed, pmfs, init_samples, new_samples){
  
  ############### Sample observations for fair comparison of methods
  X = sample_from_pop_all(truth = pmfs, size = init_samples + new_samples,
                          seed = seed, verbose = FALSE)
  # Solve MAB decision via uniform
  results_random_temp = uniform_MAB(data = X, new_samples = new_samples, 
                                    init_samples = init_samples, seed = 0)
  
  # Solve MAB decision via oracle
  results_oracle_temp = oracle_MAB(data = X, pmfs = pmfs, new_samples = new_samples,
                                  init_samples = init_samples)
  
  # Solve MAB decisions via indepDP
  results_indepDP_temp = indepDP_MAB(data = X, new_samples = new_samples, 
                                     init_samples = init_samples, seed = 0)
  
  # Solve MAB decisions via indepPY 
  results_indepPY_temp = indepPY_MAB(data = X, new_samples = new_samples,
                                     init_samples = init_samples, seed = 0)
  
  # Solve MAB decisions via plusDP
  results_plusDP_temp = plusDP_MAB(data = X, new_samples = new_samples, 
                                   init_samples = init_samples, seed = 0)
  
  # Solve MAB decisions via plusPY
  results_plusPY_temp = plusPY_MAB(data = X, new_samples = new_samples, 
                                   init_samples = init_samples, seed = 0)
  
  # Solve MAB decisions via HPY
  results_HPY_temp = HPY_MAB(data = X, new_samples = new_samples, 
                             init_samples = init_samples, seed = 0)
  
  # Solve MAB decisions via HDP
  results_HDP_temp = HDP_MAB(data = X, new_samples = new_samples, 
                             init_samples = init_samples, seed = 0)

  # True discovery probabilities along each strategy's OWN sampling history.
  selected_arms = list(indepDP = results_indepDP_temp$selected_arms,
                       indepPY = results_indepPY_temp$selected_arms,
                       plusDP = results_plusDP_temp$selected_arms,
                       plusPY = results_plusPY_temp$selected_arms,
                       HDP = results_HDP_temp$selected_arms,
                       HPY = results_HPY_temp$selected_arms)
  true_prob_new = lapply(selected_arms, function(arms)
    true_discovery_probs(X, pmfs, init_samples, arms))
  list(random = results_random_temp, oracle = results_oracle_temp,
       indepDP = results_indepDP_temp, indepPY = results_indepPY_temp,
       plusDP = results_plusDP_temp, plusPY = results_plusPY_temp,
       HPY = results_HPY_temp, HDP = results_HDP_temp,
       true_prob_new = true_prob_new)
}

replica_results = run_mab_replicas(seed_replicas, run_simulation_replica,
                                  pmfs = pmfs, init_samples = init_samples,
                                  new_samples = new_samples, workers = n_workers,
                                  progress = TRUE)
# Collect in the original replica order; leave all subsequent summaries unchanged.
for(replica in seq_len(tot_replica)){
  result = replica_results[[replica]]
  results_random[,replica] = result$random$discoveries
  results_oracle[,replica] = result$oracle$discoveries
  est_prob_new_oracle[[replica]] = result$oracle$probs
  results_indepDP[,replica] = result$indepDP$discoveries
  est_prob_new_indepDP[[replica]] = result$indepDP$probs
  results_indepPY[,replica] = result$indepPY$discoveries
  est_prob_new_indepPY[[replica]] = result$indepPY$probs
  results_plusDP[,replica] = result$plusDP$discoveries
  est_prob_new_plusDP[[replica]] = result$plusDP$probs
  results_plusPY[,replica] = result$plusPY$discoveries
  est_prob_new_plusPY[[replica]] = result$plusPY$probs
  results_HPY[,replica] = result$HPY$discoveries
  est_prob_new_HPY[[replica]] = result$HPY$probs
  results_HDP[,replica] = result$HDP$discoveries
  est_prob_new_HDP[[replica]] = result$HDP$probs
  true_prob_new[[replica]] = result$true_prob_new
}

# Compute average cumulative discoveries across replica
results_plusDP_mean = rowMeans( results_plusDP, na.rm = TRUE )
results_plusPY_mean = rowMeans( results_plusPY, na.rm = TRUE  )
results_indepDP_mean = rowMeans( results_indepDP, na.rm = TRUE  )
results_indepPY_mean = rowMeans( results_indepPY, na.rm = TRUE  )
results_random_mean  = rowMeans( results_random, na.rm = TRUE  )
results_oracle_mean  = rowMeans( results_oracle, na.rm = TRUE  )
results_HDP_mean = rowMeans( results_HDP, na.rm = TRUE  )
results_HPY_mean = rowMeans( results_HPY, na.rm = TRUE  )

################################################################################
# Plot results 
# INDEPENDENT MODELS PLOT #######################################################
# Prepare data matrix
names = c("Uniform", "Ind DP", "Ind PY", "Oracle")
num_model_to_compare = length(names)
model = c()
for(mm in 1:num_model_to_compare){
  model = c(model, rep(names[mm], new_samples))
}
data_plot <- data.frame(
  time = rep(1:new_samples, num_model_to_compare),
  model = model,
  value = c(results_random_mean, results_indepDP_mean, results_indepPY_mean,
            results_oracle_mean))

# Plotting
ggplot(data_plot, aes(x = time, y = value, color = as.factor(model)) )+
  geom_line(aes(linetype = as.factor(model)), size=1.2) +
  theme_minimal() +  # Use minimal theme for polished look
  labs(x = "Additional Samples", y = "Discoveries") +  # Set axis labels
  scale_color_brewer(palette = "Dark2") +  # Choose color palette
  theme(text = element_text(size = 20),
    legend.position = "right",  # Position legend
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # Center plot title
  ) +
  ggtitle("Independent Processes")  # Set plot title

# ADDITIVE MODELS PLOT #######################################################
# Prepare data matrix
names = c("Uniform", "+DP", "+PY", "Oracle")
num_model_to_compare = length(names)
model = c()
for(mm in 1:num_model_to_compare){
  model = c(model, rep(names[mm], new_samples))
}
data_plot <- data.frame(
  time = rep(1:new_samples, num_model_to_compare),
  model = model,
  value = c(results_random_mean,
            results_plusDP_mean, results_plusPY_mean,
            results_oracle_mean))

# Plotting
ggplot(data_plot, aes(x = time, y = value, color = as.factor(model)) )+
  geom_line(aes(linetype = as.factor(model)), size=1.2) +
  theme_minimal() +  # Use minimal theme for polished look
  labs(x = "Additional Samples", y = "Discoveries") +  # Set axis labels
  scale_color_brewer(palette = "Dark2") +  # Choose color palette
  theme(text = element_text(size = 20),
    legend.position = "right",  # Position legend
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # Center plot title
  ) +
  ggtitle("Additive processes")  # Set plot title

# HIERARCHICAL MODELS PLOT #######################################################
# Prepare data matrix
names = c("Uniform", "HDP", "HPY", "Oracle")
num_model_to_compare = length(names)
model = c()
for(mm in 1:num_model_to_compare){
  model = c(model, rep(names[mm], new_samples))
}
data_plot <- data.frame(
  time = rep(1:new_samples, num_model_to_compare),
  model = model,
  value = c(results_random_mean, results_HDP_mean, results_HPY_mean,
            results_oracle_mean))

# Plotting
ggplot(data_plot, aes(x = time, y = value, color = as.factor(model)) )+
  geom_line(aes(linetype = as.factor(model)), size=1.2) +
  theme_minimal() +  # Use minimal theme for polished look
  labs(x = "Additional Samples", y = "Discoveries") +  # Set axis labels
  scale_color_brewer(palette = "Dark2") +  # Choose color palette
  theme(text = element_text(size = 20),
        legend.position = "right",  # Position legend
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)  # Center plot title
  ) +
  ggtitle("Hierarchical Processes")  # Set plot title

# Average number of species discovered 
mean(results_plusDP[nrow(results_plusDP), ] / new_samples)
mean(results_plusPY[nrow(results_plusPY), ] / new_samples)
mean(results_indepDP[nrow(results_indepDP), ] / new_samples)
mean(results_indepPY[nrow(results_indepPY), ] / new_samples)
mean(results_random[nrow(results_random), ] / new_samples)
mean(results_oracle[nrow(results_oracle), ] / new_samples)
mean(results_HDP[nrow(results_HDP), ] / new_samples)
mean(results_HPY[nrow(results_HPY), ] / new_samples)

#MSE relative to the true unseen mass on each strategy's own trajectory

MSE_DP = 0; MSE_PY = 0;
MSE_plusDP = 0; MSE_plusPY = 0;
MSE_HDP = 0; MSE_HPY = 0;
for (replica in 1:tot_replica){
  MSE_DP = MSE_DP + sum((est_prob_new_indepDP[[replica]] - true_prob_new[[replica]]$indepDP)**2)
  MSE_PY = MSE_PY + sum((est_prob_new_indepPY[[replica]] - true_prob_new[[replica]]$indepPY)**2)
  MSE_plusDP = MSE_plusDP + sum((est_prob_new_plusDP[[replica]] - true_prob_new[[replica]]$plusDP)**2)
  MSE_plusPY = MSE_plusPY + sum((est_prob_new_plusPY[[replica]] - true_prob_new[[replica]]$plusPY)**2)
  MSE_HDP = MSE_HDP + sum((est_prob_new_HDP[[replica]] - true_prob_new[[replica]]$HDP)**2)
  MSE_HPY = MSE_HPY + sum((est_prob_new_HPY[[replica]] - true_prob_new[[replica]]$HPY)**2)
}
MSE_DP = MSE_DP / (new_samples*tot_replica*J); RMSE_DP = sqrt(MSE_DP)
MSE_PY = MSE_PY / (new_samples*tot_replica*J); RMSE_PY = sqrt(MSE_PY)
MSE_plusDP = MSE_plusDP / (new_samples*tot_replica*J); RMSE_plusDP = sqrt(MSE_plusDP)
MSE_plusPY = MSE_plusPY / (new_samples*tot_replica*J); RMSE_plusPY = sqrt(MSE_plusPY)
MSE_HDP = MSE_HDP / (new_samples*tot_replica*J); RMSE_HDP = sqrt(MSE_HDP)
MSE_HPY = MSE_HPY / (new_samples*tot_replica*J); RMSE_HPY = sqrt(MSE_HPY)

####################Additional Plots############################################

################"Plot Empirical Bayes prob of tie across"
# ptie matrix
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

# Plot prob tie matrix
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

# Temp
c(
  Ind_DP = RMSE_DP,
  Ind_PY = RMSE_PY,
  Add_DP = RMSE_plusDP,
  Add_PY = RMSE_plusPY,
  HDP = RMSE_HDP,
  HPY = RMSE_HPY
)
save.image(file = "Data-and-Results/MAB_simul.RData")
