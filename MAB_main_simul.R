#multiarmed bandit for species discovery via mSSP - simulation study
rm(list = ls())
library(rstudioapi) # version 0.15.0
library(ggplot2) # version 3.5.0 
library(readxl) # version 1.4.3 

#set working directory to Source file directory
#code to set the working directory to the current folder from RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("mSSPmab.R")

###############which true pmf? 
J = 8
pmfs = generate_zipf(param = c(rep(1.3, 4), rep(2, 4)), 
                     tot_species = 3000, j_species = 2500, seed = 0)

################"Plot prob of tie"
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







###############how many initial and new samples? 
init_samples = 30 #per each pop
new_samples = 300

###############how many replicas?
seed_replicas = seq(1,20)
tot_replica = length(seed_replicas)

###############initialize matrices and list to save for more replicas
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
###############gibbs samplers
replica = 0
for(seed in seed_replicas){
  
  replica = replica + 1
  
  cat("\nReplica", replica, "out of", tot_replica, "\n")
  
  ###############sample observations for fair comparison of methods
  X = sample_from_pop_all(truth = pmfs, size = init_samples + new_samples,
                          seed = seed, verbose = FALSE)
  #solve MAB decision via uniform
  results_random_temp = uniform_MAB(data = X, new_samples = new_samples, 
                                    seed = 0)
  results_random[,replica] = results_random_temp$discoveries
  
  #solve MAB decision via oracle
  results_oracle_temp = oracle_MAB(data = X, pmfs = pmfs, new_samples = new_samples)
  results_oracle[,replica] = results_oracle_temp$discoveries
  est_prob_new_oracle[[replica]] = results_oracle_temp$probs
  
  #solve MAB decisions via indepDP
  results_indepDP_temp = indepDP_MAB(data = X, new_samples = new_samples, 
                                     seed = 0)
  results_indepDP[,replica] = results_indepDP_temp$discoveries
  est_prob_new_indepDP[[replica]] = results_indepDP_temp$probs
  
  #solve MAB decisions via indepPY 
  results_indepPY_temp = indepPY_MAB(data = X, new_samples = new_samples, seed = 0)
  results_indepPY[,replica] = results_indepPY_temp$discoveries
  est_prob_new_indepPY[[replica]] = results_indepPY_temp$probs
  
  #solve MAB decisions via plusDP
  results_plusDP_temp = plusDP_MAB(data = X, new_samples = new_samples, 
                                   seed = 0)
  results_plusDP[,replica] = results_plusDP_temp$discoveries
  est_prob_new_plusDP[[replica]] = results_plusDP_temp$probs
  
  #solve MAB decisions via plusPY
  results_plusPY_temp = plusPY_MAB(data = X, new_samples = new_samples, 
                                   seed = 0)
  results_plusPY[,replica] = results_plusPY_temp$discoveries
  est_prob_new_plusPY[[replica]] = results_plusPY_temp$probs
  
  #solve MAB decisions via plusMD (not available)
  #results_plusMD = plusMD_MAB(data, new_samples = new_samples, seed = 0)
  
  #solve MAB decisions via HPY
  results_HPY_temp = HPY_MAB(data = X, new_samples = new_samples, 
                             seed = 0)
  results_HPY[,replica] = results_HPY_temp$discoveries
  est_prob_new_HPY[[replica]] = results_HPY_temp$probs   
  
  #solve MAB decisions via HDP
  results_HDP_temp = HDP_MAB(data = X, new_samples = new_samples, 
                             seed = 0)
  results_HDP[,replica] = results_HDP_temp$discoveries
  est_prob_new_HDP[[replica]] = results_HDP_temp$probs  
}

#compute average cumulative discoveries across replica
results_plusDP_mean = rowMeans( results_plusDP, na.rm = TRUE )
results_plusPY_mean = rowMeans( results_plusPY, na.rm = TRUE  )
results_indepDP_mean = rowMeans( results_indepDP, na.rm = TRUE  )
results_indepPY_mean = rowMeans( results_indepPY, na.rm = TRUE  )
results_random_mean  = rowMeans( results_random, na.rm = TRUE  )
results_oracle_mean  = rowMeans( results_oracle, na.rm = TRUE  )
results_HDP_mean = rowMeans( results_HDP, na.rm = TRUE  )
results_HPY_mean = rowMeans( results_HPY, na.rm = TRUE  )

################################################################################
#plot results 
#INDEPENDENT MODELS PLOT #######################################################
# prepare data matrix
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

#ADDITIVE MODELS PLOT #######################################################
# prepare data matrix
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

#HIERARCHICAL MODELS PLOT #######################################################
# prepare data matrix
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

#average number of species discovered 
sum(diff(results_plusDP_mean)) / new_samples
sum(diff(results_plusPY_mean)) / new_samples
sum(diff(results_indepDP_mean)) / new_samples
sum(diff(results_indepPY_mean)) / new_samples
sum(diff(results_random_mean)) / new_samples
sum(diff(results_oracle_mean)) / new_samples
sum(diff(results_HDP_mean)) / new_samples
sum(diff(results_HPY_mean)) / new_samples

#MSE

MSE_DP = 0; MSE_PY = 0;
MSE_plusDP = 0; MSE_plusPY = 0;
MSE_HDP = 0; MSE_HPY = 0;
for (replica in 1:tot_replica){
  MSE_DP = MSE_DP + sum((est_prob_new_indepDP[[replica]] - est_prob_new_oracle[[replica]])**2)
  MSE_PY = MSE_PY + sum((est_prob_new_indepPY[[replica]] - est_prob_new_oracle[[replica]])**2)
  MSE_plusDP = MSE_plusDP + sum((est_prob_new_plusDP[[replica]] - est_prob_new_oracle[[replica]])**2)
  MSE_plusPY = MSE_plusPY + sum((est_prob_new_plusPY[[replica]] - est_prob_new_oracle[[replica]])**2)
  MSE_HDP = MSE_HDP + sum((est_prob_new_HDP[[replica]] - est_prob_new_oracle[[replica]])**2)
  MSE_HPY = MSE_HPY + sum((est_prob_new_HPY[[replica]] - est_prob_new_oracle[[replica]])**2)
}
MSE_DP = MSE_DP / (new_samples*tot_replica*J); RMSE_DP = sqrt(MSE_DP)
MSE_PY = MSE_PY / (new_samples*tot_replica*J); RMSE_PY = sqrt(MSE_PY)
MSE_plusDP = MSE_plusDP / (new_samples*tot_replica*J); RMSE_plusDP = sqrt(MSE_plusDP)
MSE_plusPY = MSE_plusPY / (new_samples*tot_replica*J); RMSE_plusPY = sqrt(MSE_plusPY)
MSE_HDP = MSE_HDP / (new_samples*tot_replica*J); RMSE_HDP = sqrt(MSE_HDP)
MSE_HPY = MSE_HPY / (new_samples*tot_replica*J); RMSE_HPY = sqrt(MSE_HPY)

####################Additional Plots############################################

################"Plot Empirical Bayes prob of tie across"
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

