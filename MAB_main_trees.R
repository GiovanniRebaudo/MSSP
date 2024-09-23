#multiarmed bandit for species discovery via mSSP - tree data
rm(list = ls())
library(rstudioapi) # version 0.15.0
library(ggplot2) # version 3.5.0 
library(readxl) # version 1.4.3 
library(ggpubr) # version 0.6.0  
library(viridis) # version 0.6.5

#set working directory to Source file directory
#code to set the working directory to the current folder from RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("mSSPmab.R")

#read data
data_trees = read_excel("conditwebtable.xls")

#associate a numeric unique label to each species
data_trees = cbind(seq(1:nrow(data_trees)), data_trees)
colnames(data_trees)[1] = "cod_species"

#aggregate plots
data_trees_4pop = data.frame(species = as.factor(data_trees[,1]),
                             Pop1 = rowSums(data_trees[, grep("^B", colnames(data_trees))]),
                             Pop2 = rowSums(data_trees[, grep("^P", colnames(data_trees))]),
                             Pop3 = rowSums(data_trees[, grep("^S", colnames(data_trees))]),
                             Pop4 = rowSums(data_trees[, grep("^C", colnames(data_trees))]), 
                             Pop1_rel = rowSums(data_trees[, grep("^B", colnames(data_trees))])/
                               sum(rowSums(data_trees[, grep("^B", colnames(data_trees))])),
                             Pop2_rel = rowSums(data_trees[, grep("^P", colnames(data_trees))])/
                               sum(rowSums(data_trees[, grep("^P", colnames(data_trees))])), 
                             Pop3_rel = rowSums(data_trees[, grep("^S", colnames(data_trees))])/
                               sum(rowSums(data_trees[, grep("^S", colnames(data_trees))])),
                             Pop4_rel = rowSums(data_trees[, grep("^C", colnames(data_trees))])/
                               sum(rowSums(data_trees[, grep("^C", colnames(data_trees))])) )

#reorder species based on their frequency in the second Pop for plotting them 
order_labels = order(data_trees_4pop$Pop2, decreasing = TRUE)
data_trees_4pop_reorder = data.frame(species = as.factor(data_trees[,1]),
                                     Pop1_rel = data_trees_4pop$Pop1_rel[order_labels],
                                     Pop2_rel = data_trees_4pop$Pop2_rel[order_labels],
                                     Pop3_rel = data_trees_4pop$Pop3_rel[order_labels],
                                     Pop4_rel = data_trees_4pop$Pop4_rel[order_labels] )
#plot empirical distributions
p1 = ggbarplot(data_trees_4pop_reorder, x = "species", y = "Pop1_rel",
          color = viridis(4)[1],  
          x.text.angle = 90 )
p2 = ggbarplot(data_trees_4pop_reorder, x = "species", y = "Pop2_rel",
               color = viridis(4)[2], 
               x.text.angle = 90 )
p3 = ggbarplot(data_trees_4pop_reorder, x = "species", y = "Pop3_rel",
               color = viridis(4)[3], 
               x.text.angle = 90 )
p4 = ggbarplot(data_trees_4pop_reorder, x = "species", y = "Pop4_rel",
               color = viridis(4)[4], 
               x.text.angle = 90 )

ggarrange(p1 + rremove("x.text") + rremove("ylab"), 
          p2 + rremove("x.text") + rremove("ylab"), 
          p3 + rremove("x.text") + rremove("ylab"), 
          p4 + rremove("x.text") + rremove("ylab"),  
          labels = c("   Group 1", "   Group 2",
                     "   Group 3", "   Group 4"),
          ncol = 2, nrow = 2)

sample1 = rep(1:nrow(data_trees_4pop), data_trees_4pop$Pop1)
sample2 = rep(1:nrow(data_trees_4pop), data_trees_4pop$Pop2)
sample3 = rep(1:nrow(data_trees_4pop), data_trees_4pop$Pop3)
sample4 = rep(1:nrow(data_trees_4pop), data_trees_4pop$Pop4)

J = 4
pmfs = vector("list", J) 
pmfs[[1]] = data_trees_4pop$Pop1_rel
pmfs[[2]] = data_trees_4pop$Pop2_rel
pmfs[[3]] = data_trees_4pop$Pop3_rel
pmfs[[4]] = data_trees_4pop$Pop4_rel
################"Plot empirical prob of tie"
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
  geom_text(aes(label = format(ptie, scientific = FALSE, digits = 1) ),
            color = "black")+
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank()) + 
  guides(fill=guide_legend(title="Prob. tie")) 

J = 4
init_samples = 30 #per each pop
new_samples = 300

###############how many replicas?
tot_replica = 20

###############initialize for more replicas
results_plusDP_real = matrix(NA, nrow = new_samples, ncol = tot_replica) 
results_plusPY_real = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_indepDP_real = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_indepPY_real = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_HDP_real = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_HPY_real = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_random_real  = matrix(NA, nrow = new_samples, ncol = tot_replica)
results_oracle_real  = matrix(NA, nrow = new_samples, ncol = tot_replica)

est_prob_new_plusDP_real = vector("list", tot_replica)
est_prob_new_plusPY_real = vector("list", tot_replica)
est_prob_new_HDP_real = vector("list", tot_replica)
est_prob_new_HPY_real = vector("list", tot_replica)
est_prob_new_indepDP_real = vector("list", tot_replica)
est_prob_new_indepPY_real = vector("list", tot_replica)
est_prob_new_oracle_real = vector("list", tot_replica)

###############gibbs samplers
replica = 0
for(seed in 1:tot_replica){
  
  replica = replica + 1
  
  cat("\nReplica", replica, "out of", tot_replica, "\n")
  
  ###############sample observations for fair comparison of methods
  set.seed(seed)
  X = matrix(NA, nrow = J, ncol = init_samples+new_samples)
  X[1,] = sample(sample1, init_samples+new_samples, replace = FALSE)
  X[2,] = sample(sample2, init_samples+new_samples, replace = FALSE)
  X[3,] = sample(sample3, init_samples+new_samples, replace = FALSE)
  X[4,] = sample(sample4, init_samples+new_samples, replace = FALSE)
  
  #solve MAB decisions via plusDP
  results_plusDP_temp = plusDP_MAB(data = X, new_samples = new_samples, 
                                   seed = 0)
  results_plusDP_real[,replica] = results_plusDP_temp$discoveries
  est_prob_new_plusDP_real[[replica]] = results_plusDP_temp$probs
  
  #solve MAB decisions via plusPY
  results_plusPY_temp = plusPY_MAB(data = X, new_samples = new_samples, 
                                   seed = 0)
  results_plusPY_real[,replica] = results_plusPY_temp$discoveries
  est_prob_new_plusPY_real[[replica]] = results_plusPY_temp$probs
  
  #solve MAB decisions via indepDP
  results_indepDP_temp = indepDP_MAB(data = X, new_samples = new_samples, 
                                     seed = 0)
  results_indepDP_real[,replica] = results_indepDP_temp$discoveries
  est_prob_new_indepDP_real[[replica]] = results_indepDP_temp$probs
  
  #solve MAB decisions via indepPY 
  results_indepPY_temp = indepPY_MAB(data = X, new_samples = new_samples, 
                                     seed = 0)
  results_indepPY_real[,replica] = results_indepPY_temp$discoveries
  est_prob_new_indepPY_real[[replica]] = results_indepPY_temp$probs
  
  #solve MAB decisions via plusMD (not available)
  #results_plusMD = plusMD_MAB(data, new_samples = new_samples, seed = 0)
  
  #solve MAB decision via uniform
  results_random_temp = uniform_MAB(data = X, new_samples = new_samples, 
                                    seed = 0)
  results_random_real[,replica] = results_random_temp$discoveries
  
  #solve MAB decisions via HPY
  results_HPY_temp = HPY_MAB(data = X, new_samples = new_samples, 
                             seed = 0)
  results_HPY_real[,replica] = results_HPY_temp$discoveries
  est_prob_new_HPY_real[[replica]] = results_HPY_temp$probs   #t(prob)
  
  #solve MAB decisions via HDP
  results_HDP_temp = HDP_MAB(data = X, new_samples = new_samples, 
                             seed = 0)
  results_HDP_real[,replica] = results_HDP_temp$discoveries
  est_prob_new_HDP_real[[replica]] = results_HDP_temp$probs  #t(prob)
}

results_plusDP_mean = rowMeans( results_plusDP_real, na.rm = TRUE )
results_plusPY_mean = rowMeans( results_plusPY_real, na.rm = TRUE )
results_HDP_mean = rowMeans( results_HDP_real, na.rm = TRUE )
results_HPY_mean = rowMeans( results_HPY_real, na.rm = TRUE )
results_indepDP_mean = rowMeans( results_indepDP_real, na.rm = TRUE )
results_indepPY_mean = rowMeans( results_indepPY_real, na.rm = TRUE )
results_random_mean  = rowMeans( results_random_real, na.rm = TRUE )


################################################################################
#plot results 
#INDEPENDENT MODELS PLOT #######################################################
# prepare data matrix

names = c("Uniform", "Ind DP", "Ind PY")
num_model_to_compare = length(names)
model = c()
for(mm in 1:num_model_to_compare){
  model = c(model, rep(names[mm], new_samples))
}
data_plot <- data.frame(
  time = rep(1:new_samples, num_model_to_compare),
  model = model,
  value = c(results_random_mean, results_indepDP_mean, results_indepPY_mean))

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
  )  + scale_y_continuous(limits = c(0, 100)) +
  ggtitle("Independent Processes")  # Set plot title

#ADDITIVE MODELS PLOT #######################################################
# prepare data matrix
names = c("Uniform", "+DP", "+PY")
num_model_to_compare = length(names)
model = c()
for(mm in 1:num_model_to_compare){
  model = c(model, rep(names[mm], new_samples))
}
data_plot <- data.frame(
  time = rep(1:new_samples, num_model_to_compare),
  model = model,
  value = c(results_random_mean,
            results_plusDP_mean, results_plusPY_mean))

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
  )  + scale_y_continuous(limits = c(0, 100)) +
  ggtitle("Additive processes")  # Set plot title

#HIERARCHICAL MODELS PLOT #######################################################
# prepare data matrix
names = c("Uniform", "HDP", "HPY")
num_model_to_compare = length(names)
model = c()
for(mm in 1:num_model_to_compare){
  model = c(model, rep(names[mm], new_samples))
}
data_plot <- data.frame(
  time = rep(1:new_samples, num_model_to_compare),
  model = model,
  value = c(results_random_mean, results_HDP_mean, results_HPY_mean))

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
  )  + scale_y_continuous(limits = c(0, 100)) +
  ggtitle("Hierarchical Processes")  # Set plot title

#average number of species discovered 
sum(diff(results_plusDP_mean)) / new_samples
sum(diff(results_plusPY_mean)) / new_samples
sum(diff(results_indepDP_mean)) / new_samples
sum(diff(results_random_mean)) / new_samples
sum(diff(results_indepPY_mean)) / new_samples
sum(diff(results_HDP_mean)) / new_samples
sum(diff(results_HPY_mean)) / new_samples
