library(ggplot2)
#A priori plots 
J = 2
corr_vec = c(0, 0.25, 0.5, 0.75, 1)#when alpha param equal, it equals eps
n1 = 10000
data = matrix(NA, nrow = length(corr_vec)*n1, ncol = 2)
tot_sim = 1000

#HDP ###########################################################################
ptiewith = 0.5

count = 0

for(corr in corr_vec){
  print(corr)
  #label
  data[(count+1):(count+n1), 1] = rep(corr, n1)
  
  if(corr == 0){
    
    p_new = rep(1, n1)
    
  }else if (corr == 1){
    
    alpha0 = 1
    
    p_new = alpha0 / (alpha0 + seq(1, n1) )
    
  }else{
  
    alpha = (1 - ptiewith) / ptiewith / (1 - corr)
    alpha0 = (1 + alpha) * (1 - corr) / corr
    p_new =  matrix(NA, nrow = n1, ncol = tot_sim)
    # Initializes the progress bar
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = tot_sim, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
    for(sim in 1:tot_sim){
      r = 1 
      p_new_temp = alpha0 / (alpha0 + 1)
      for(i in 2:n1){
        temp_r = sample(c(r,max(r)+1), 1, prob = c(rep(1, i-1), alpha))
        r = c(r, temp_r)
        p_new_temp = c(p_new_temp, alpha0 / (alpha0 + length(unique(r))))
      }
      p_new[,sim] = p_new_temp
      setTxtProgressBar(pb, sim)
    }
    p_new = apply(p_new, 1, mean)
    
  }
  
  data[(count+1):(count+n1), 2] = p_new
  
  count = count + n1
}

data_plot_HDP <- data.frame(
  samplesize = rep( seq(1,n1), length(corr_vec)),
  correlation = as.factor(data[, 1]),
  p_new = data[, 2])

ggplot(data_plot_HDP, aes(x = samplesize, y = p_new, color = correlation) )+
  geom_line(aes(linetype=correlation), size=1.2) +
  theme_minimal() +  # Use minimal theme for polished look
  labs(x = "# of subjects observed from pop. k", y = "prob of new species in pop. j") +  # Set axis labels
  scale_color_brewer(palette = "Dark2") +  # Choose color palette
  theme(text = element_text(size = 20),
        legend.position = "right",  # Position legend
        plot.title = element_text(hjust = 0.5)  # Center plot title
  ) + scale_x_continuous(trans='log10') +
  ggtitle("HDP \n with prob. of tie within = 0.5") 

#NDP ###########################################################################
ptiewith = 0.5
beta = ptiewith / (1 - ptiewith)

data = matrix(NA, nrow = length(corr_vec)*n1, ncol = 2)
count = 0

for(corr in corr_vec){
  print(corr)
  #label
  data[(count+1):(count+n1), 1] = rep(corr, n1)
  
  if(corr == 0){
    
    data[(count+1):(count+n1), 2] = rep(1, n1)
    
  }else{
    
    alpha = (1 - corr) / corr
    
    #prob
    data[(count+1):(count+n1), 2] = beta / (beta + seq(1,n1) ) / (alpha + 1) + 
      alpha / (alpha + 1)
  }
  
  
  #counter
  count = count + n1
}

data_plot_NDP <- data.frame(
  samplesize = rep( seq(1,n1), length(corr_vec)),
  correlation = as.factor(data[, 1]),
  p_new = data[, 2])

ggplot(data_plot_NDP, aes(x = samplesize, y = p_new, color = correlation) )+
  geom_line(aes(linetype=correlation), size=1.2) +
  theme_minimal() +  # Use minimal theme for polished look
  labs(x = "# of subjects observed from pop. k", y = "prob of new species in pop. j") +  # Set axis labels
  scale_color_brewer(palette = "Dark2") +  # Choose color palette
  theme(text = element_text(size = 20),
    legend.position = "right",  # Position legend
    plot.title = element_text(hjust = 0.5)  # Center plot title
  ) + scale_x_continuous(trans='log10') +
  ggtitle("NDP \n with prob. of tie within = 0.5") 




#+DP ###########################################################################
ptiewith = 0.5

eps_from_rho <- function(rho) {
  sqrt(rho) / (sqrt(rho) + sqrt(1 - rho))
}

alpha_from_eps <- function(eps, ptiewith = 0.5) {
  (eps^2 + (1 - eps)^2) / ptiewith - 1
}

count = 0

for (rho in corr_vec) {
  
  print(rho)
  
  eps   = eps_from_rho(rho)
  alpha = alpha_from_eps(eps, ptiewith)
  
  # Indicator that an observation from population k comes from Q_0
  common_indicator <- matrix(
    rbinom(n1 * tot_sim, size = 1, prob = eps),
    nrow = n1,
    ncol = tot_sim
  )
  
  n_common = apply(common_indicator, 2, cumsum)
  
  # Probability that a draw from Q_0 is new
  if (alpha < 1e-12) {
    
    # Limit as alpha -> 0
    p_new_Q0 = (n_common == 0)
    
  } else {
    
    p_new_Q0 = alpha / (alpha + n_common)
    
  }
  
  # A draw from Q_j is always new relative to population k
  p_new = rowMeans(
    (1 - eps) + eps * p_new_Q0
  )
  
  index = count + seq_len(n1)
  
  data[index, 1] = rho
  data[index, 2] = p_new
  
  count = count + n1
}

data_plot_plusDP <- data.frame(
  samplesize = rep( seq(1,n1), length(corr_vec)),
  correlation = as.factor(data[, 1]),
  p_new = data[, 2])

ggplot(data_plot_plusDP, aes(x = samplesize, y = p_new, color = correlation) )+
  geom_line(aes(linetype=correlation), size=1.2) +
  theme_minimal() +  # Use minimal theme for polished look
  labs(x = "# of subjects observed from pop. k", y = "prob of new species in pop. j") +  # Set axis labels
  scale_color_brewer(palette = "Dark2") +  # Choose color palette
  theme(text = element_text(size = 20),
        legend.position = "right",  # Position legend
        plot.title = element_text(hjust = 0.5)  # Center plot title
  ) + scale_x_continuous(trans='log10') +
  ggtitle("+DP \n with prob. of tie within = 0.5") 
