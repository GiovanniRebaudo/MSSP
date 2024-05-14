# Prior elicitation
# HDP
alpha0Samples = rgamma(2000, shape=1, rate=1/3)
alpha1Samples = rgamma(2000, shape=1, rate=1/2)

probTiesWithinSamples = 
  (1+alpha0Samples+ alpha1Samples)/((1+alpha0Samples)*(1+alpha1Samples))

hist(probTiesWithinSamples)
mean(probTiesWithinSamples);var(probTiesWithinSamples)

# DP
alphaSamples = rgamma(2000, shape=1, rate=1)
probTiesWithinSamples = 1/(alphaSamples+1)
hist(probTiesWithinSamples)
mean(probTiesWithinSamples);var(probTiesWithinSamples)

# HPYP

