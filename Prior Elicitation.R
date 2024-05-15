# Prior elicitation
nSamples = 2000
# HDP
alpha0Samples = rgamma(nSamples, shape=1, rate=1/3)
alpha1Samples = rgamma(nSamples, shape=1, rate=1/2)

probTiesWithinSamples = 
  (1+alpha0Samples+ alpha1Samples)/((1+alpha0Samples)*(1+alpha1Samples))

hist(probTiesWithinSamples)
mean(probTiesWithinSamples);var(probTiesWithinSamples)

# DP
alphaSamples = rgamma(nSamples, shape=1, rate=1)
probTiesWithinSamples = 1/(alphaSamples+1)
hist(probTiesWithinSamples)
mean(probTiesWithinSamples);var(probTiesWithinSamples)

# HPYP

((1 - sigma0Samples * sigma1Samples)+
    alpha1Samples *(1-sigma1Samples) + alpha0Samples*(1-sigma0Samples ))/
  ((1+alpha1Samples)*(1+alpha0Samples)) 
