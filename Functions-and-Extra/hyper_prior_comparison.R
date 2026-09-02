#fix hyperprior 

set.seed(0)
mc_samples = 10000

## DP #################################################################
alpha1Samples = rgamma(mc_samples, shape = 0.75, rate = 1)

probTieWithinSamples_PY = 
  1  / (1 + alpha1Samples)

mean(probTieWithinSamples_PY)
var(probTieWithinSamples_PY)

## PY #################################################################
alpha1Samples = rgamma(mc_samples, shape = 0.2, rate = 1)
sigma1Samples = rbeta(mc_samples, 1, 3)

probTieWithinSamples_PY = 
  (1 - sigma1Samples)  / (1 + alpha1Samples)

mean(probTieWithinSamples_PY)
var(probTieWithinSamples_PY)

## +DP #################################################################
eps1 = sample(c(0,1,2), mc_samples, replace = TRUE, prob = c(0.02,0.28,0.7))
eps1[eps1 == 2] = rbeta(sum(eps1 == 2), 1, 1)
eps2 = sample(c(0,1,2), mc_samples, replace = TRUE, prob = c(0.02,0.28,0.7))
eps2[eps2 == 2] = rbeta(sum(eps2 == 2), 1, 1)
alpha0Samples = rgamma(mc_samples, shape = 1, rate = 20)
alpha1Samples = rgamma(mc_samples, shape = 6, rate = 10)

probTieAcrossSamples_plusDP = 
  (eps1 * eps2) / (1 + alpha0Samples)

mean(probTieAcrossSamples_plusDP)
var(probTieAcrossSamples_plusDP)

probTieWithinSamples_plusDP = 
  eps1^2 /  (1 + alpha0Samples) + (1 - eps1)^2 / (1 + alpha1Samples)

mean(probTieWithinSamples_plusDP)
var(probTieWithinSamples_plusDP)

## +PY #################################################################
eps1 = sample(c(0,1,2), mc_samples, replace = TRUE, 
              prob = c(0.02,0.28,0.7))
eps1[eps1 == 2] = rbeta(sum(eps1 == 2), 1, 1)

eps2 = sample(c(0,1,2), mc_samples, replace = TRUE, 
              prob = c(0.02,0.28,0.7))
eps2[eps2 == 2] = rbeta(sum(eps2 == 2), 1, 1)

alpha0Samples = rgamma(mc_samples, shape = 1, rate = 20)
alpha1Samples = rgamma(mc_samples, shape = 1, rate = 10)
sigma0Samples = rbeta(mc_samples, .5, 20)
sigma1Samples = rbeta(mc_samples, 1, 3)

probTieAcrossSamples_plusPY = 
  (eps1 * eps2) * (1 - sigma0Samples) / (1 + alpha0Samples)

mean(probTieAcrossSamples_plusPY)
var(probTieAcrossSamples_plusPY)

probTieWithinSamples_plusPY = 
  eps1^2 * (1 - sigma0Samples)  /  (1 + alpha0Samples) + 
  (1 - eps1)^2 * (1 - sigma1Samples)  / (1 + alpha1Samples)

mean(probTieWithinSamples_plusPY)
var(probTieWithinSamples_plusPY)

## HDP #################################################################
alpha0Samples = rgamma(mc_samples, shape = 1, rate = 1/3)
alpha1Samples = rgamma(mc_samples, shape = 1, rate = 1/2)

probTieAcrossSamples_HDP = 1 / (1 + alpha0Samples)

mean(probTieAcrossSamples_HDP)
var(probTieAcrossSamples_HDP)

probTieWithinSamples_HDP = 
  (1 +  alpha1Samples + alpha0Samples) / 
  (1 + alpha1Samples) / (1 + alpha0Samples)

mean(probTieWithinSamples_HDP)
var(probTieWithinSamples_HDP)

## HPY #################################################################
alpha0Samples = rgamma(mc_samples, shape = 1, rate = 1)
alpha1Samples = rgamma(mc_samples, shape = 1, rate = 1)
sigma0Samples = rbeta(mc_samples, 1, 2)
sigma1Samples = rbeta(mc_samples, 1, 2)

probTieAcrossSamples_HPY = 
  (1 - sigma0Samples) / (1 + alpha0Samples)

mean(probTieAcrossSamples_HPY)
var(probTieAcrossSamples_HPY)

probTieWithinSamples_HPY = 
  ( (1 - sigma1Samples*sigma0Samples) + 
      alpha1Samples*(1 - sigma0Samples) + 
      alpha0Samples*(1 - sigma1Samples) ) / 
  (1 + alpha1Samples) / (1 + alpha0Samples)

mean(probTieWithinSamples_HPY)
var(probTieWithinSamples_HPY)
