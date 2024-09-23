# Check and evaluate efficiency R code
# rm(list=ls())
library(microbenchmark)
library(orthopolynom)

pochhammer(1.9, 10)

gamma(11.9)/gamma(1.9)

microbenchmark(pochhammer(1.9, 10),
               gamma(11.9)/gamma(1.9))

library(boot)
inv.logit(logit_sig_prop)

