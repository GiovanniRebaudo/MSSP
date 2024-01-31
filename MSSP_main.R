# Codes accompanying "Multivariate Species Sampling Process (MSSP)"
# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.14
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(salso)   # version 0.3.29
library(ggplot2) # version 3.4.2

# Load functions
source("MSSP_fcts.R")
Save_Plot = F


