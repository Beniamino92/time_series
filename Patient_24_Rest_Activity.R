################# Spectral Analysis Rest Activity ##################
#################       Patient 24                 ##################

require("astsa")


# setwd("C:/Users/Beniamino/Desktop/Project_8")
setwd("/homes/hadjamar/Documents/Project_8")
# setwd("/home/hadjamar/Desktop/Project_8")

source("analysis_healthy_patients.R")

# Setting times and frequencies

T <- length(RA24) # 97
t <- 1:T
freq <- 0:((T-1)/2 - 1)/T

# Adding a little bit of (positive) noise
# N.zeros <- length(which(RA8 == 0))
# noise <- abs(rnorm(N.zeros, 0, 1))
# RA8[which(RA8 == 0)] <- RA8[which(RA8 == 0)] + noise
