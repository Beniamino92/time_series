################# Spectral Analysis Rest Activity ##################
#################       Patient 8                 ##################

require("astsa")
library("astsa")

# setwd("C:/Users/Beniamino/Desktop/Project_8")
# setwd("/homes/hadjamar/Documents/Project_8")
setwd("/home/hadjamar/Desktop/Project_8")

source("analysis_healthy_patients.R")




# Setting times and frequencies
T <- length(RA8) # 97
t <- 1:T
freq <- 0:((T-1)/2 - 1)/T

# Adding a little bit of (positive) noise
noise <- rexp(T, 10)
RA8 <- RA8 + noise

# Plot Time Series
plot.ts(RA8, type = "o", pch = 19)

# It's better to do spectral analysis on the residuals (i.e considering a 
# zero mean process)
res.RA8 <- lm(RA8 ~ t)$residuals
plot.ts(res.RA8, ylab = "Residual Rest Activity", 
        main = "Residual Rest Activity [8]")
plot.ts(res.RA8, ylab = "Residual Rest Activity", 
        main = "Residual Rest Activity [8]", 
        type = "o", pch = 19)


# Getting and plotting spectrum of periodogram, smoothed periodogram (unif weights),
# smoothed periodogram (daniell weights))
spectrums <- make_periodograms(res.RA8, freq, plot = TRUE)

periodogram <- spectrums$periodogram
fhat.unif <- spectrums$fhat.unif
fhat.daniel <- spectrums$fhat.daniel

# Now I want to find the frequency for which we have 
# maximum value on the periodogram

driving.freq <- freq[which(fhat.daniel
                           == max(fhat.daniel))+ 1] # 0.04123711
driving.freq



### Exploring frequency domain trough AR(p)

# Look at PACF
acf(res.RA8, lag = 31, main = "ACFc")
pacf(res.RA8, lag = 31, main = "PACF")     # AR(10) ?

# Plotting AIC & BIC given by the fit of AR(p), as p increases
performance.ar.p(res.RA8, 20, plot = TRUE)





