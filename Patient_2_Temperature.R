################# Spectral Analysis Temperature ##################
#################       Patient 2               ##################

require("astsa")

# setwd("C:/Users/Beniamino/Desktop/Project_8")
setwd("/home/hadjamar/Desktop/Project_8")
#setwd("/homes/hadjamar/Documents/Project_8/")

source("analysis_healthy_patients.R")

### Let's work with Patient 8 ###

T <- length(Temp2) # 128
t <- 1:T

freq <- 0:((T/2 - 1))/T

# It's better to do spectral analysis on the residuals (i.e considering a 
# zero mean process)

res.Temp2 <- lm(Temp2 ~ t)$residuals
plot.ts(res.Temp2, ylab = "Residual Temperature", 
        main = "Residual Temperature [2]")

M <- 3

# Getting and plotting spectrum of periodogram, smoothed periodogram (unif weights),
# smoothed periodogram (daniell weights))
spectrums <- make_periodograms(res.Temp2, freq, plot = TRUE)

periodogram <- spectrums$periodogram
fhat.unif <- spectrums$fhat.unif
fhat.daniel <- spectrums$fhat.daniel

# Now I want to find the frequency for which we have 
# maximum value on the periodogram

driving.freq <- freq[which(periodogram
                           == max(periodogram))+ 1] # 0.0390625
driving.freq


### Exploring frequency domain trough AR(p)

# Look at PACF
acf(res.Temp2, lag = 31, main = "ACFc") 
pacf(res.Temp2, lag = 31, main = "PACF") # AR(1)

# Plotting AIC & BIC given by the fit of AR(p), as p increases
performance.ar.p(res.Temp2, 20, plot = TRUE)

# It looks AR(1)




#### Comparing periodogram with AR spectrum, by choosing AR(1) ##############

p <- 1

# Fitting the time series to an AR(1) and getting the coefficients
fit1 <- arima(res.Temp2, method = "CSS-ML", order = c(p, 0, 0))
coef1 <- fit1$coef[1:p]

# Getting the true spectrum of an AR(1) with that coefficients
ar1.fit <- arma.spec(ar = coef1, n.freq = length(freq) + 1,
                     var.noise = fit1$sigma2)$spec

# Remove the estimate at frequency 0
ar1.fit <- ar1.fit[-1]



# Plotting Periodogram, AR(1)
plot(freq, periodogram, type = "l", ylab = "Power spectrum", xlab = "Frequency",
     ylim = c(0, max(max(periodogram), max(periodogram))))
lines(freq, ar1.fit, col = "blue", lwd = 3)
legend("topright", c("Periodogram", "Spectrum AR(1)"), col = c("black", "blue"),
       lty = 1, lwd = 3)


#################### FOURIER REGRESSION ###################

# Selecting the four frequencies, in order to obtain the four harmonics later

driving.frequencies <- c()

number.harmonics <- 5

# Ranking the frequencies
rank.freq <- sort(periodogram, decreasing = TRUE )

for(i in 1:number.harmonics) {
  driving.frequencies[i] <- freq[which(periodogram == rank.freq[i]) + 1]
}

# Getting the harmonics 
harmonics <- list()

for(i in 1:number.harmonics) {
  harmonics[[i]] <- get_harmonic(res.Temp2, driving.frequencies[i])
}

# Getting the trend
t <- 1:T
trend.Temp2 <- as.vector(fitted(lm(Temp2 ~ t)))

# Final Model

model.Temp2 <- get_model(harmonics, T) + trend.Temp2
# model.RA8 <- harmonics[[1]] + harmonics[[2]] + trend.Temp2

# Plotting time series, adding single harmonics + final model
plot.ts(Temp2, type = "o", pch = 19, main = "Temperature Patient 2",
        ylab = "Temperature")
# lines(1:T, harmonics[[1]] + mean.RA8, col = "blue", lwd = 3, lty = 3)
# lines(1:T, harmonics[[2]] + mean.RA8, col = "chartreuse3", lwd = 3, lty = 3)
# lines(1:T, harmonics[[4]] + mean.RA8, col = "grey", lwd = 3, lty = 3)
lines(1:T, model.Temp2, col = "red", lwd = 5)




### Residual Analysis
residuals <- Temp2 - model.Temp2
std.residuals <- residuals/sd(residuals)

# Standardized residuals plot
plot(1:T, std.residuals, pch = 19, col = "red")
abline(h = 0, col = "black", lty = 2)

# QQPlot

grid <- -300:300/100
qqnorm(std.residuals, ylim = c(-3, 3), xlim = c(-3, 3))
lines(grid, grid)

# The residuals are not exactly Normal...but that's what we have.




#### Testing the model:

# This consist of fitting the model (using harmonic regression) on 3 days,
# and testing on the 4th.

training.set <- Temp2[1:72]
test.set <- Temp2[72:95]

T <- length(training.set)
t <- 1:T
freq <- 0:((T-1)/2 - 1)/T

res.Temp2 <- lm(training.set ~ t)$residuals
trend.Temp2 <- as.vector(fitted(lm(training.set ~ t)))

# Harmonic Regression (using the frequencies ) 
harmonics <- list()
for(i in 1:5) {
  harmonics[[i]] <- get_harmonic(res.Temp2, driving.frequencies[i])
}
model.Temp2 <- get_model(harmonics, T) + trend.Temp2
forecast24h <- (model.Temp2[1:24] + model.Temp2[25:48] + model.Temp2[49:72])/3 

sum((forecast24h - test.set)^2)/24 # 1.919622






