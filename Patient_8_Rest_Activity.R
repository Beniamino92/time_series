################# Spectral Analysis Rest Activity ##################
#################       Patient 8                 ##################

require("astsa")


# setwd("C:/Users/Beniamino/Desktop/Project_8")
# setwd("/homes/hadjamar/Documents/Project_8")
setwd("/home/hadjamar/Desktop/Project_8")

source("analysis_healthy_patients.R")




# Setting times and frequencies
T <- length(RA8) # 97
t <- 1:T
freq <- 0:((T-1)/2 - 1)/T

# Adding a little bit of (positive) noise
# N.zeros <- length(which(RA8 == 0))
# noise <- abs(rnorm(N.zeros, 0, 1))
# RA8[which(RA8 == 0)] <- RA8[which(RA8 == 0)] + noise

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
performance.ar.p(res.RA8, 30, plot = TRUE)

# It looks possibile either AR(1) or AR(26)


#### Comparing periodogram with AR spectrum, by choosing AR(1) ##############

p <- 1

# Fitting the time series to an AR(1) and getting the coefficients
fit1 <- arima(res.RA8, method = "CSS-ML", order = c(p, 0, 0))
coef1 <- fit1$coef[1:p]

# Getting the true spectrum of an AR(1) with that coefficients
ar1.fit <- arma.spec(ar = coef1, n.freq = length(freq) + 1,
                     var.noise = fit1$sigma2)$spec

# Remove the estimate at frequency 0
ar1.fit <- ar1.fit[-1]

# Plotting periodogram and spectrum AR(1)
plot(freq, periodogram, type = "l", ylab = "Power spectrum", xlab = "Frequency")
lines(freq, ar1.fit, col = "red", lwd = 2)
legend("topright", c("Periodogram", "Spectrum AR(1)"), col = c("black", "red"),
       lty = 1)


##### Comparing periodogram with AR spectrum, by choosing AR(9) ##############

p <- 26

# Fitting the time series to an AR(26) and getting the coefficients
fit26 <- arima(res.RA8, method = "CSS-ML", order = c(p, 0, 0))
coef26 <- fit26$coef[1:p]

# Getting the true spectrum of an AR(26) with that coefficients
ar26.fit <- arma.spec(ar = coef26, n.freq = length(freq) + 1,
                     var.noise = fit26$sigma2)$spec

# Remove the estimate at frequency 0
ar26.fit <- ar26.fit[-1]

# Plotting periodogram and spectrum AR(26)
plot(freq, periodogram, type = "l", ylab = "Power spectrum", xlab = "Frequency", 
     ylim = c(0, max(max(periodogram), max(ar26.fit))))
lines(freq, ar26.fit, col = "red", lwd = 2)
legend("topright", c("Periodogram", "Spectrum AR(26)"), col = c("black", "red"),
       lty = 1)




# Plotting Periodogram, AR(1), AR(26)
plot(freq, periodogram, type = "l", ylab = "Power spectrum", xlab = "Frequency",
     ylim = c(0, max(max(periodogram), max(ar26.fit))))
lines(freq, ar1.fit, col = "blue", lwd = 3)
lines(freq, ar26.fit, col = "red", lwd = 3)
legend("topright", c("Periodogram", "Spectrum AR(1)", "Spectrum AR(26)"), col = c("black", "blue", "red"),
       lty = 1, lwd = 3)





#################### FOURIER REGRESSION ###################

# Selecting the four frequencies, in order to obtain the four harmonics later

driving.frequencies <- c()

number.harmonics <- 5

# Ranking the frequencies
rank.freq <- sort(ar26.fit, decreasing = TRUE )

for(i in 1:number.harmonics) {
  driving.frequencies[i] <- freq[which(ar26.fit == rank.freq[i]) + 1]
}

driving.frequencies

# Getting the harmonics 
harmonics <- list()

for(i in 1:number.harmonics) {
  harmonics[[i]] <- get_harmonic(RA8, driving.frequencies[i])
}

# Getting the trend
t <- 1:T
trend.RA8 <- as.vector(fitted(lm(RA8 ~ t)))

# Mean 
mean.RA8 <- mean(RA8[-(which(RA8 == max(RA8)))])

# Final Model
# model.RA8 <- get_model(harmonics, T) + trend.RA8
model.RA8 <- harmonics[[1]] + harmonics[[2]] + harmonics[[4]]  + mean.RA8

# Plotting time series, adding single harmonics + final model
plot.ts(RA8, type = "o", pch = 19, main = "Rest Activty Patient 8", ylab = "Temperature", 
        ylim = c(-4, max(RA8)))
# lines(1:T, harmonics[[1]] + mean.RA8, col = "blue", lwd = 3, lty = 3)
# lines(1:T, harmonics[[2]] + mean.RA8, col = "chartreuse3", lwd = 3, lty = 3)
# lines(1:T, harmonics[[4]] + mean.RA8, col = "grey", lwd = 3, lty = 3)
lines(1:T, model.RA8, col = "red", lwd = 5)

# legend("topright", 
#       c(expression(paste(omega, " = 1/24 ")), 
#         expression(paste(omega, " = 1/12 ")),
#         expression(paste(omega, " = 1/8 "))), 
#       lty = 3, col = c("blue", "chartreuse3", 
#                        "darkorange1"), lwd = 3)
# legend("topleft", "Total Model", lwd = 5, col = "red", lty = 1)



### Residual Analysis
residuals <- RA8 - model.RA8
std.residuals <- residuals/sd(residuals)

# Standardized residuals plot
plot(1:T, std.residuals, pch = 19, col = "red")
abline(h = 0, col = "black", lty = 2)

# QQPlot

grid <- -300:300/100
qqnorm(std.residuals, ylim = c(-3, 3), xlim = c(-3, 3))
lines(grid, grid)




########### FORECASTING #############


# My forecast for the next 24 hours, it's made by an average
# of the final model we fitted, over 24 hours.
# That is:

Forecast24h <- (model.RA8[1:24] + model.RA8[25:48] + 
                  model.RA8[49:72] + model.RA8[73:96])/4

Residuals24h <- (residuals[1:24] + residuals[25:48] + 
                   residuals[49:72] + residuals[73:96])/4

upper.CI <- Forecast24h + 1.96*sd(Residuals24h)
lower.CI <- Forecast24h - 1.96*sd(Residuals24h)


# Plotting forecasting and relative confidence intervals
plot.ts(RA8, lwd = 1, type = "o", pch = 19,
        xlim = c(1, T + 24), ylab = "Rest Activity", ylim = c(min(lower.CI), max(RA8)))
lines(1:T, model.RA8, col = "red", lwd = 4)
lines(T:(T + 24), c(model.RA8[T], Forecast24h), col = "blue", lwd = 4)
lines(T:(T + 24), c(model.RA8[T], upper.CI),
      col = "blue", lwd = 2, lty = 3)
lines(T:(T + 24), c(model.RA8[T], lower.CI),
      col = "blue", lwd = 2, lty = 3)






#### Testing the model:

# This consist of fitting the model (using harmonic regression) on 3 days,
# and testing on the 4th.

training.set <- RA8[1:72]
test.set <- RA8[72:95]

T <- length(training.set)
t <- 1:T
freq <- 0:((T-1)/2 - 1)/T

res.RA8 <- lm(training.set ~ t)$residuals
trend.RA8 <- as.vector(fitted(lm(training.set ~ t)))

# Harmonic Regression (using the frequencies ) 
harmonics <- list()
for(i in 1:number.harmonics) {
  harmonics[[i]] <- get_harmonic(res.RA8, driving.frequencies[i])
}

model.RA8 <- harmonics[[1]] + harmonics[[2]] + harmonics[[4]]  + trend.RA8

forecast24h <- (model.RA8[1:24] + model.RA8[25:48] + model.RA8[49:72])/3 

sum((forecast24h - test.set)^2)/24 # 141.3878

# Fitting the model just by using 3 days, and testing to the 4th, doesn't 
# look good. 


# Just to see what's going on
plot.ts(RA8, type = "o", pch = 19, ylim = c(-10, 70))
lines(1:T, model.RA8, lwd = 4, col = "red")
lines(72:95, forecast24h, col = "blue", lwd = 4)
