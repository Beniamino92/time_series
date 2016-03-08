################# Spectral Analysis Temperature ##################
#################       Patient 24               ##################


require("astsa")

# setwd("C:/Users/Beniamino/Desktop/Project_8")
#setwd("/home/hadjamar/Desktop/Project_8")
setwd("/homes/hadjamar/Documents/Project_8/")

source("analysis_healthy_patients.R")

### Let's work with Patient 24 ###

T <- length(Temp26) # 97
t <- 1:T
freq <- 0:((T-1)/2 - 1)/T

# It's better to do spectral analysis on the residuals (i.e considering a 
# zero mean process)

res.Temp26 <- lm(Temp26 ~ t)$residuals
plot.ts(res.Temp26, ylab = "Residual Temperature", 
        main = "Residual Temperature [2]")

M <- 3

# Getting and plotting spectrum of periodogram, smoothed periodogram (unif weights),
# smoothed periodogram (daniell weights))
spectrums <- make_periodograms(res.Temp26, freq, plot = TRUE)

periodogram <- spectrums$periodogram
fhat.unif <- spectrums$fhat.unif
fhat.daniel <- spectrums$fhat.daniel


# Now I want to find the frequency for which we have 
# maximum value on the periodogram

driving.freq <- freq[which(periodogram
                           == max(periodogram))+ 1] # 0.04123711
driving.freq



### Exploring frequency domain trough AR(p)

# Look at PACF
acf(res.Temp26, lag = 31, main = "ACFc") 
pacf(res.Temp26, lag = 31, main = "PACF") # AR(1) or AR(22)

# Plotting AIC & BIC given by the fit of AR(p), as p increases
performance.ar.p(res.Temp26, 20, plot = TRUE)

# It looks AR(2)


############ Comparing periodogram with AR spectrum, by choosing AR(3) ##############

p <- 2

# Fitting the time series to an AR(3) and getting the coefficients
fit2 <- arima(res.Temp26, method = "CSS-ML", order = c(p, 0, 0))
coef2 <- fit2$coef[1:p]

# Getting the true spectrum of an AR(3) with that coefficients
ar2.fit <- arma.spec(ar = coef2, n.freq = length(freq) + 1,
                     var.noise = fit2$sigma2)$spec

# Remove the estimate at frequency 0]
ar2.fit <- ar2.fit[-1]

# Plotting periodogram and spectrum AR(2)
plot(freq, periodogram, type = "l", ylab = "Power spectrum", xlab = "Frequency", lwd  = 2)
lines(freq, ar2.fit, col = "red", lwd = 2)
legend("topright", c("Periodogram", "Spectrum AR(2)"), col = c("black", "red"),
       lty = 1)







#################### FOURIER REGRESSION ###################

# Selecting the four frequencies, in order to obtain the four

driving.frequencies <- c()

number.harmonics <- 3

# Ranking the frequencies
rank.freq <- sort(periodogram, decreasing = TRUE )
rank.freq

for(i in 1:number.harmonics) {
  driving.frequencies[i] <- freq[which(periodogram == rank.freq[i]) + 1]
}

# Getting the harmonics 
harmonics <- list()

for(i in 1:number.harmonics) {
  harmonics[[i]] <- get_harmonic(res.Temp26, driving.frequencies[i])
}

# Getting the trend
t <- 1:T
trend.Temp26 <- as.vector(fitted(lm(Temp26 ~ t)))


# Final Model

model.Temp26 <- get_model(harmonics, T) + trend.Temp26

# Plotting time series, adding single harmonics + final model
plot.ts(Temp26, type = "o", pch = 19, main = "Temperature Patient 26",
        ylab = "Temperature")
# lines(1:T, harmonics[[1]] + mean.RA8, col = "blue", lwd = 3, lty = 3)
# lines(1:T, harmonics[[2]] + mean.RA8, col = "chartreuse3", lwd = 3, lty = 3)
# lines(1:T, harmonics[[4]] + mean.RA8, col = "grey", lwd = 3, lty = 3)
lines(1:T, model.Temp26, col = "red", lwd = 5)



### Residual Analysis

residuals <- Temp26 - model.Temp26 
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

Forecast24h <- (model.Temp26[1:24] + model.Temp26[25:48] + 
                  model.Temp26[49:72] + model.Temp26[73:96])/4

Residuals24h <- (residuals[1:24] + residuals[25:48] + 
                   residuals[49:72] + residuals[73:96])/4

upper.CI <- Forecast24h + 1.96*sd(Residuals24h)
lower.CI <- Forecast24h - 1.96*sd(Residuals24h)


# Plotting forecasting and relative confidence intervals
plot.ts(Temp26, lwd = 1, type = "o", pch = 19,
        xlim = c(1, T + 24), ylab = "Temperature")
lines(1:T, model.Temp26, col = "red", lwd = 4, type = "o")
lines(T:(T + 24), c(model.Temp26[T], Forecast24h), col = "blue", lwd = 4, type = "o")
lines((T+1):(T + 24), upper.CI,
      col = "blue", lwd = 2, lty = 3)
lines((T+1):(T + 24), lower.CI,
      col = "blue", lwd = 2, lty = 3)






#### Testing the model:

# This consist of fitting the model (using harmonic regression) on 3 days,
# and testing on the 4th.

training.set <- Temp26[1:72]
test.set <- Temp26[73:96]

T <- length(training.set)
t <- 1:T
freq <- 0:((T/2 - 1))/T

res.Temp26 <- lm(training.set ~ t)$residuals
trend.Temp26 <- as.vector(fitted(lm(training.set ~ t)))

# Harmonic Regression (using the frequencies ) 
harmonics <- list()
for(i in 1:number.harmonics) {
  harmonics[[i]] <- get_harmonic(res.Temp26, driving.frequencies[i])
}
model.Temp26 <- get_model(harmonics, T) + trend.Temp26
forecast24h <- (model.Temp26[1:24] + model.Temp26[25:48] + model.Temp26[49:72])/3 

sum((forecast24h - test.set)^2)/24 # 0.4278687

# Just to see what's going on
plot.ts(Temp26[1:96], type = "o", pch = 19)
lines(1:T, model.Temp26, lwd = 4, col = "red", type = "o")
lines(73:96, forecast24h, col = "blue", lwd = 4, type = "o")
