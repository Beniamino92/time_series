################# Spectral Analysis Rest Activity ##################
#################       Patient 24                 ##################

require("astsa")


# setwd("C:/Users/Beniamino/Desktop/Project_8")
setwd("/homes/hadjamar/Documents/Project_8")
# setwd("/home/hadjamar/Desktop/Project_8")

source("analysis_healthy_patients.R")

# Setting times and frequencies

T <- length(RA24); T # 97
t <- 1:T
freq <- 0:((T-1)/2 - 1)/T

# Adding a little bit of (positive) noise
# N.zeros <- length(which(RA8 == 0))
# noise <- abs(rnorm(N.zeros, 0, 1))
# RA8[which(RA8 == 0)] <- RA8[which(RA8 == 0)] + noise

# Plot Time Series
plot.ts(RA24, type = "o", pch = 19)


# It's better to do spectral analysis on the residuals (i.e considering a 
# zero mean process)
res.RA24 <- lm(RA24 ~ t)$residuals
plot.ts(res.RA24, ylab = "Residual Rest Activity", 
        main = "Residual Rest Activity [24]")
plot.ts(res.RA24, ylab = "Residual Rest Activity", 
        main = "Residual Rest Activity [24]", 
        type = "o", pch = 19)



# Getting and plotting spectrum of periodogram, smoothed periodogram (unif weights),
# smoothed periodogram (daniell weights))
spectrums <- make_periodograms(res.RA24, freq, plot = TRUE)

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
acf(res.RA24, lag = 31, main = "ACFc")
pacf(res.RA24, lag = 31, main = "PACF")     # AR(1) AR(6) AR(7)?

# Plotting AIC & BIC given by the fit of AR(p), as p increases
performance.ar.p(res.RA24, 30, plot = TRUE)

# It looks possibile AR(1) 





#### Comparing periodogram with AR spectrum, by choosing AR(1) ##############

p <- 1

# Fitting the time series to an AR(1) and getting the coefficients
fit1 <- arima(res.RA24, method = "CSS-ML", order = c(p, 0, 0))
coef1 <- fit1$coef[1:p]

# Getting the true spectrum of an AR(1) with that coefficients
ar1.fit <- arma.spec(ar = coef1, n.freq = length(freq) + 1,
                     var.noise = fit1$sigma2)$spec

# Remove the estimate at frequency 0
ar1.fit <- ar1.fit[-1]

# Plotting periodogram and spectrum AR(1)
plot(freq, periodogram, type = "l", ylab = "Power spectrum", xlab = "Frequency",
     lwd = 2)
lines(freq, ar1.fit, col = "red", lwd = 2)
legend("topright", c("Periodogram", "Spectrum AR(1)"), col = c("black", "red"),
       lty = 1)





#################### FOURIER REGRESSION ###################

# Selecting the four frequencies, in order to obtain the four

driving.frequencies <- c()

number.harmonics <- 6

# Ranking the frequencies
rank.freq <- sort(fhat.daniel, decreasing = TRUE )
rank.freq

for(i in 1:number.harmonics) {
  driving.frequencies[i] <- freq[which(fhat.daniel == rank.freq[i]) + 1]
}

# Getting the harmonics 
harmonics <- list()

for(i in 1:number.harmonics) {
  harmonics[[i]] <- get_harmonic(res.RA24, driving.frequencies[i])
}

# Getting the trend
t <- 1:T
trend.RA24<- as.vector(fitted(lm(RA24 ~ t)))


########## Final Model

model.RA24 <- get_model(harmonics, T) + trend.RA24

# Setting zero when the fit is negative:
model.RA24[which(model.RA24 < 0)] <- 0


# Plotting time series, adding single harmonics + final model
plot.ts(RA24, type = "o", pch = 19, main = "Rest Activity Patient 24",
        ylab = "Temperature", ylim = c(min(model.RA24), max(RA24)))
lines(1:T, model.RA24, col = "red", lwd = 5, type = "o")



### Residual Analysis

residuals <- RA24 - model.RA24
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

Forecast24h <- (model.RA24[1:24] + model.RA24[25:48] + 
                  model.RA24[49:72] + model.RA24[73:96])/4

Residuals24h <- (residuals[1:24] + residuals[25:48] + 
                   residuals[49:72] + residuals[73:96])/4

upper.CI <- Forecast24h + 1.96*sd(Residuals24h)
lower.CI <- Forecast24h - 1.96*sd(Residuals24h)

lower.CI[which(lower.CI < 0)] <- 0

# Plotting forecasting and relative confidence intervals
plot.ts(RA24, lwd = 1, type = "o", pch = 19,
        xlim = c(1, T + 24), ylab = "Temperature", ylim = c(-5, 90))
lines(1:T, model.RA24, col = "red", lwd = 4, type = "o")
lines((T+1):(T + 24), Forecast24h, col = "blue", lwd = 4, type = "l")
lines((T+1):(T + 24), upper.CI,
      col = "blue", lwd = 2, lty = 3)
lines((T+1):(T + 24), lower.CI,
      col = "blue", lwd = 2, lty = 3)





#### Testing the model:

# This consist of fitting the model (using harmonic regression) on 3 days,
# and testing on the 4th.

training.set <- RA24[1:72]
test.set <- RA24[73:96]

T <- length(training.set)
t <- 1:T
freq <- 0:((T-1)/2 - 1)/T

res.RA24 <- lm(training.set ~ t)$residuals
trend.RA24 <- as.vector(fitted(lm(training.set ~ t)))

# Harmonic Regression (using the frequencies ) 
harmonics <- list()
for(i in 1:number.harmonics) {
  harmonics[[i]] <- get_harmonic(res.RA24, driving.frequencies[i])
}


# I decided to use first and second harmonics, otherwise it's gonna
# over fit. (i.e., I get a bigger testing error)
model.RA24 <- get_model(harmonics, T)+ trend.RA24

model.RA24[which(model.RA24 < 0)] <- 0

forecast24h <- (model.RA24[1:24] + model.RA24[25:48] + model.RA24[49:72])/3 

sum((forecast24h - test.set)^2)/24 # 430.8416



# Just to see what's going on
plot.ts(RA24, type = "o", pch = 19)
lines(1:T, model.RA24, lwd = 4, col = "red")
lines(73:96, forecast24h, col = "blue", lwd = 4, type = "o")
