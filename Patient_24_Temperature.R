################# Spectral Analysis Temperature ##################
#################       Patient 24               ##################


require("astsa")

# setwd("C:/Users/Beniamino/Desktop/Project_8")
#setwd("/home/hadjamar/Desktop/Project_8")
setwd("/homes/hadjamar/Documents/Project_8/")

source("analysis_healthy_patients.R")

### Let's work with Patient 24 ###

T <- length(Temp24) # 97
t <- 1:T
freq <- 0:((T-1)/2 - 1)/T

# It's better to do spectral analysis on the residuals (i.e considering a 
# zero mean process)

res.Temp24 <- lm(Temp24 ~ t)$residuals
plot.ts(res.Temp24, ylab = "Residual Temperature", 
        main = "Residual Temperature [2]")

M <- 3

# Getting and plotting spectrum of periodogram, smoothed periodogram (unif weights),
# smoothed periodogram (daniell weights))
spectrums <- make_periodograms(res.Temp24, freq, plot = TRUE)

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
acf(res.Temp24, lag = 31, main = "ACFc") 
pacf(res.Temp24, lag = 31, main = "PACF") # AR(11)

# Plotting AIC & BIC given by the fit of AR(p), as p increases
performance.ar.p(res.Temp24, 20, plot = TRUE)

# It looks AR(3) or AR(16)


############ Comparing periodogram with AR spectrum, by choosing AR(3) ##############

p <- 3

# Fitting the time series to an AR(3) and getting the coefficients
fit3 <- arima(res.Temp24, method = "CSS-ML", order = c(p, 0, 0))
coef3 <- fit3$coef[1:p]

# Getting the true spectrum of an AR(3) with that coefficients
ar3.fit <- arma.spec(ar = coef3, n.freq = length(freq) + 1,
                     var.noise = fit3$sigma2)$spec

# Remove the estimate at frequency 0]
ar3.fit <- ar3.fit[-1]

# Plotting periodogram and spectrum AR(2)
plot(freq, periodogram, type = "l", ylab = "Power spectrum", xlab = "Frequency", lwd  = 2)
lines(freq, ar3.fit, col = "red", lwd = 2)
legend("topright", c("Periodogram", "Spectrum AR(3)"), col = c("black", "red"),
       lty = 1)


############ Comparing periodogram with AR spectrum, by choosing AR(16) ##############

p <- 16

# Fitting the time series to an AR(16) and getting the coefficients
fit16 <- arima(res.Temp24, method = "CSS-ML", order = c(p, 0, 0))
coef16 <- fit16$coef[1:p]

# Getting the true spectrum of an AR(16) with that coefficients
ar16.fit <- arma.spec(ar = coef16, n.freq = length(freq) + 1,
                     var.noise = fit16$sigma2)$spec

# Remove the estimate at frequency 0
ar16.fit <- ar16.fit[-1]

# Plotting periodogram and spectrum AR(16)
plot(freq, periodogram, type = "l", ylab = "Power spectrum", xlab = "Frequency", lwd = 2, ylim = c(0, 10))
lines(freq, ar16.fit, col = "red", lwd = 2)
legend("topright", c("Periodogram", "Spectrum AR(16)"), col = c("black", "red"),
       lty = 1)




# Plotting Periodogram, AR(3), AR(16)
plot(freq, periodogram, type = "l", ylab = "Power spectrum", xlab = "Frequency", 
     ylim = c(0, 9), lwd = 2)
lines(freq, ar3.fit, col = "blue", lwd = 3)
lines(freq, ar16.fit, col = "red", lwd = 3)
legend("topright", c("Periodogram", "Spectrum AR(3)", "Spectrum AR(16)"), col = c("black", "blue", "red"),
       lty = 1, lwd = 3)





#################### FOURIER REGRESSION ###################

# Selecting the four frequencies, in order to obtain the four

driving.frequencies <- c()

number.harmonics <- 7

# Ranking the frequencies
rank.freq <- sort(fhat.daniel, decreasing = TRUE )
rank.freq

for(i in 1:number.harmonics) {
  driving.frequencies[i] <- freq[which(fhat.daniel == rank.freq[i]) + 1]
}

# Getting the harmonics 
harmonics <- list()

for(i in 1:number.harmonics) {
  harmonics[[i]] <- get_harmonic(res.Temp24, driving.frequencies[i])
}

# Getting the trend
t <- 1:T
trend.Temp24 <- as.vector(fitted(lm(Temp24 ~ t)))


# Final Model

model.Temp24 <- get_model(harmonics, T) + trend.Temp24

# Plotting time series, adding single harmonics + final model
plot.ts(Temp24, type = "o", pch = 19, main = "Temperature Patient 24",
        ylab = "Temperature")
# lines(1:T, harmonics[[1]] + mean.RA8, col = "blue", lwd = 3, lty = 3)
# lines(1:T, harmonics[[2]] + mean.RA8, col = "chartreuse3", lwd = 3, lty = 3)
# lines(1:T, harmonics[[4]] + mean.RA8, col = "grey", lwd = 3, lty = 3)
lines(1:T, model.Temp24, col = "red", lwd = 5)



