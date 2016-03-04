################# Spectral Analysis Temperature ##################
#################       Patient 8               ##################

require("astsa")

# setwd("C:/Users/Beniamino/Desktop/Project_8")
setwd("/homes/hadjamar/Documents/Project_8")

source("analysis_healthy_patients.R")

### Let's work with Patient 8 ###

T <- length(Temp8) # 95
t <- 1:T

freq <- 0:((T-1)/2 - 1)/T

# It's better to do spectral analysis on the residuals (i.e considering a 
# zero mean process)

res.Temp8 <- lm(Temp8 ~ t)$residuals
plot.ts(res.Temp8, ylab = "Residual Temperature", 
        main = "Residual Temperature [2]")

M <- 2

# Periodogram 
periodogram <- spec.pgram(res.Temp8, fast = FALSE, taper=0,
                          log="no", detrend=TRUE, plot=FALSE)

# Smoothed Periodogram (Uniform weights)
fhat.unif <- spec.pgram(res.Temp8, fast = FALSE, taper=0,
                        log="no", detrend=TRUE, plot=FALSE,
                        kernel = kernel(rep(1, M +1)/(2*M + 1)))

# Smoothed Periodogram (Daniell modified weights)
fhat.daniel <- spec.pgram(res.Temp8, fast = FALSE, taper=0,
                          log="no", detrend=TRUE, plot=FALSE,
                          kernel("modified.daniell", c(2,2)))


# Plotting periodogram, and smoothed periodograms:
plot(freq, periodogram$spec, type = "l", lwd = 3, col = "black",
     xlab = "Frequency", ylab = "Power Spectrum", 
     main = "Power Spectrum - Temperature [2] ")
lines(freq, fhat.unif$spec, col = "blue", lwd = 3)
lines(freq, fhat.daniel$spec, col = "red", lwd = 3)
lines(freq, ar2.fit, col = "green", lwd = 3)
lines(freq, ar9.fit, col = "orange", lwd = 3)
legend("topright",lty=1,col=c("black","blue","red"),
       c("Periodogram","Uniform weights","Daniell weights"),
       lwd = 2)


# Now I want to find the frequency for which we have 
# maximum value on the periodogram

driving.freq <- freq[which(fhat.daniel$spec
                           == max(fhat.daniel$spec))+ 1] # 0.03157895
driving.freq


############## AR Approach ##########

# Look at PACF
acf(res.Temp8, lag = 31, main = "ACFc")
pacf(res.Temp8, lag = 31, main = "PACF")     # AR(9) ?

# Let's fit several AR(p) for different p
AIC <- rep(0, 20) 
AICc <- rep(0, 20)
BIC <- rep(0, 20)

for(k in 1:20) {
  fit <- ar(res.Temp8, order = k, aic = FALSE)
  sigma2 <- var(fit$resid, na.rm = TRUE)
  BIC[k] <- log(sigma2) + (k*log(T)/T)
  AICc[k] <- log(sigma2) + ((T+k)/(T-k-2))
  AIC[k] <- log(sigma2) + ((T+2*k)/T)
}

IC <- cbind(AIC, BIC)

ts.plot(IC, type = "o", xlab = "p",
        ylab = "AIC and BIC")
text(18.1, -0.48, "BIC", lwd = 2, col = "deeppink3", cex = 1.5)
text(18.1, 0.06, "AIC", lwd = 2, col = "deeppink3", cex = 1.5)
points(2, BIC[2], col = "blue", pch = 19, lwd = 4)
points(9, AIC[9], col = "red", pch = 19, lwd = 4)

legend("bottomright", c("AR(2)", "AR(9)"), pch = 19,
       col = c("blue", "red"))




############ Comparing periodogram with AR spectrum, by choosing AR(2) ##############

p <- 2

# Fitting the time series to an AR(2) and getting the coefficients
fit2 <- arima(res.Temp8, method = "CSS-ML", order = c(p, 0, 0))
coef2 <- fit2$coef[1:p]

# Getting the true spectrum of an AR(2) with that coefficients
ar2.fit <- arma.spec(ar = coef2, n.freq = length(freq) + 1,
                    var.noise = fit2$sigma2)$spec

# Remove the estimate at frequency 0]
ar2.fit <- ar2.fit[-1]

# Plotting periodogram and spectrum AR(2)
plot(freq, periodogram$spec, type = "l", ylab = "Power spectrum", xlab = "Frequency")
lines(freq, ar2.fit, col = "red", lwd = 2)
legend("topright", c("Periodogram", "Spectrum AR(2)"), col = c("black", "red"),
       lty = 1)


############ Comparing periodogram with AR spectrum, by choosing AR(9) ##############

p <- 9

# Fitting the time series to an AR(9) and getting the coefficients
fit9 <- arima(res.Temp8, method = "CSS-ML", order = c(p, 0, 0))
coef9 <- fit9$coef[1:p]

# Getting the true spectrum of an AR(9) with that coefficients
ar9.fit <- arma.spec(ar = coef9, n.freq = length(freq) + 1,
                     var.noise = fit9$sigma2)$spec

# Remove the estimate at frequency 0
ar9.fit <- ar9.fit[-1]

# Plotting periodogram and spectrum AR(9)
plot(freq, periodogram$spec, type = "l", ylab = "Power spectrum", xlab = "Frequency")
lines(freq, ar9.fit, col = "red", lwd = 2)
legend("topright", c("Periodogram", "Spectrum AR(9)"), col = c("black", "red"),
       lty = 1)




# Plotting Periodogram, AR(2), AR(9)
plot(freq, periodogram$spec, type = "l", ylab = "Power spectrum", xlab = "Frequency")
lines(freq, ar2.fit, col = "blue", lwd = 3)
lines(freq, ar9.fit, col = "red", lwd = 3)
legend("topright", c("Periodogram", "Spectrum AR(2)", "Spectrum AR(9)"), col = c("black", "blue", "red"),
       lty = 1, lwd = 3)




freq[which(ar9.fit == max(ar9.fit)) + 1]

# Alright, the driving frequency is 0.04210526 # forst first harmonic









#################### FOURIER REGRESSION ###################


# Functiont to obtain the harmonic, given a specific frequency omega
get_harmonic <- function(x, omega) {
  
  T <- length(x)
  cosine <- cos(2*pi*1:T*omega)
  sine <- sin(2*pi*1:T*omega)
  
  fit <- lm(x ~ 0 + cosine + sine)
  
  U1 <- as.numeric(fit$coefficients[1])
  U2 <- as.numeric(fit$coefficients[2])
  phi <- atan2(U2, U1)
  A <- sqrt(U1^2 + U2^2)
  
  harmonic <- A * cos(2*pi*1:T*omega + phi)
  return(harmonic)
  
}

# Selecting the four frequencies, in order to obtain the four harmonics later
driving.frequencies <- c(freq[which(fhat.daniel$spec
                                    == max(fhat.daniel$spec))+ 1],
                         freq[5 + 1],
                         freq[2 + 1],
                         freq[7 + 1])


# Getting the harmonics and adding the trend
harmonics <- list()

for(i in 1:4) {
  harmonics[[i]] <- get_harmonic(res.Temp8, driving.frequencies[i])
}

# Getting the trend
t <- 1:T
trend.temp8 <- as.vector(fitted(lm(Temp8 ~ t)))

# Final Model
model.temp8 <- harmonics[[1]] + harmonics[[2]] + harmonics[[3]] + harmonics[[4]] + trend.temp8

# Plotting time series, adding single harmonics + final model
plot.ts(Temp8, type = "o", pch = 19, main = "Temperature Patient 8", ylab = "Temperature")
lines(1:T, harmonics[[1]] + trend.temp8, col = "blue", lwd = 3, lty = 3)
lines(1:T, harmonics[[2]] + trend.temp8, col = "chartreuse3", lwd = 3, lty = 3)
lines(1:T, harmonics[[3]] + trend.temp8, col = "grey", lwd = 3, lty = 3)
lines(1:T, harmonics[[4]] + trend.temp8, col = "darkorange1", lwd = 3, lty = 3)
lines(1:T, model.temp8, col = "red", lwd = 5)

legend("bottomright", 
       c(expression(paste(omega, " = 1/24")), 
         expression(paste(omega, " = 1/19")),
         expression(paste(omega, " = 1/47")),
         expression(paste(omega, " = 1/13"))), 
       lty = 3, col = c("blue", "chartreuse3", "grey", "darkorange1"), lwd = 3)
legend("bottomleft", "Total Model", lwd = 5, col = "red", lty = 1)





# 




