################# Spectral Analysis Temperature ##################
#################       Patient 8               ##################

setwd("C:/Users/Beniamino/Desktop/Project_8")
source("spectral_analysis_healthy_patients.R")

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
legend("topright",lty=1,col=c("black","blue","red"),
       c("Raw Periodogram","Uniform weights","Daniell weights"),
       lwd = 2)


# Now I want to find the frequency for which we have 
# maximum value on the periodogram

driving.freq <- freq[which(fhat.daniel$spec
                           == max(fhat.daniel$spec))] # 0.03157895



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



