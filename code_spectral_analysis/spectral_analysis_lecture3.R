########## LECTURE 3 ###########

## MA spectrum
ma.spec = function(theta,sigma2=1,num.freq=1000){ 
  freq = seq(0,0.5-0.5/num.freq,length.out=num.freq)
  fw = sigma2*(1+theta^2+2*theta*cos(2*pi*freq))
  return(fw)
}

## AR spectrum
ar1.spec = function(phi,sigma2=1,num.freq=1000){ 
  freq = seq(0,0.5-0.5/num.freq,length.out=num.freq)
  fw = sigma2/(1+phi^2-2*phi*cos(2*pi*freq))
  return(fw)
}


# AR1, phi = 0.5
T <- 128
xt <- arima.sim(n = T, model = list(ar = c(0.5)))
freq <- 0:(T/2-1)/T
plot.ts(xt)

I <- abs(fft(xt))^2/T
plot(freq, I[1:(T/2)], type = "l", ylab = "Power spectrum", xlab = "Frequency")
lines(freq, ar1.spec(phi = 0.5, num.freq = T/2), col = "red", lwd = 2)


# AR1, phi = -0.5
T <- 128
xt <- arima.sim(n = T, model = list(ar = c(-0.5)))
freq <- 0:(T/2-1)/T
plot.ts(xt)

I <- abs(fft(xt))^2/T
plot(freq, I[1:(T/2)], type = "l", ylab = "Power spectrum", xlab = "Frequency")
lines(freq, ar1.spec(phi = -0.5, num.freq = T/2), col = "red", lwd = 2)


# Averaging (?)
MC <- 1000
T <- 8000
period.mc <- matrix(0, nrow = 1000, ncol = T)
for(mc in 1:MC) {
  period.mc[mc, ] <- abs(fft(arima.sim(n = T, model = list(ar = c(-0.5)), sd = 1)))^2/T
}
freq <- 0:(T/2 - 1)/T

I.mc <- colMeans(period.mc)
plot(freq, I.mc[1:(T/2)], type = "l", col = "red",
     lwd = 2)
lines(freq, I.mc[1:(T/2)], col = "black", lwd = 2)


# Now we put different T
Tstar <- 512
freq <- 0:(Tstar/2 - 1)/(Tstar)
xt <- arima.sim(Tstar, model = list(ar = c(0.5)), sd = 1)
plot.ts(xt)

I = abs(fft(xt))^2/Tstar
plot(freq, I[1:(Tstar/2)], type = "l")
lines(freq, ar1.spec(phi = 0.5, num.freq = Tstar/2), col = "red", lwd = 2)


Tstar <- 2048
freq <- 0:(Tstar/2 - 1)/(Tstar)
xt <- arima.sim(Tstar, model = list(ar = c(0.5)), sd = 1)
plot.ts(xt)

I = abs(fft(xt))^2/Tstar
plot(freq, I[1:(Tstar/2)], type = "l")
lines(freq, ar1.spec(phi = 0.5, num.freq = Tstar/2), col = "red", lwd = 2)


Tstar <- 8192
freq <- 0:(Tstar/2 - 1)/(Tstar)
xt <- arima.sim(Tstar, model = list(ar = c(0.5)), sd = 1)
plot.ts(xt)

I = abs(fft(xt))^2/Tstar
plot(freq, I[1:(Tstar/2)], type = "l")
lines(freq, ar1.spec(phi = 0.5, num.freq = Tstar/2), col = "red", lwd = 2)


### So the periodogram it's not a consistent estimator at all!!!



########## SMOOTHED PERIODOGRAM (UNIFORMS WEIGHTS)

T <- 512
freq <- seq(0, 0.5, len = T/2)
xt <- arima.sim(n = T, model = list(ar = c(0.5)))
I <- abs(fft(xt))^2/T
plot(freq, I[1:(T/2)], lwd = 1, col="black", type="l",
     main = "No Smoothing", ylab = "Power Spectrum", xlab = "Frequency")
lines(freq,ar1.spec(0.5,sigma2=1,num.freq=T/2),col="red",lwd = 3)


M <- 11
# Creating kernel coefficients with equal weights given by 1/(2M + 1)
k <- kernel(rep(1, M + 1)/(2*M + 1))

plot(freq, ar1.spec(0.5, sigma2 = 1, num.freq = T/2),
     lwd = 2, main = "T = 512, 2M + 1 = 11",  
     ylab = "Power Spectrum", xlab = "Frequency", 
     col = "red", type = "l")
spec.pgram(xt, kernel = k, add = TRUE, lwd = 2)


M <- 31
# Creating kernel coefficients with equal weights given by 1/(2M + 1)
k <- kernel(rep(1, M + 1)/(2*M + 1))

plot(freq, ar1.spec(0.5, sigma2 = 1, num.freq = T/2),
     lwd = 2, main = "T = 512, 2M + 1 = 31",  
     ylab = "Power Spectrum", xlab = "Frequency", 
     col = "red", type = "l")
spec.pgram(xt, kernel = k, add = TRUE, lwd = 2)


M <- 75
# Creating kernel coefficients with equal weights given by 1/(2M + 1)
k <- kernel(rep(1, M + 1)/(2*M + 1))

plot(freq, ar1.spec(0.5, sigma2 = 1, num.freq = T/2),
     lwd = 2, main = "T = 512, 2M + 1 = 75",  
     ylab = "Power Spectrum", xlab = "Frequency", 
     col = "red", type = "l")
spec.pgram(xt, kernel = k, add = TRUE, lwd = 2)


# It looks that by using M = 11 we get the best fit. 
