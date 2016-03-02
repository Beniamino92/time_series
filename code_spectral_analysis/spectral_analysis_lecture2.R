##### LECTURE 2 ######

T <- 256
signal <- 2*cos(2*pi*1:T/50 + .6*pi)
noise <- rnorm(T, 0, 1)
x <- signal + noise

par(mfrow = c(3, 1))
plot.ts(signal, main = "signal")
plot.ts(noise, main = "noise")
plot.ts(x, main = "signal + noise")

cosine <- cos(2*pi*1:T/50)
sine <- sin(2*pi*1:T/50)
fit <- lm(x ~ 0 + cosine + sine)
summary(fit)

par(mfrow = c(1, 1))
plot.ts(x, lty = 2, col = "grey")
lines(fitted(fit), col = "red", lwd = 2)

I <- abs(fft(x))^2/500 # the periodogram
P <- (4/500) * I[1:250] # the scaled periodogram
f <- 0:249/500 # frequencies
plot(f, P, type = "l", xlab = "frequency", ylab = "periodogram")
length(f)

I.mod <- abs(fft(x))^2
P.mod <- I.mod[1:250]
f.mod <- 0:249
plot(f.mod, P.mod, type = "l")

f[which(P ==max(P))]
length(P)


### Example
T <- 100

x1 <- 2*cos(2*pi*1:T*6/T) + 3*sin(2*pi*1:T*6/T)
x2 <- 4*cos(2*pi*1:T*10/T) + 5*sin(2*pi*1:T*10/T)
x3 <- 6*cos(2*pi*1:T*40/T) + 7*sin(2*pi*1:T*40/T)

par(mfrow = c(3, 1))
plot.ts(x1)
plot.ts(x2)
plot.ts(x3)

x <- x1 + x2 + x3

par(mfrow = c(1, 1))
plot.ts(x)

I <- (abs(fft(x))^2)/T # the periodogram
f <- 0:(T/2 - 1)/T


par(mfrow = c(1, 1))
plot(f, I[1:(T/2)], type = "l",
     xlab = "Frequency", ylab = "Periodogram")



## MA Spectrum
ma.spec <- function(theta, sigma2 = 1, num.freq = 1000) {
  freq <- seq(from = 0, to = 0.5 - 0.5/num.freq, len = num.freq)
  fw <- sigma2*(1 + theta^2 + 2*theta*cos(2*pi*freq))
  return(fw)
}

## AR Spectrum
ar1.spec <- function(phi, sigma2 = 1, num.freq = 1000) {
  freq <- seq(from = 0, to = 0.5 - 0.5/num.freq, length.out = num.freq)
  fw <- sigma2/(1 + phi^2 - 2*phi*cos(2*pi*freq))
  return(fw)
}



# Example

T <- 128
freq <- seq(0, 0.5, len = T/2)

# Plot MA(1) spectrum, theta = 0.5
plot(freq, ma.spec(0.5, num.freq = T/2), type = "l")
# Plot MA(1) spectrum, theta = -0.5
plot(freq, ma.spec(-0.5, num.freq = T/2), type = "l")


par(mfrow = c(2, 1))
x <- arima.sim(n = T, list(ar = c(0.9)))
plot.ts(x, main = expression(phi == 0.9), ylab = expression(X[t]))
x <- arima.sim(n = T, list(ar = c(-0.9)))
plot.ts(x, main = expression(phi == -0.9), ylab = expression(X[t]))


par(mfrow = c(1, 1))
plot(freq, ar1.spec(0.9, num.freq = T/2), type = "l")
plot(freq, ar1.spec(-0.9, num.freq = T/2), type = "l")




# PERIODOGRAM AND SPECTRUM FOR AR AND MA MODELS

T <- 128
phi <- 0.5
x <- arima.sim(n = T, list(ar = c(phi)))
plot.ts(x, main = expression(phi == 0.5), ylab = expression(X[t]))

I <- abs(fft(x))^2/T
fx <- ar1.spec(phi, sigma2 = 1, num.freq = T/2)
freq <- 0:(T/2 - 1)/T
plot(freq, I[1:(T/2)], type="l", xlab="Frequency", ylab="Periodogram")
lines(freq,fx,col="red",lwd=2)

phi <- -0.5
x <- arima.sim(n = T, list(ar = c(phi)))
plot.ts(x, main = expression(phi == -0.5), ylab = expression(X[t]))

I <- abs(fft(x))^2/T
fx <- ar1.spec(phi, sigma2 = 1, num.freq = T/2)
plot(freq, I[1:(T/2)], type = "l", xlab = "Frequency", ylab = "Periodogram")
lines(freq, fx, col = "red", lwd = 2)

