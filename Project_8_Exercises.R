##### Project 8, Interpreting the Autocorrelation Function


########### P1: Function to evaluate ACF ##########

my.acf <- function(y, lag) {
  
  N <- length(y)
  y.mean <- mean(y)
  
  out <- numeric(lag + 1)
  out[1] <- var(y)
  
  # Evaluating ACF for each lag
  for(h in 2:(lag + 1)) {
    temp <- 0.0
    for(t in 1:(N-(h-1))) {
      temp <- temp + (y[t] - y.mean)*(y[t+(h-1)] - y.mean)
    }
    out[h] <- temp/N
  }
   
  # Plotting ACF and time series
  par(mfrow = c(2, 1))
  grid <- seq(from = 0, to = lag, len = length(out))
  plot.ts(y, type = "l", ylab = expression(Y[t]))
  plot(grid, out, type = "h", main = "ACF", xlab = "Lag", ylab = "ACF")
  abline(h = 0)
  par(mfrow = c(1, 1))
  
}

# From now on I will use the acf available in R

############# P2: Take Y_t = e_t (Gaussian White Noise) ##########

T <- 200
Y <- rnorm(T, 0, 1)
acf(Y)



########### P3: Introducing outliers ##############

T <- 200
Y.outliers <- rnorm(T, 0, 1)

# Percentege of outliers
n.outliers <- 20
outliers.indexes <- sample(T, n.outliers)

# Simulating outliers from normal (0, 10)
Y.outliers[outliers.indexes] <- replicate(n.outliers, rnorm(1, 0, 20))

par(mfrow = c(2, 1))  
plot(Y)
plot(Y.outliers)
par(mfrow = c(1, 1))

#ACF
acf(Y)
acf(Y.outliers)


########## P4: Random Walk ##########

T <- 200
random.walk <- numeric(T)
errors <- rnorm(T, 0, 1)

random.walk <- errors[1]
for(t in 2:T) {
  random.walk[t] <- random.walk[t-1] + errors[t]
}

#ACF
acf(random.walk)
pacf(random.walk)




########## P5: AR1 ############

T <- 200

ar1 <- function(alpha, T = 200) {
  errors <- rnorm(T)
  Y <- numeric(T)
  Y[1] <- rnorm(1, 0, sd = sqrt((1/(1-alpha^2))))
  for(t in 2:T) {
    Y[t] <- alpha * Y[t-1] + errors[t]
  }
  Y
}

# Testing behaviour for different alpha.

# Positive
acf(ar1(0.01))
pacf(ar1(0.01))

acf(ar1(0.10))
pacf(ar1(0.10))

acf(ar1(0.50))
pacf(ar1(0.50))

acf(ar1(0.75))
pacf(ar1(0.75))

acf(ar1(0.90))
pacf(ar1(0.90))

# Negative
acf(ar1(-0.01))
pacf(ar1(-0.01))

acf(ar1(-0.10))
pacf(ar1(-0.10))

acf(ar1(-0.50))
pacf(ar1(0.50))

acf(ar1(-0.75))
pacf(ar1(0.75))

acf(ar1(-0.90))
pacf(ar1(-0.90))


# True autocovariance function of ACF
acf.ar1 <- function(alpha, h) {
  (alpha^h)/(1-alpha^2)
}

# Testing
my.ar1 <- acf(ar1(alpha = -0.8, T = 1e5), type = "covariance")$acf
plot(1:20, acf.ar1(alpha = -0.8, 1:20), type = "l", col = "red", ylim = c(-5, 5),
    lwd = 2)
lines(1:20, my.ar1[1:20], type = "h")


##### P6: Deterministic Logistic Process:
T <- 200
Y <- numeric(T)
Y[1] <- runif(1)
lambda <- 1

for(t in 2:T) {
  Y[t] <- lambda*Y[t-1]*(1 - Y[t-1])
}

acf(Y)
pacf(Y)


##
