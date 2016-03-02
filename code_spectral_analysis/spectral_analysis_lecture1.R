################# Spectral Analysis of Time Series ############

###############
## LECTURE 1 ##
###############

# White noise
xt <- rnorm(1024)
gamma <- acf(xt, type = "covariance")
rho <- acf(xt, type = "correlation")
plot.ts(xt, main = expression(sigma^2==1), ylab = "X(t)")

plot(gamma, main = "Sample Autocovariance Function", ylab = expression(gamma(h)))
plot(rho,main = "Sample Autocorrelation Function")

# X_{t} = Z_{t}
plot.ts(rnorm(256), main = expression(sigma^2==1))
# MA(1)  X_{t} = Z_{t} + theta_{1} * Z_{t-1}
plot.ts(arima.sim(n = 256, list(ma = c(0.5)), sd = 1),
                  main = as.expression(expression(theta==0.5)))
plot.ts(arima.sim(n = 256, list(ma = c(-0.5)), sd = 1),
        main = as.expression(expression(theta==-0.5)))


# MA(1)
xt <- arima.sim(n = 256, list(ma = c(0.5)), sd = 1)
plot.ts(xt) 
gamma <- acf(xt, type = "covariance", main = expression(theta==0.5),
             ylab = expression(gamma(h)))
rho <- acf(xt, type = "correlation")
plot(gamma,main = "Sample Autocovariance Function")
plot(rho,main = "Sample Autocorrelation Function")



