# Question 1 Part 7
# Function to generate an AR process
ar2_sim <- function(r1, r2, len = 200) {
  out <- rep(0, len)
  a1 <- (1 / r1) + (1 / r2)
  a2 <- - 1 / (r1 * r2)
  epsilon <- rnorm(len)
  for(t in 3:len) {
    out[t] <- out[t - 2] * a1 + out[t - 1] * a2 + epsilon[t]
  }
  out
}

ar2_sim(2.5, 2.5)


acf(ar2_sim(5, 2))
pacf(ar2_sim(5, 2))

emp_acf <- function(r1, r2, h) {
  R1 <- 1 / r1
  R2 <- 1 / r2
  rho <- ((R1 / (1 - R1 ^ 2) - R2 / (1 - R2 ^ 2)) ^ - 1) * ((R1 ^ (h + 1)) / (1 - R1 ^ 2) - ((R2 ^ (h + 1)) / (1 - R2 ^ 2)))
  rho
}

h <- seq(1, 30)


r1 <- 5; r2 <- 2
acf(ar2_sim(r1, r2, len = 1000))
pacf(ar2_sim(r1, r2))
true_acf <- sapply(h, function(x) Re(emp_acf(r1, r2, x)))
plot(h, true_acf, type = "h"); abline(h = 0)



r1 <- 10; r2 <- 2
acf(ar2_sim(r1, r2, len = 1000))
pacf(ar2_sim(r1, r2))
true_acf <- sapply(h, function(x) Re(emp_acf(r1, r2, x)))
plot(h, true_acf, type = "h"); abline(h = 0)


r1 <- - 10; r2 <- 2
acf(ar2_sim(r1, r2, len = 1000))
pacf(ar2_sim(r1, r2))
true_acf <- sapply(h, function(x) Re(emp_acf(r1, r2, x)))
plot(h, true_acf, type = "h"); abline(h = 0)

r1 <- 2 /3 * (1 + sqrt(3) * 1i); r2 <- 2 / 3 * (1 - sqrt(3) * 1i)
acf(Re(ar2_sim(r1, r2, len = 1000)))
pacf(Re(ar2_sim(r1, r2)))
true_acf <- sapply(h, function(x) Re(emp_acf(r1, r2, x)))
plot(h, true_acf, type = "h"); abline(h = 0)


r1 <- (1 - 0025 + 0.3258 * 1i); r2 <- (1 - 0025 - 0.3258 * 1i)
acf(Re(ar2_sim(r1, r2, len = 1000)))
pacf(Re(ar2_sim(r1, r2)))
true_acf <- sapply(h, function(x) Re(emp_acf(r1, r2, x)))
plot(h, true_acf, type = "h"); abline(h = 0)
