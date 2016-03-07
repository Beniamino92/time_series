######## PATIENT DATA & FUNCTIONS ############

# setwd("C:/Users/Beniamino/Desktop/Data_Healthy_Patients")
setwd("/homes/hadjamar/Documents/Project_8/time_series")
setwd("/homes/hadjamar/Documents/Project_8/")
# setwd("/home/hadjamar/Desktop/Project_8")


# Rest activity
RA2.data <- read.csv("RA2.csv", header = T)         # 172 obs
RA8.data <- read.csv("RA8.csv", header = T)         # 97 obs
RA24.data <- read.csv("RA24.csv", header = T)       # 98 obs
RA26.data <- read.csv("RA26.csv", header = T)       # 97 obs

# Temperature
Temp2.data <- read.csv("Temp2.csv", header = T)     # 128 obs
Temp8.data <- read.csv("Temp8.csv", header = T)     # 95 obs
Temp24.data <- read.csv("Temp24.csv", header = T)   # 97 obs
Temp26.data <- read.csv("Temp26.csv", header = T)   # 95 obs

# Getting Time Series

RA2 <- RA2.data[, 1]
RA8 <- RA8.data[, 1]
RA24 <- RA24.data[, 1]
RA26 <- RA26.data[, 1]

Temp2 <- Temp2.data[, 1]
Temp8 <- Temp8.data[, 1]
Temp24 <- Temp24.data[, 1]
Temp26 <- Temp26.data[, 1]

# Getting Measuraments Time

time.RA2 <- RA2.data[, 2]
time.RA8 <- RA8.data[, 2]
time.RA24 <- RA24.data[, 2]
time.RA26 <- RA26.data[, 2]

time.Temp2 <- Temp2.data[, 2]
time.Temp8 <- Temp8.data[, 2]
time.Temp24 <- Temp24.data[, 2]
time.Temp26 <- Temp26.data[, 2]


##### FUNCTIONS:


# Function to get and/or plot of -->
#    [periodogram, smoothed periodogram (uniform weights),
#     smoothed periodogram (daniell weights)]

make_periodograms <- function(residuals.data, freq, plot = TRUE, M = 2) {
  
  # Periodogram 
  periodogram <- spec.pgram(residuals.data, fast = FALSE, taper=0,
                            log="no", detrend=TRUE, plot=FALSE)
  
  # Smoothed Periodogram (Uniform weights)
  fhat.unif <- spec.pgram(residuals.data, fast = FALSE, taper=0,
                          log="no", detrend=TRUE, plot=FALSE,
                          kernel = kernel(rep(1, M +1)/(2*M + 1)))
  
  # Smoothed Periodogram (Daniell modified weights)
  fhat.daniel <- spec.pgram(residuals.data, fast = FALSE, taper=0,
                            log="no", detrend=TRUE, plot=FALSE,
                            kernel("modified.daniell", c(2,2)))
  
  # Plotting periodogram, and smoothed periodograms:
  if( plot == TRUE) {
    par(mfrow = c(1, 1))
    plot(freq, periodogram$spec, type = "l", lwd = 3, col = "black",
         xlab = "Frequency", ylab = "Power Spectrum")
    lines(freq, fhat.unif$spec, col = "blue", lwd = 3)
    lines(freq, fhat.daniel$spec, col = "red", lwd = 3)
    legend("topright",lty=1,col=c("black","blue","red"),
           c("Periodogram","Uniform weights","Daniell weights"),
           lwd = 2)
  }
  
  return(list(periodogram = periodogram$spec,
              fhat.unif = fhat.unif$spec,
              fhat.daniel = fhat.daniel$spec))
}


# Function to obtain the harmonic, given a specific frequency omega

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



# Function to obtain AIC and BIC for the fitting to an
# AR model of order p, for i in 1:p.

performance.ar.p <- function(residuals.data, p, plot = TRUE) {
  
  T <- length(residuals.data)
  
  AIC <- rep(0, p) 
  BIC <- rep(0, p)
  
  for(k in 1:p) {
    fit <- ar(residuals.data, order = k, aic = FALSE)
    sigma2 <- var(fit$resid, na.rm = TRUE)
    BIC[k] <- log(sigma2) + (k*log(T)/T)
    AIC[k] <- log(sigma2) + ((T+2*k)/T)
  }
  
  if (plot == TRUE) {
    IC <- cbind(AIC, BIC)
    ts.plot(IC, type = "o", xlab = "p",
            ylab = "AIC and BIC", col = c("red", "blue"),
            lwd = 2)
    points(which(BIC == min(BIC)), BIC[which(BIC == min(BIC))], 
           col = "black", pch = 19, lwd = 4)
    points(which(AIC == min(AIC)), AIC[which(AIC == min(AIC))], 
           col = "black", pch = 19, lwd = 4)
    legend("bottomright", c("AIC", "BIC"), lty = 1, lwd = 2,
           col = c("red", "blue"))
  }
  
  return(list(AIC = AIC, BIC = BIC))
}






