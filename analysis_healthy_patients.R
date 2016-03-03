# Getting the data
setwd("C:/Users/Beniamino/Desktop/Data_Healthy_Patients")

# Rest activity
RA2.data <- read.csv("RA2.csv", header = T)         # 172 obs
RA8.data <- read.csv("RA8.csv", header = T)         # 97 obs
RA24.data <- read.csv("RA24.csv", header = T)       # 98 obs
RA26.data <- read.csv("RA26.csv", header = T)       # 97 obs

# Temperature
Temp2.data <- read.csv("Temp2.csv", header = T)     # 128 obs
Temp8.data <- read.csv("Temp8.csv", header = T)     # 97 obs
Temp24.data <- read.csv("Temp24.csv", header = T)   # 97 obs
Temp26.data <- read.csv("Temp26.csv", header = T)   # 95 obs



# If I a want to compare different patients, I need to be careful
# regarding the starting time when they start the measurments. Especially
# for coherence analysis. 


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


# Plotting Rest Activity
plot.ts(RA2, type = "l", lwd = 2, col = "blue", ylab = "Rest Activity",
        main = "Rest Activity [2]")
plot.ts(RA8, type = "l", lwd = 2, col = "blue", ylab = "Rest Activity",
        main = "Rest Activity [8]")
plot.ts(RA24, type = "l", lwd = 2, col = "blue", ylab = "Rest Activity",
        main = "Rest Activity [24]")
plot.ts(RA26, type = "l", lwd = 2, col = "blue", ylab = "Rest Activity", 
        main = "Rest Activity [26]")

# Plotting Temperature
plot.ts(Temp2, type = "l", lwd = 2, col = "blue", ylab = "Temperature",
        main = "Temperature [2]")
plot.ts(Temp8, type = "l", lwd = 2, col = "blue", ylab = "Temperature",
        main = "Temperature [8]")
plot.ts(Temp24, type = "l", lwd = 2, col = "blue", ylab = "Temperature",
        main = "Temperature [24]")
plot.ts(Temp26, type = "l", lwd = 2, col = "blue", ylab = "Temperature", 
        main = "Temperature [26]")










