library(dplyr)
library(ggplot2)
library(plotly)

data <- read.csv("data/table1.csv")
bondEnergy <- read.csv("data/table2.csv")
gasLaw <- read.csv("data/table3.csv")
conversion <- read.csv("data/table4.csv")

# Problem 1

# a. (Units are kJ/mol)

source("scripts/partA.R")

deltaDRxn1 <- disMeOH + disWat - disDiox - 3 * disHyd # kJ/mol
print(deltaDRxn1)

# b.i.

source("scripts/partitions.R")
source("scripts/partB.R")

# Kp for Reaction 1

R <- 0.00814 # kJ/mol
k2 <- (8.206e-5) / (6.022e23) # R/N Units are m^3 atm / K
KpRxn1 <- function(temp) {
  (k2 * temp)^-2 * qTransRxn1(temp) * qVibRxn1(temp) * qRotRxn1(temp) *
    qEleRxn1 * exp(deltaDRxn1 / (R * temp))
}

# b.ii.

ggplotly(KpPlotRxn1)


# b.iii.

ggplotly(HoffPlotRxn1)

print(deltaHRxn1)
  
# b.iV.

print(deltaHPrimeRxn1)

# c.i.

source("scripts/partC.R")

print(deltaDRxn2)

# c.ii.

# Kp for Reaction 2

R <- 0.00814 # kJ/mol
k2 <- (8.206e-5) / (6.022e23) # R/N Units are m^3 atm / K
KpRxn1 <- function(T) {
  (k2 * T)^-2 * qTransRxn1(T) * qVibRxn1(T) * qRotRxn1(T) *
    qEleRxn1 * exp(deltaDRxn1 / (R * T))
}

# c.iii.

ggplotly(KpPlotRxn2)

# c.iv.

ggplotly(HoffPlotRxn2)

print(deltaHRxn2)

# c.v.

print(deltaHPrimeRxn2)

# c.vi.

ggplotly(hPlotRxn2)

ggplotly(uPlotRxn2)

ggplotly(sPlotRxn2)

# c.vii.

# At equilibrium deltaU is equal to zero. That point is at 510 degrees C.

# d.i

source("scripts/partD.R")

ggplotly(KpPlot)

# d.ii.

ggplotly(KPlot)

# d.iv.

ggplotly(lnKPlot)

print(dlnKdP)






