# q Translational (units are m^-3)

h <- 6.63e-34 # Js
k1 <- 1.38e-23 # J/K

qTrans <- function(temp, m) {
  (2*pi*m*k1*temp / (h^2))^(3/2)
}

# q vibrational (unitless)

qVib <- function(temp, theta) {
  1 / (1 - exp(-theta / temp))
}

# q linear rotation

qRotlin <- function(temp, I, sigma) {
  (8*pi^2*I*k1*temp) / (sigma*h^2)
}

# q nonlinear rotation

qRotNonlin <- function(temp, I, sigma) {
  (pi*I)^(1/2) * (8*pi^2*k1*temp/(h^2))^(3/2) / sigma
}

# q electronic
# theta_electronic is large and T is relativly small so qElec = g0
# g0 = 1 for all the molecules we are dealing with