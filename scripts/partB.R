# Part B

# b.i.

N <- 6.022e23
mDiox <- 44.1 / N / 1000 # g/mol to kg
mHyd <- 1.008 / N / 1000 # g/mol to kg
mMeOH <- 32.04 / N / 1000 # g/mol to kg
mWat <- 18.015 / N / 1000 # g/mol to kg

vibDiox <- data$vibration[1] # K
vibHyd <- data$vibration[2] # K
vibMeOH <- data$vibration[3] # K
vibWat <- data$vibration[4] # K

iDiox <- data$inertia[1] * 1e-7 # g*cm^2 to kg*m^2
iHyd <- data$inertia[2] * 1e-7 # g*cm^2 to kg*m^2
iMeOH <- data$inertia[3] * 1e-21 # g^3*cm^6 to kg^3*m^6
iWat <- data$inertia[4] * 1e-21 # g^3*cm^6 to kg^3*m^6

sigmaDiox <- data$sigma[1]
sigmaHyd <- data$sigma[2]
sigmaMeOH <- data$sigma[3]
sigmaWat <- data$sigma[4]

# q translational

qTransRxn1 <- function(temp) {
  (qTrans(temp, mMeOH) * qTrans(temp, mWat)) / 
    (qTrans(temp, mDiox) * qTrans(temp, mHyd)^3)
}

# q vibrational

qVibRxn1 <- function(temp) {
  (qVib(temp, vibMeOH) * qVib(temp, vibWat)) / 
    (qVib(temp, vibDiox) * qVib(temp, vibHyd)^3)
}

#q nonlinear rotation

qRotRxn1 <- function(temp) {
  (qRotNonlin(temp, iMeOH, sigmaMeOH) * qRotNonlin(temp, iWat, sigmaWat)) / 
    (qRotlin(temp, iDiox, sigmaDiox) * qRotlin(temp, iHyd, sigmaHyd)^3)
}

# q electronic
# theta_electronic is large and T is relativly small so qElec = g0
# g0 = 1 for all the molecules we are dealing with

qEleRxn1 <- 1

# Kp for Reaction 1

R <- 0.00814 # kJ/mol*K
k2 <- 8.206e-5 / N # R/N Units are m^3 atm / K
KpRxn1 <- function(temp) {
  (k2 * temp)^-2 * qTransRxn1(temp) * qVibRxn1(temp) * qRotRxn1(temp) *
    qEleRxn1 * exp(deltaDRxn1 / (R * temp))
}

# b.ii.

tempC <- seq(from = 50, to = 1500, by = 1)
tempK <- tempC + 273
inverseTempK <- 1 / tempK

KpReaction1 <- sapply(tempK, KpRxn1)
log10KpRxn1 <- log10(KpReaction1)
lnKpRxn1 <- log(KpReaction1)

ReactionOne <- data.frame(
  "temp" = tempC,
  "inverseTemp" = inverseTempK,
  "Kp" = KpReaction1,
  "log10Kp" = log10KpRxn1,
  "lnKp" = lnKpRxn1
)

KpPlotRxn1 <- ggplot(data = ReactionOne,
                 mapping = aes(
                   x = temp,
                   y = log10Kp)) +
  geom_line() +
  labs(title = "Figure 1",
       subtitle = "Reaction One Equilibrium Constant") +
  xlab("Temperature (Celcius)") +
  ylab("Kp (atm^-2) Log10 Scale")

ggsave("pics/Figure1.pdf")

# b.iii.

HoffPlotRxn1 <- ggplot(data = ReactionOne,
                   mapping = aes(
                     x = inverseTemp,
                     y = lnKp)) +
  geom_line() +
  labs(title = "Figure 2",
       subtitle = "Reaction One van't Hoff Plot") +
  xlab("Temperature (Celcius^-1)") +
  ylab("Kp (atm^-2) Natural Log Scale") +
  labs(caption = "Figure")

ggsave("pics/Figure2.pdf")

deltaHRxn1 <- -R * (ReactionOne$lnKp[500] - ReactionOne$lnKp[600]) /
  (ReactionOne$inverseTemp[500] - ReactionOne$inverseTemp[600]) # kJ/mol

# b.iv.

deltaHDiox <- data$heatFormation[1]
deltaHHyd <- data$heatFormation[2]
deltaHMeOH <- data$heatFormation[3]
deltaHWat <- data$heatFormation[4]

deltaHPrimeRxn1 <- deltaHMeOH + deltaHWat - deltaHDiox - 3 * deltaHHyd