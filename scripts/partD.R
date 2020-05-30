# Part D

# d.i.

# calculate ratios based on a feed on 100:300 PCO2:PH2

conversion <- conversion %>%
  mutate("Diox" = 100 - conversion) %>%
  mutate("Hyd" = 300 - conversion) %>%
  mutate("MeOH" = conversion) %>%
  mutate("Wat" = conversion) %>%
  mutate("Total" = Diox + Hyd + MeOH + Wat) %>%
  mutate("pDiox" = Diox / Total) %>%
  mutate("pHyd" = Hyd / Total) %>%
  mutate("pMeOH" = MeOH / Total) %>%
  mutate("pWat" = Wat / Total) %>%
  mutate("Kp" = (pMeOH * pWat) / (pDiox * pHyd^3)) %>%
  mutate("log10Kp" = log10(Kp)) %>%
  mutate("lnKp" = log(Kp))

KpPlot <- ggplot(data = conversion) +
  geom_point(mapping = aes(x = temp, y = log10Kp, color = pressure)) +
  geom_point(data = ReactionOne,
             mapping = aes(x = temp, y = log10Kp)) +
  ggtitle("Equilibrium Constant (Pressure) for Reaction One from Conversion Data") +
  xlab("Temperature (Celcius)") +
  ylab("Kp (atm^-2) Log10 Scale")

# d.i.

k3 <- (8.206e-2) / (6.022e23)

conversion <- conversion %>%
  mutate("K" = Kp * (k3 * (temp+273))^2) %>%
  mutate("log10K" = log10(K)) %>%
  mutate("lnK" = log(K))

KPlot <- ggplot(data = conversion) +
  geom_point(mapping = aes(x = temp, y = log10K, color = pressure)) +
  ggtitle("Equilibrium Constant for Reaction One from Conversion Data") +
  xlab("Temperature (Celcius)") +
  ylab("K (atm^-2) Log10 Scale")

