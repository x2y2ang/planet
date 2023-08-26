
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# library
library(vegan)
library(FactoMineR)
library(lavaan)
library(data.table)
library(MASS)

all.plant = read.csv("../2.DX2013_145sites/data/plant1314.cov.network.mass.rich.for.SEM.csv")
colnames(all.plant)

fit <- lm(Veg.mass ~ Plant.net + pre.AI + Plant.rich,
          data=all.plant)

summary(fit)

plot(all.plant$Plant.rich, all.plant$Veg.mass)
plot(all.plant$Plant.net, all.plant$Veg.mass)


# aridity <= 0.69
all.below = subset(all.plant, Elevation > 4700)

fit.low <- lm(Veg.mass ~ Plant.net + pre.AI + Plant.rich,
          data=all.below)

summary(fit.low)

plot(all.below$Plant.rich, all.below$Veg.mass)
plot(all.below$Plant.net, all.below$Veg.mass)


# aridity >= 0.69
all.high = subset(all.plant, Elevation <= 4700)

fit.high <- lm(Veg.mass ~ Plant.net + pre.AI + Plant.rich,
              data=all.high)

summary(fit.high)

plot(all.high$Plant.rich, all.high$Veg.mass)
plot(all.high$Plant.net, all.high$Veg.mass)

a = all.high[order(all.high$Plant.net, decreasing = F), ]




0.000065
0.00001417
0.0000002
