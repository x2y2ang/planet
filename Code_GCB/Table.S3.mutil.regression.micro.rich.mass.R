
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


all = read.csv("../data/GCB/micro.rich.inter.biomass.combine.for.SEM.csv", row.names = 1)
colnames(all)

fit <- lm(Micro.mass ~ Micro.net + pre.AI + Micro.rich,
          data=all)

summary(fit)

plot(all$Micro.rich, all$Micro.mass)
plot(all$Micro.net, all$Micro.mass)


# aridity < 0.78

all.below = subset(all, pre.AI < 0.78)


fit.low <- lm(Micro.mass ~ Micro.net + pre.AI + Micro.rich,
              data=all.below)

summary(fit.low)

plot(all.below$Micro.rich, all.below$Micro.mass)
plot(all.below$Micro.net, all.below$Micro.mass)



# aridity >= 0.78
all.high = subset(all, pre.AI >= 0.78)

fit.high <- lm(Micro.mass ~ Micro.net + pre.AI + Micro.rich,
               data=all.high)

summary(fit.high)

plot(all.high$Micro.rich, all.high$Micro.mass)
plot(all.high$Micro.net, all.high$Micro.mass)



