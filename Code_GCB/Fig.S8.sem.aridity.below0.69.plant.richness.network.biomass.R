rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# library
library(vegan)
library(FactoMineR)
library(lavaan)
library(semPlot)
library(AICcmodavg)
library(data.table)
library(MASS)

all.plant = read.csv("../data/GCB/plant1314.cov.network.mass.rich.for.SEM.mean.csv",
                     row.names = 1)


# aridity < 0.69

all.below = subset(all.plant, pre.AI < 0.69)


net.below.data = all.below[, c(2:5)]

net.below.SEM.md = na.omit(net.below.data)
net.below.SEM = scale(net.below.SEM.md,center=T,scale=T) 
net.below.SEM = as.data.frame(net.below.SEM)
net.below.SEM = data.table(net.below.SEM)
colnames(net.below.SEM)


lvmod.1 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Plant.net ~ pre.AI + Plant.rich
 Veg.mass ~ pre.AI + Plant.rich + Plant.net
'
lvmod.1.fit <- sem(lvmod.1, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Plant.net ~ pre.AI + Plant.rich
 Veg.mass ~ Plant.rich + Plant.net # pre.AI +
'
lvmod.2.fit <- sem(lvmod.2, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


source("../function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.1.fit,lvmod.2.fit),
              c("model1","model2")) 


semPaths(lvmod.2.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)


summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.2.fit, c("cfi", "rmsea", "srmr"))



