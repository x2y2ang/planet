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


all = read.csv("micro.rich.inter.biomass.combine.csv", row.names = 1)

# aridity < 0.78

all.below = subset(all, Elevation >= 4700)


net.below.data = all.below

net.below.SEM.md = na.omit(net.below.data)
net.below.SEM = scale(net.below.SEM.md,center=T,scale=T) 
net.below.SEM = as.data.frame(net.below.SEM)
net.below.SEM = data.table(net.below.SEM)
colnames(net.below.SEM)


lvmod.1 <- ' 
# Regressions
 Micro.rich ~ pre.AI 
 Micro.net ~ pre.AI + Micro.rich
 Micro.mass ~ pre.AI + Micro.rich + Micro.net
'
lvmod.1.fit <- sem(lvmod.1, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Micro.rich ~ pre.AI 
 Micro.net ~ pre.AI + Micro.rich
 Micro.mass ~ pre.AI  + Micro.rich #  + Micro.net
'
lvmod.2.fit <- sem(lvmod.2, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.3 <- ' 
# Regressions
 Micro.rich ~ pre.AI 
 Micro.net ~ pre.AI  + Micro.rich
 Micro.mass ~ pre.AI + Micro.net # + Micro.rich 
'
lvmod.3.fit <- sem(lvmod.3, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)

source("../2.DX2013_145sites/function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.1.fit,lvmod.2.fit,lvmod.3.fit),
              c("model1","model2","model3")) 


semPaths(lvmod.3.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)


summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.3.fit, c("cfi", "rmsea", "srmr"))




