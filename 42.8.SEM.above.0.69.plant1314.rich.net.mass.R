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

all.plant = read.csv("../2.DX2013_145sites/data/plant1314.cov.network.mass.rich.for.SEM.csv",
                     row.names = 1)


# aridity <= 0.69

all.above = subset(all.plant, Elevation <= 4700)


net.above.data = all.above

net.above.SEM.md = na.omit(net.above.data)
net.above.SEM = scale(net.above.SEM.md,center=T,scale=T) 
net.above.SEM = as.data.frame(net.above.SEM)
net.above.SEM = data.table(net.above.SEM)
colnames(net.above.SEM)


lvmod.1 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Plant.net ~ pre.AI + Plant.rich
 Veg.mass ~ pre.AI + Plant.rich + Plant.net
'
lvmod.1.fit <- sem(lvmod.1, data=net.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Plant.net ~ pre.AI # + Plant.rich
 Veg.mass ~ pre.AI + Plant.rich + Plant.net
'
lvmod.2.fit <- sem(lvmod.2, data=net.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.3 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Plant.net ~ pre.AI # + Plant.rich
 Veg.mass ~ pre.AI + Plant.rich # + Plant.net
'
lvmod.3.fit <- sem(lvmod.3, data=net.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)


source("../2.DX2013_145sites/function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.1.fit,lvmod.2.fit,lvmod.3.fit),
              c("model1","model2","model3")) 


semPaths(lvmod.2.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)


summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.2.fit, c("cfi", "rmsea", "srmr"))
