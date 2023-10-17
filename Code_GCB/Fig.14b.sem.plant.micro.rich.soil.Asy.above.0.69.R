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

Asy.all = read.csv("../data/GCB/plant.micro.soil.Asy.100sites.csv")
colnames(Asy.all)[1] = c("Elevation")

# aridity >= 0.69

Asy.above = subset(Asy.all, pre.AI >= 0.69)

colnames(Asy.above)

Asy.above.data = Asy.above

Asy.above.SEM.md = na.omit(Asy.above.data)
Asy.above.SEM = scale(Asy.above.SEM.md,center=T,scale=T) 
Asy.above.SEM = as.data.frame(Asy.above.SEM)
Asy.above.SEM = data.table(Asy.above.SEM)
colnames(Asy.above.SEM)


lvmod.1 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy + Micro.Rich.Asy
'
lvmod.1.fit <- sem(lvmod.1, data=Asy.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy # + Micro.Rich.Asy
'
lvmod.2.fit <- sem(lvmod.2, data=Asy.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.3 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy # + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy # + Micro.Rich.Asy
'
lvmod.3.fit <- sem(lvmod.3, data=Asy.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.4 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy # + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy # + Micro.Rich.Asy
'
lvmod.4.fit <- sem(lvmod.4, data=Asy.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.4.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.5 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy # + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy # + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy # + Micro.Rich.Asy
'
lvmod.5.fit <- sem(lvmod.5, data=Asy.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.5.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.6 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy # + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy # + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy + Micro.Rich.Asy
'
lvmod.6.fit <- sem(lvmod.6, data=Asy.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.6.fit, rsq=T, standardized=T,fit.measures = TRUE)


source("../function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.2.fit,lvmod.3.fit,lvmod.4.fit,lvmod.5.fit, lvmod.6.fit),
              c("model2","model3","model4","model5", "model6")) 

semPaths(lvmod.3.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)

summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.3.fit, c("cfi", "rmsea", "srmr"))



