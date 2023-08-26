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

Asy.all = read.csv("../2.DX2013_145sites/data/plant.micro.soil.Asy.100sites.csv")
colnames(Asy.all)[1] = c("pre.AI")

### 4300-4700 ###
Asy.low = Asy.all[c(80:132), ]
colnames(Asy.low)

Asy.low.data = Asy.low

Asy.low.SEM.md = na.omit(Asy.low.data)
Asy.low.SEM = scale(Asy.low.SEM.md,center=T,scale=T) 
Asy.low.SEM = as.data.frame(Asy.low.SEM)
Asy.low.SEM = data.table(Asy.low.SEM)
colnames(Asy.low.SEM)


lvmod.1 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy + Micro.Rich.Asy
'
lvmod.1.fit <- sem(lvmod.1, data=Asy.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy + Micro.Rich.Asy
'
lvmod.2.fit <- sem(lvmod.2, data=Asy.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.3 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy # + Micro.Rich.Asy
'
lvmod.3.fit <- sem(lvmod.3, data=Asy.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.4 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy # + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy # + Micro.Rich.Asy
'
lvmod.4.fit <- sem(lvmod.4, data=Asy.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.4.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.5 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy # + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy # + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy # + Micro.Rich.Asy
'
lvmod.5.fit <- sem(lvmod.5, data=Asy.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.5.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.6 <- ' 
# Regressions
 Plant.Rich.Asy ~ pre.AI 
 Micro.Rich.Asy ~ pre.AI + Plant.Rich.Asy # + Soil.PCA1.Asy
 Soil.PCA1.Asy ~ pre.AI + Plant.Rich.Asy # + Micro.Rich.Asy
 EF.Asy ~ pre.AI + Plant.Rich.Asy + Soil.PCA1.Asy + Micro.Rich.Asy
'
lvmod.6.fit <- sem(lvmod.6, data=Asy.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.6.fit, rsq=T, standardized=T,fit.measures = TRUE)


source("../2.DX2013_145sites/function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.2.fit,lvmod.3.fit,lvmod.4.fit,lvmod.5.fit, lvmod.6.fit),
              c("model2","model3","model4","model5", "model6")) 

semPaths(lvmod.6.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)

summary(lvmod.6.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.6.fit, c("cfi", "rmsea", "srmr"))



