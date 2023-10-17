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

CV.all = read.csv("../data/GCB/plant.micro.soil.CV.100sites.csv")
colnames(CV.all)[1] = c("Elecvation")

# aridity < 0.69

cv.below = subset(CV.all, pre.AI < 0.69)

cv.below.data = cv.below

cv.below.SEM.md = na.omit(cv.below.data)
cv.below.SEM = scale(cv.below.SEM.md,center=T,scale=T) 
cv.below.SEM = as.data.frame(cv.below.SEM)
cv.below.SEM = data.table(cv.below.SEM)
colnames(cv.below.SEM)


lvmod.1 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI
 EF.CV ~ Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV + pre.AI
'
lvmod.1.fit <- sem(lvmod.1, data=cv.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI
 EF.CV ~ Soil.PCA1.CV + Micro.Rich.CV + pre.AI # Plant.Rich.CV
'
lvmod.2.fit <- sem(lvmod.2, data=cv.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


source("../function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.1.fit,lvmod.2.fit),
              c("model1","model2")) 


semPaths(lvmod.1.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)

summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.1.fit, c("cfi", "rmsea", "srmr"))



