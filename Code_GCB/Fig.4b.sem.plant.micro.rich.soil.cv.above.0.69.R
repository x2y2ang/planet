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

# aridity >= 0.69

cv.above = subset(CV.all, pre.AI >= 0.69)

colnames(cv.above)

cv.above.data = cv.above

cv.above.SEM.md = na.omit(cv.above.data)
cv.above.SEM = scale(cv.above.SEM.md,center=T,scale=T) 
cv.above.SEM = as.data.frame(cv.above.SEM)
cv.above.SEM = data.table(cv.above.SEM)
colnames(cv.above.SEM)


lvmod.1 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV + Soil.PCA1.CV
 Soil.PCA1.CV ~ pre.AI + Plant.Rich.CV + Micro.Rich.CV
 EF.CV ~ pre.AI + Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV 
'
lvmod.1.fit <- sem(lvmod.1, data=cv.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV + Soil.PCA1.CV
 Soil.PCA1.CV ~ pre.AI + Plant.Rich.CV # + Micro.Rich.CV
 EF.CV ~ pre.AI + Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV 
'
lvmod.2.fit <- sem(lvmod.2, data=cv.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.3 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV + Soil.PCA1.CV
 Soil.PCA1.CV ~ pre.AI + Plant.Rich.CV # + Micro.Rich.CV
 EF.CV ~ pre.AI + Plant.Rich.CV + Soil.PCA1.CV # + Micro.Rich.CV
'
lvmod.3.fit <- sem(lvmod.3, data=cv.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.4 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV + Soil.PCA1.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI 
 EF.CV ~ Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV + pre.AI
'
lvmod.4.fit <- sem(lvmod.4, data=cv.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.4.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.5 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Soil.PCA1.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI 
 EF.CV ~ Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV + pre.AI
'
lvmod.5.fit <- sem(lvmod.5, data=cv.above.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.5.fit, rsq=T, standardized=T,fit.measures = TRUE)


source("../function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.2.fit,lvmod.3.fit,lvmod.4.fit,lvmod.5.fit),
              c("model2","model3","model4","model5")) 

semPaths(lvmod.5.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)

summary(lvmod.5.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.5.fit, c("cfi", "rmsea", "srmr"))



