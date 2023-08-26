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

CV.all = read.csv("../2.DX2013_145sites/data/plant.micro.soil.CV.100sites.csv")
colnames(CV.all)[1] = c("Elecvation")

### 4300-4700 ###
cv.low = CV.all[c(80:132), ]
colnames(cv.low)

cv.low.data = cv.low

cv.low.SEM.md = na.omit(cv.low.data)
cv.low.SEM = scale(cv.low.SEM.md,center=T,scale=T) 
cv.low.SEM = as.data.frame(cv.low.SEM)
cv.low.SEM = data.table(cv.low.SEM)
colnames(cv.low.SEM)


lvmod.1 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI
 EF.CV ~ Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV + pre.AI
'
lvmod.1.fit <- sem(lvmod.1, data=cv.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI
 EF.CV ~ Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV + pre.AI
'
lvmod.2.fit <- sem(lvmod.2, data=cv.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.3 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV + Soil.PCA1.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI 
 EF.CV ~ Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV + pre.AI
'
lvmod.3.fit <- sem(lvmod.3, data=cv.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.4 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV + Soil.PCA1.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI 
 EF.CV ~ Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV + pre.AI
'
lvmod.4.fit <- sem(lvmod.4, data=cv.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.4.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.5 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Soil.PCA1.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI 
 EF.CV ~ Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV + pre.AI
'
lvmod.5.fit <- sem(lvmod.5, data=cv.low.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.5.fit, rsq=T, standardized=T,fit.measures = TRUE)


source("../2.DX2013_145sites/function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.3.fit,lvmod.4.fit,lvmod.5.fit),
              c("model3","model4","model5")) 

semPaths(lvmod.5.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)

summary(lvmod.5.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.5.fit, c("cfi", "rmsea", "srmr"))



