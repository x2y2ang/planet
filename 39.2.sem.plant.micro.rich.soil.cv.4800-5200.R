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

### 4800-5200 ###
cv.high = CV.all[c(1:79), ]
colnames(cv.high)

cv.high.data = cv.high

cv.high.SEM.md = na.omit(cv.high.data)
cv.high.SEM = scale(cv.high.SEM.md,center=T,scale=T) 
cv.high.SEM = as.data.frame(cv.high.SEM)
cv.high.SEM = data.table(cv.high.SEM)
colnames(cv.high.SEM)


lvmod.1 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI
 EF.CV ~ Plant.Rich.CV + Soil.PCA1.CV + Micro.Rich.CV + pre.AI
'
lvmod.1.fit <- sem(lvmod.1, data=cv.high.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.3 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI
 EF.CV ~ Soil.PCA1.CV + Micro.Rich.CV + pre.AI
'
lvmod.3.fit <- sem(lvmod.3, data=cv.high.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.4 <- ' 
# Regressions
 Plant.Rich.CV ~ pre.AI 
 Micro.Rich.CV ~ pre.AI + Plant.Rich.CV
 Soil.PCA1.CV ~ Plant.Rich.CV + pre.AI
 EF.CV ~ Plant.Rich.CV + Soil.PCA1.CV + pre.AI
'
lvmod.4.fit <- sem(lvmod.4, data=cv.high.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.4.fit, rsq=T, standardized=T,fit.measures = TRUE)



source("../2.DX2013_145sites/function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.1.fit,lvmod.3.fit,lvmod.4.fit),
              c("model1","model3","model4")) 


semPaths(lvmod.1.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)

summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.1.fit, c("cfi", "rmsea", "srmr"))



