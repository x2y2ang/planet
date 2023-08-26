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

enviall = read.csv("../2.DX2013_145sites/data/micro.plant1314.alp.mass.envi.all.last.add.plant.rich.csv",
                     row.names = 1)

# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")

load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")


cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)


colnames(envi1314)
soil = envi1314[, c(11:19)]

pca = PCA(soil, scale.unit=T, graph=F)
summary(pca)
pca.dd <- pca$ind$coord


############# EMF ################
EMF = cbind(enviall[, c("TOC", "TN", "DOC", "DON", "TP", "Veg.mass", "PLFA")])
# z-score
EMF.tmp = scale(EMF)
# mean
EMF.la = data.frame(apply(EMF.tmp, 1, mean))
colnames(EMF.la) = c("EMF")

colnames(enviall)
sem.all = cbind(enviall[, c("Elevation", "pre.AI", "Micro.rich", "Plant.rich")],
                pca.dd[, 1], EMF.la)

colnames(sem.all)[5] = c("Soil.pca1")

write.csv(sem.all, "../2.DX2013_145sites/data/last.sem.plant.micro.rich.csv")

# aridity <= 0.69

all.below = subset(sem.all, Elevation > 4700)


net.below.data = all.below

net.below.SEM.md = na.omit(net.below.data)
net.below.SEM = scale(net.below.SEM.md,center=T,scale=T) 
net.below.SEM = as.data.frame(net.below.SEM)
net.below.SEM = data.table(net.below.SEM)
colnames(net.below.SEM)


lvmod.1 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Micro.rich ~ pre.AI + Plant.rich + Soil.pca1
 Soil.pca1 ~ pre.AI + Plant.rich + Micro.rich
 EMF ~ pre.AI + Plant.rich + Soil.pca1 + Micro.rich
'
lvmod.1.fit <- sem(lvmod.1, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Micro.rich ~ pre.AI + Plant.rich + Soil.pca1
 Soil.pca1 ~ pre.AI + Plant.rich + Micro.rich
 EMF ~ Plant.rich + Soil.pca1 + Micro.rich #pre.AI + 
'
lvmod.2.fit <- sem(lvmod.2, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.3 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Micro.rich ~ pre.AI + Soil.pca1 # + Plant.rich
 Soil.pca1 ~ pre.AI + Plant.rich + Micro.rich
 EMF ~ Plant.rich + Soil.pca1 + Micro.rich #pre.AI + 
'
lvmod.3.fit <- sem(lvmod.3, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.4 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Micro.rich ~ pre.AI + Soil.pca1 # + Plant.rich
 Soil.pca1 ~ pre.AI + Plant.rich # + Micro.rich
 EMF ~ Plant.rich + Soil.pca1 + Micro.rich #pre.AI + 
'
lvmod.4.fit <- sem(lvmod.4, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.4.fit, rsq=T, standardized=T,fit.measures = TRUE)


source("../2.DX2013_145sites/function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.2.fit, lvmod.3.fit, lvmod.4.fit),
              c("model2","model3","model4")) 


semPaths(lvmod.4.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)


summary(lvmod.4.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.4.fit, c("cfi", "rmsea", "srmr"))




