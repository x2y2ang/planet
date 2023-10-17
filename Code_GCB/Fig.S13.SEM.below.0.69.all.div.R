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

enviall = read.csv("../data/GCB/micro.plant1314.alp.mass.envi.all.last.add.plant.rich.csv",
                     row.names = 1)

# load data
load("../data/GCB/resample.envi1314.comm1314.taxa1314.Rdata")

load("../data/GCB/1.comm.envi.div.all.Rdata")


cli1314.last = read.csv("../data/GCB/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)


colnames(envi1314)
soil = envi1314[, c("TOC", "TN", "DOC", "DON", "TP")]

pca = PCA(soil, scale.unit=T, graph=F)
summary(pca)
pca.dd <- pca$ind$coord


############# EF ################
EF = cbind(enviall[, c("TOC", "TN", "DOC", "DON", "TP", "Veg.mass", "PLFA")])
# z-score
EF.tmp = scale(EF)
# mean
EF.la = data.frame(apply(EF.tmp, 1, mean))
colnames(EF.la) = c("EF")

colnames(enviall)
sem.all = cbind(enviall[, c("Elevation", "pre.AI", "Micro.rich", "Plant.rich")],
                pca.dd[, 1], EF.la)

colnames(sem.all)[5] = c("Soil.pca1")

# aridity < 0.69

all.below = subset(sem.all, pre.AI < 0.69)


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
 EF ~ pre.AI + Plant.rich + Soil.pca1 + Micro.rich
'
lvmod.1.fit <- sem(lvmod.1, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Micro.rich ~ pre.AI + Soil.pca1 # + Plant.rich
 Soil.pca1 ~ pre.AI + Plant.rich + Micro.rich
 EF ~ pre.AI + Plant.rich + Soil.pca1 + Micro.rich
'
lvmod.2.fit <- sem(lvmod.2, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.3 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Micro.rich ~ pre.AI + Soil.pca1 # + Plant.rich
 Soil.pca1 ~ pre.AI + Plant.rich + Micro.rich
 EF ~ pre.AI + Plant.rich + Soil.pca1 # + Micro.rich
'
lvmod.3.fit <- sem(lvmod.3, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)

source("../function/lavaan.modavg.R")

aictab.lavaan(list(lvmod.2.fit, lvmod.3.fit),
              c("model2","model3")) 

semPaths(lvmod.3.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)

summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.3.fit, c("cfi", "rmsea", "srmr"))




