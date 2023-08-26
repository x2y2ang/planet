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

# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")

load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")

cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

veg14 = read.csv("../2.DX2013_145sites/data/plant.2014.species.number.csv",
                 row.names = 1)

envi.mass = read.csv("../2.DX2013_145sites/data/envi.1314.add.mass.csv",
                     row.names = 1)


cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

source("../2.DX2013_145sites/function/diversity.functions.R")

micro.alp1314 = quick.diversity(comm1314.samp)

sem.all = cbind(cli1314.last[, 3], micro.alp1314, envi.mass)
colnames(sem.all)[1] = c("pre.AI")

write.csv(sem.all, "../2.DX2013_145sites/data/micro.plant1314.alp.mass.envi.all.last.csv")

veg14 = round((veg14/4))

# na value equal 0
veg14[is.na(veg14)] <- 0

# to match the data 2014
sp36.45 = matrix(data=0, nrow = 110, ncol = 10, byrow = FALSE, dimnames = NULL)
colnames(sp36.45) = c("sp36", "sp37", "sp38", "sp39","sp40",
                      "sp41", "sp42", "sp43", "sp44","sp45")

veg13 = cbind(Veg.sp, sp36.45)


# combine veg data 2013 and 2014
veg1314 = rbind(veg13, veg14)



plant1314.rich = quick.diversity(veg1314)

plant1314.rich = cbind(rownames(plant1314.rich), plant1314.rich)
colnames(plant1314.rich)[1] = c("ID")

net.plant = read.csv("../2.DX2013_145sites/data/plant1314.188sites.coverage.cohesion.all.csv",
                     row.names = 1)
net.plant = cbind(rownames(net.plant), net.plant)
colnames(net.plant)[1] = c("ID")

net.plant.last = merge(net.plant, plant1314.rich, by = c("ID"))
rownames(net.plant.last) = net.plant.last[, 1]

mass.plant = envi.mass[rownames(net.plant.last), ]

all.plant = cbind(mass.plant[, c(1:2)], net.plant.last[, c("Total.coh", "pre.AI", "comm.richness")])

colnames(all.plant)[5] = c("Plant.rich")
colnames(all.plant)[3] = c("Plant.net")

write.csv(all.plant, "../2.DX2013_145sites/data/plant1314.cov.network.mass.rich.for.SEM.csv")

# aridity <= 0.69

all.below = subset(all.plant, Elevation > 4700)


net.below.data = all.below

net.below.SEM.md = na.omit(net.below.data)
net.below.SEM = scale(net.below.SEM.md,center=T,scale=T) 
net.below.SEM = as.data.frame(net.below.SEM)
net.below.SEM = data.table(net.below.SEM)
colnames(net.below.SEM)


lvmod.1 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Plant.net ~ pre.AI + Plant.rich
 Veg.mass ~ pre.AI + Plant.rich + Plant.net
'
lvmod.1.fit <- sem(lvmod.1, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Plant.rich ~ pre.AI 
 Plant.net ~ pre.AI + Plant.rich
 Veg.mass ~  Plant.rich + Plant.net #pre.AI +
'
lvmod.2.fit <- sem(lvmod.2, data=net.below.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


source("../2.DX2013_145sites/function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.1.fit,lvmod.2.fit),
              c("model1","model2")) 


semPaths(lvmod.2.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)


summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.2.fit, c("cfi", "rmsea", "srmr"))




