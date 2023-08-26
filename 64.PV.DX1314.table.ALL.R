
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# library
library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(splitstackshape)
library(FactoMineR)

# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")


load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")


source("../2.DX2013_145sites/function/diversity.functions.R")
dim(comm1314.samp)
micro1314 = quick.diversity(comm1314.samp)

cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

cli1314.last = cbind(cli1314.last, envi1314[, c("SWC")])

colnames(cli1314.last)[11] = c("SWC")


veg14 = read.csv("../2.DX2013_145sites/data/plant.2014.species.number.csv",
                 row.names = 1)

#### 2014 plant quadrat areas 1 m2 #####

#### standard with 2020 0.25 m2 ########
veg14 = round((veg14/4))

# na value equal 0
veg14[is.na(veg14)] <- 0

source("../2.DX2013_145sites/function/diversity.functions.R")

veg14.div = quick.diversity(veg14)
veg13.div = quick.diversity(Veg.sp)

# combine 1314
veg.div1314 = rbind(veg13.div, veg14.div)

envi188 = envi1314[rownames(veg.div1314), ]
cli188 = cli1314.last[rownames(veg.div1314), ]

################    loop for 10 eles    ######################
veg188.cli.tmp = cbind(cli188[, 1:6], envi188[, c(6:10)], veg.div1314$comm.richness)
colnames(veg188.cli.tmp)[12] = c("Plant.rich")

write.csv(veg188.cli.tmp, "../2.DX2013_145sites/data/DX1314.Plant.rich.AI.for.PV.csv")




colnames(envi1314)
soil = envi1314[, c(11:19)]

pca = PCA(soil, scale.unit=T, graph=F)
summary(pca)
pca.dd <- pca$ind$coord

################    loop for 10 eles    ######################
Soil.pca1.cli = cbind(cli1314.last, pca.dd[, 1])

colnames(Soil.pca1.cli)[12] = c("Soil.pca1")



############# EMF ################

envi.mass = read.csv("../2.DX2013_145sites/data/envi.1314.add.mass.csv",
                     row.names = 1)

EMF = cbind(envi1314[, c("TOC", "TN", "DOC", "DON", "TP")], 
            envi.mass[,c("Veg.mass", "PLFA")])
# z-score
EMF.tmp = scale(EMF)
# mean
EMF.la = data.frame(apply(EMF.tmp, 1, mean))
colnames(EMF.la) = c("EMF")

EMF.AI = data.frame(envi1314$Elevation,  EMF.la) 

colnames(EMF.AI) = c("Elevation", "EMF")



Plant.micro.soil.EMF = cbind(cli1314.last[, 1:6], envi1314[, c(6:10)],  micro1314$comm.richness, EMF.AI$EMF, envi.mass[, 2:3],
                             Soil.pca1.cli$Soil.pca1)

colnames(Plant.micro.soil.EMF)[12:16] = c("Micro.rich", "EMF", "AGB", "PLFA", "Soil.pca1")

write.csv(Plant.micro.soil.EMF, "../2.DX2013_145sites/data/DX1314.Plant.micro.soil.EMF.AI.for.PV.csv")


