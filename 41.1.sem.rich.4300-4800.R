
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
library(reshape2)
library(RColorBrewer)
library(dplyr)

# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")

load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")


cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

veg14 = read.csv("../2.DX2013_145sites/data/plant.2014.species.number.csv",
                 row.names = 1)

# na value equal 0
veg14[is.na(veg14)] <- 0

source("../2.DX2013_145sites/function/diversity.functions.R")

veg14.div = quick.diversity(veg14)
veg13.div = quick.diversity(Veg.sp)

micro1314.div = quick.diversity(comm1314.samp)

# combine 1314
veg.div1314 = rbind(veg13.div, veg14.div)

envi188 = envi1314[rownames(veg.div1314), ]
veg.div1314.ele = cbind(envi188$Elevation, veg.div1314)
colnames(veg.div1314.ele)[1] = c("Elevation")

# cli188 = cli1314.last[rownames(veg.div1314), ]

?left_join
A = left_join(envi1314, veg.div1314.ele, by=c("Elevation"))

envi1314 %>% left_join(veg.div1314.ele)

### merge the two table by contig name ###
envi1314.tmp = cbind(rownames(envi1314), envi1314)
colnames(envi1314.tmp)[1] = c("ID")
veg.div1314.ele.tmp = cbind(rownames(veg.div1314.ele), veg.div1314.ele)
colnames(veg.div1314.ele.tmp)[1] = c("ID")
seq.taxa = merge(envi1314.tmp, veg.div1314.ele.tmp, by = 'ID') 

# soil PCA
colnames(envi1314)
soil = envi1314[, c(11:19)]
pca = PCA(soil, scale.unit=T, graph=F)
summary(pca)
pca.dd <- pca$ind$coord

# for microbial richness
vec.tmp = cbind(cli1314.last[, c("Elevation")], 
                micro1314.div[, c("comm.richness"),],
                pca.dd[,1]) 
vec.tmp = data.frame(vec.tmp)

colnames(vec.tmp) = c("Elevation", "Micro.rich", "Soil.PCA1")
