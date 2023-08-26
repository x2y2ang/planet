
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# library
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(MASS)
library(plyr)
library(dplyr)

# plant species 2014
sp2014 = read.csv("../2.DX2013_145sites/data/plant78sites.species.abundance.2014.csv")

colnames(sp2014)[1] = c("Sites")

sp2014.me = melt(sp2014, id = c("Sites"))


# plant species height 2014
height2014 = read.csv("../2.DX2013_145sites/data/plant78sites.species.height.csv")

colnames(height2014)[1] = c("Sites")

height2014.me = melt(height2014, id = c("Sites"))

# plant species coverage 2014
coverage2014 = read.csv("../2.DX2013_145sites/data/plant78sites.species.coverage.2014.csv")

colnames(coverage2014)[1] = c("Sites")

coverage2014.me = melt(coverage2014, id = c("Sites"))

num = data.frame(rep(1:45, each = 78))


# plant biomass lm regression parameters
lm.re.para = read.csv("../2.DX2013_145sites/data/biomass.lm.regression.parameter.add.45species.2014.78sites.csv")
colnames(lm.re.para)[1] = c("species")

# combine all parameter
plant2014.all = cbind(num, height2014.me$value, sp2014.me$value,
                      coverage2014.me$value)

colnames(plant2014.all) = c("species", "spe.hei", "spe.num", "spe.cov")


plant.lm.re = merge(plant2014.all, lm.re.para, by = c("species"))

# spe.cov 
plant.lm.re1 = cbind(sp2014.me$Sites, plant.lm.re)
plant.lm.re2 = subset(plant.lm.re1, spe.hei != "0" & spe.num != "0" 
                     & spe.cov != "0")
colnames(plant.lm.re2)[1] = c("Sites")

# calculate standing crop(vegetation volume)
veg.vol = plant.lm.re2$a + (plant.lm.re2$b)*(plant.lm.re2$spe.cov)*(plant.lm.re2$spe.hei)
veg.vol = data.frame(veg.vol)

# cbind 
plant.vol.2014 = cbind(plant.lm.re2, veg.vol)

# calculate biomass
biomass2014 = (plant.vol.2014$c)*(plant.vol.2014$veg.vol)+plant.vol.2014$d


plant.mass2014 = cbind(plant.lm.re2, biomass2014)


mass2014 = plant.mass2014[order(plant.mass2014$Sites), ]


bioma = ddply(mass2014, c("Sites"), summarise, sum = sum(biomass2014))

spe = ddply(mass2014, c("Sites"), summarise, sum = sum(spe.num))

mass.spe = cbind(bioma, spe[,2])

# "Sites", "Biomass", "Species.number" of each site
colnames(mass.spe) = c("Sites", "Biomass", "Spe.num")

#### 2014 plant quadrat areas 1 m2 #####

#### standard with 2020 0.25 m2 ########
mass.spe.0.25m2 = (mass.spe[, 2:3]/4)

mass.spe.0.25m2 = cbind(mass.spe[, c("Sites")], mass.spe.0.25m2)

colnames(mass.spe.0.25m2)[1] =  c("Sites")


spe.num = ddply(mass2014.tmp, c("Sites", "species"), summarise, 
                sum = sum(spe.num))


spe.num0.25m2 = data.frame(round((spe.num[, 3]/4)))

spe.num0.25m2 = cbind(spe.num[, c("Sites", "species")], spe.num0.25m2)

colnames(spe.num0.25m2)[3] =  c("spe.num")


# save data
write.csv(mass.spe.0.25m2, "../2.DX2013_145sites/data/last.DX2014.plant.biomass.total.species.number.csv")
write.csv(spe.num0.25m2, "../2.DX2013_145sites/data/last.DX2014.plant.species.number.each.sites.csv")


