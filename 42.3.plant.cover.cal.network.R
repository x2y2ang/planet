
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(vegan)
library(reshape2)
library(plyr)
library(dplyr)
library(openxlsx)
library(WGCNA)
library(igraph)
library(psych)
library(vegan)
library(FactoMineR)
library(plyr)
library(igraph)
library(impute)
library(GO.db)
library(preprocessCore)
library(AnnotationDbi)
library(psych)
library(Hmisc)
library(splitstackshape)


plant <- read.csv("../2.DX2013_145sites/data/Plant1314.each.188sites.interactions.saiz.method.csv")

colnames(plant)

plant1314 <- melt(plant, id=c("Sites", "row"))

#split character to multiple columns
split1 <- strsplit(as.character(plant1314$value), " ")
split1

# split species coverage
spe.cov <- sapply(split1, "[", 2)

spe <- sapply(split1, "[", 1)
spe = data.frame(spe)

split2 <- strsplit(as.character(spe$spe), "/")
split2

# split species 
species <- sapply(split2, "[", 1)

# split species height
spe.hei <- sapply(split2, "[", 2)

# split species number
spe.num <- sapply(split2, "[", 3)

# plant properties
plant.tmp <- cbind(plant1314$Sites, species, spe.hei, spe.num, spe.cov)

plant.tmp2 = na.omit(plant.tmp)

colnames(plant.tmp2)[1] = c("Sites")

plant.tmp2 = data.frame(plant.tmp2)

dim(plant.tmp2)

# write.xlsx(plant.tmp2, "../2.DX2013_145sites/data/Split.Plant1314.each.188sites.interactions.saiz.method.xlsx")

write.csv(plant.tmp2, "../2.DX2013_145sites/data/Split.Plant1314.each.188sites.interactions.saiz.method.csv")


plant.tmp3 = as.data.frame(lapply(plant.tmp2[, c(2:5)], as.numeric))

plant.tmp4 = cbind(plant.tmp2[,1], plant.tmp3)
colnames(plant.tmp4)[1] = c("Sites") 


plant.tmp5 = na.omit(plant.tmp4)
  
spe.num = ddply(plant.tmp5, c("Sites", "species"), summarise, 
                sum = sum(spe.num))

spe.cov = ddply(plant.tmp5, c("Sites", "species"), summarise, 
                sum = sum(spe.cov))

spe.site.num = as.data.frame(table(plant.tmp5$Sites))

plant.fre = data.frame(c(rep(1, 24550)))

plant.last = cbind(plant.tmp5, plant.fre)
colnames(plant.last)[6] = c("fre")

spe.num.fre = ddply(plant.last, c("Sites", "species"), summarise, 
                sum = sum(fre))

real.cov = data.frame((spe.cov$sum)/(spe.num.fre$sum))

Plant.RA.cov = cbind(spe.num.fre, real.cov)

colnames(Plant.RA.cov)[3:4] = c("Number", "Coverage")

write.csv(Plant.RA.cov, "../2.DX2013_145sites/data/Plant1314.188sites.species.number.coverage.saiz.method.csv")

plant = dcast(Plant.RA.cov, Sites + species ~ Number)



# # max(spe.num.fre$sum)
# 
# # subset(spe.num.fre, sum >= 100)
# 
# plant.tmp3 = plant.tmp2[, c(1,2,5)]
# 
# HDX007 = subset(data.frame(plant.tmp3), Sites == "HDX007")
# 
# HDX007 = HDX007[order(HDX007$species, decreasing = ), ]
# 
# colnames(HDX007)
# 
# 
# sp.ii = 12
# for (sp.ii in unique(HDX007$species)) {
#   colnames(HDX007)
#   tp.clim.tmp = subset(HDX007, species == sp.ii)
#   print(sp.ii)
#   
#   rownames(tp.clim.tmp) = rep(c(1:nrow(tp.clim.tmp)))
#   rownames(tp.clim.tmp)
#   
#   tp.clim.tmp.t = cbind(rownames(tp.clim.tmp), tp.clim.tmp)
#   colnames(tp.clim.tmp.t)[1] = c("id")
#  
#       if (sp.ii == 12) {
#     env.out2 = tp.clim.tmp.t
#   } else {
#     env.out2 = merge(env.out2, tp.clim.tmp.t, by = c("id"), all = TRUE)
#     
#   }
# }
# 
# colnames(env.out2)
# cov007 = env.out2[, c("spe.cov.x", "spe.cov.y")]
#   
# 
# ### two loop ###
# plant.tmp3 = data.frame(plant.tmp3)
# unique(plant.tmp3$Sites)
# 
# site.jj = "HDX007"
# 
# for (site.jj in unique(plant.tmp3$Sites)) {
# 
# sp.ii = 12
# for (sp.ii in unique(plant.tmp3$species)) {
#   colnames(plant.tmp3)
#   tp.clim.tmp = subset(plant.tmp3, species == sp.ii)
#   print(sp.ii)
#   
#   rownames(tp.clim.tmp) = rep(c(1:nrow(tp.clim.tmp)))
#   rownames(tp.clim.tmp)
#   
#   tp.clim.tmp.t = cbind(rownames(tp.clim.tmp), tp.clim.tmp)
#   colnames(tp.clim.tmp.t)[1] = c("id")
#   
#   if (sp.ii == 12) {
#     env.out2 = tp.clim.tmp.t
#   } else {
#     env.out2 = merge(env.out2, tp.clim.tmp.t, by = c("id"), all = TRUE)
#     
#   }
# }
# 
#   if (site.jj == "HDX007"){
#     env.out.site = env.out2
#   } else {
#     env.out.site = cbind(env.out.site, env.out2)
#     
#   }
#     
# }



# sp12 = subset(HDX007, species == "12")
# sp18 = subset(HDX007, species == "18")
# 
# rownames(sp12) = rep(c(1:nrow(sp12)))
# rownames(sp12)
# 
# sp12.t = cbind(rownames(sp12), sp12)
# colnames(sp12.t)[1] = c("id")
# 
# rownames(sp18) = rep(c(1:nrow(sp18)))
# rownames(sp18)
# 
# sp18.t = cbind(rownames(sp18), sp18)
# colnames(sp18.t)[1] = c("id")
# 
# env.out1 = merge(sp12.t, sp18.t, by = c("id"), all = TRUE)
# 
# 
