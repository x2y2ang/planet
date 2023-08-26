
rm(list=ls())

library(vegan)
library(SpatialEpi)
library(SciViews)
library(segmented)
library(ggplot2)
library(MASS)
library(reshape2)
library(dplyr)
library(plyr)

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")

load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")


#### micro comm ####
geo = envi1314[,c("Longtitude", "Latitude")]
ele <- unique(envi1314[,"Elevation"])
length(ele)

i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  
  otu.tmp <- comm1314.samp[which(envi1314$Elevation==ele.tmp),]
  otu.dist.tmp = as.dist(vegdist(otu.tmp, method = "bray")) 
  
  geo.tmp <- geo[which(envi1314$Elevation==ele.tmp),]
  
  spatial.transform = latlong2grid(geo.tmp[,c("Longtitude","Latitude")])
  spatial.dist.tmp = as.dist(vegdist(spatial.transform,method = "euclid"))
  
  spa.df = data.frame(otu.d = as.vector(otu.dist.tmp),
                      spa.d = as.vector(spatial.dist.tmp))
  
  # all sites
  geo.lm.tmp = summary(lm(otu.d ~ ln(spa.d), spa.df))
  geo.all.lm = data.frame(ele = ele[i],
                          slope = geo.lm.tmp$coefficients[2,1],
                          Ine = geo.lm.tmp$coefficients[1,1],
                          p = geo.lm.tmp$coefficients[2,4],
                          r2 = geo.lm.tmp$adj.r.squared)
  
  table.sig.tmp = rbind(geo.all.lm) 
  
  if(i==1){
    micro.table.sig = table.sig.tmp
  }
  else{
    micro.table.sig = rbind(micro.table.sig,table.sig.tmp)
  }
}

plot(micro.table.sig$ele, micro.table.sig$slope)

### micro rich ###

source("../2.DX2013_145sites/function/diversity.functions.R")

micro1314.div = quick.diversity(comm1314.samp)

micro1314.rich = data.frame(micro1314.div[, 1])

i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  
  otu.tmp <- micro1314.div[which(envi1314$Elevation==ele.tmp),]
  otu.dist.tmp = as.dist(vegdist(otu.tmp, method = "bray")) 
  
  geo.tmp <- geo[which(envi1314$Elevation==ele.tmp),]
  
  spatial.transform = latlong2grid(geo.tmp[,c("Longtitude","Latitude")])
  spatial.dist.tmp = as.dist(vegdist(spatial.transform,method = "euclid"))
  
  spa.df = data.frame(otu.d = as.vector(otu.dist.tmp),
                      spa.d = as.vector(spatial.dist.tmp))
  
  # all sites
  geo.lm.tmp = summary(lm(otu.d ~ ln(spa.d), spa.df))
  geo.all.lm = data.frame(ele = ele[i],
                          slope = geo.lm.tmp$coefficients[2,1],
                          Ine = geo.lm.tmp$coefficients[1,1],
                          p = geo.lm.tmp$coefficients[2,4],
                          r2 = geo.lm.tmp$adj.r.squared)
  
  table.sig.tmp = rbind(geo.all.lm) 
  
  if(i==1){
    micro.rich.table.sig = table.sig.tmp
  }
  else{
    micro.rich.table.sig = rbind(micro.rich.table.sig,table.sig.tmp)
  }
}

plot(micro.rich.table.sig$ele, micro.rich.table.sig$slope)


#### for soil nutrient ###
library(FactoMineR)
soil = envi1314[, c("TOC", "TN", "DOC", "DON", "TP")]
soil.pca = PCA(soil, scale.unit=T, graph=F)
summary(soil.pca)
# Dim.1 66.681
soil.pca.dd <- data.frame(soil.pca$ind$coord)
soil.last = soil.pca.dd[, 1]
soil.last = data.frame(soil.last)


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  
  otu.tmp <- soil.last[which(envi1314$Elevation==ele.tmp),]
  otu.dist.tmp = as.dist(vegdist(otu.tmp, method = "euclidean")) 
  
  geo.tmp <- geo[which(envi1314$Elevation==ele.tmp),]
  
  spatial.transform = latlong2grid(geo.tmp[,c("Longtitude","Latitude")])
  spatial.dist.tmp = as.dist(vegdist(spatial.transform,method = "euclid"))
  
  spa.df = data.frame(otu.d = as.vector(otu.dist.tmp),
                      spa.d = as.vector(spatial.dist.tmp))
  
  # all sites
  geo.lm.tmp = summary(lm(otu.d ~ ln(spa.d), spa.df))
  geo.all.lm = data.frame(ele = ele[i],
                          slope = geo.lm.tmp$coefficients[2,1],
                          Ine = geo.lm.tmp$coefficients[1,1],
                          p = geo.lm.tmp$coefficients[2,4],
                          r2 = geo.lm.tmp$adj.r.squared)
  
  table.sig.tmp = rbind(geo.all.lm) 
  
  if(i==1){
    soil.table.sig = table.sig.tmp
  }
  else{
    soil.table.sig = rbind(soil.table.sig,table.sig.tmp)
  }
}

plot(soil.table.sig$ele, soil.table.sig$slope)


decay.micro.soil = rbind(micro.table.sig, soil.table.sig)

gp = data.frame(c(rep(c("Micro.decay"), 10), rep(c("Soil.decay"), 10)))

decay.micro.soil1 = cbind(gp, decay.micro.soil)

write.csv(decay.micro.soil1, "../2.DX2013_145sites/data/decay.micro.soil1.last.csv")

plot(envi$Elevation, decay.micro.soil1$slope)

### plant comm ###
veg14 = read.csv("../2.DX2013_145sites/data/plant.2014.species.number.csv",
                 row.names = 1)

#### standard with 2020 0.25 m2 ########
veg14 = round((veg14/4))


# na value equal 0
veg14[is.na(veg14)] <- 0

write.csv(Veg.sp, "../2.DX2013_145sites/data/plant2013.comm.csv")

veg13.sp = read.csv("../2.DX2013_145sites/data/plant2013.comm45species.csv",
                    row.names = 1)

veg1314.sp = rbind(veg13.sp, veg14)

envi1314.tmp = envi1314[rownames(veg1314.sp), ]


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  
  otu.tmp <- veg1314.sp[which(envi1314.tmp$Elevation==ele.tmp),]
  otu.dist.tmp = as.dist(vegdist(otu.tmp, method = "bray")) 
  
  geo.tmp <- geo[which(envi1314.tmp$Elevation==ele.tmp),]
  
  spatial.transform = latlong2grid(geo.tmp[,c("Longtitude","Latitude")])
  spatial.dist.tmp = as.dist(vegdist(spatial.transform,method = "euclid"))
  
  spa.df = data.frame(otu.d = as.vector(otu.dist.tmp),
                      spa.d = as.vector(spatial.dist.tmp))
  
  # all sites
  geo.lm.tmp = summary(lm(otu.d ~ ln(spa.d), spa.df))
  geo.all.lm = data.frame(ele = ele[i],
                          slope = geo.lm.tmp$coefficients[2,1],
                          Ine = geo.lm.tmp$coefficients[1,1],
                          p = geo.lm.tmp$coefficients[2,4],
                          r2 = geo.lm.tmp$adj.r.squared)
  
  table.sig.tmp = rbind(geo.all.lm) 
  
  if(i==1){
   plant.table.sig = table.sig.tmp
  }
  else{
    plant.table.sig = rbind(plant.table.sig,table.sig.tmp)
  }
}

plot(plant.table.sig$ele, plant.table.sig$slope)




### plant cover ###
envi107 = envi1314[rownames(Veg.cov), ]
plot(envi107$Elevation, Veg.cov$Veg.cov)


cov = cbind(envi107$Elevation, Veg.cov)
colnames(cov)[1] = c("Elevation")



spe.cov = ddply(cov, c("Elevation"), summarise, 
                mean = mean(Veg.cov))

plot(plant.table.sig$slope, spe.cov$mean)

plant.cov.decay = cbind(plant.table.sig$slope, spe.cov)

colnames(plant.cov.decay) = c("Plant.decay", "Elevation", "Plant.cov")

write.csv(plant.cov.decay, "../2.DX2013_145sites/data/plant.cov.decay.csv")

plant.cov.decay$Elevation = factor(plant.cov.decay$Elevation)

library(RColorBrewer)


lm.plant = lm(Plant.decay ~ Plant.cov, plant.cov.decay)
summary(lm.plant)
# Adjusted R-squared:  0.3464, p-value: 0.04305

qu.plant = lm(Plant.decay ~ Plant.cov + I(Plant.cov^2), plant.cov.decay)
summary(qu.plant)

AIC(lm.plant, qu.plant)

AI = data.frame(c("0.99", "0.94", "0.86", "0.78", "0.69", "0.58", "0.45",
                  "0.3", "0.15", "0.002"))

plant.cov.decay = cbind(AI, plant.cov.decay)
colnames(plant.cov.decay)[1] = c("AI")
plant.cov.decay$AI = factor(plant.cov.decay$AI)

p1 = ggplot(plant.cov.decay, aes(x=Plant.cov, y=Plant.decay, color = AI))+ 

  geom_point(size = 2)+
  stat_smooth(method="lm", color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+  
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(strip.background = element_blank())+
  theme(axis.text.x  = element_text(vjust=0.5))+
  theme(text = element_text(size = 12))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))+
  xlab("Plant.cov")+ylab("Plant.decay")

p1

# micro composition
comm.hel.B = decostand(comm1314.samp, "hellinger")
# bray
euc <- vegdist(comm.hel.B, method="bray")
# cmdscale
b.pcoa1 <- cmdscale(euc,k=(nrow(comm.hel.B)-1),eig=TRUE)
b.pcoa1
summary(b.pcoa1)
# The first two principal coordinates of each point
micro.pcoa <- as.data.frame(b.pcoa1$points[,1:2])
# change names
names(micro.pcoa) <- c("PCoA1","PCoA2")
pcoa.eig.b <- (b.pcoa1$eig)[1:2] / sum(b.pcoa1$eig)

micro.pcoa.ele <- cbind(envi1314[,c("Elevation")], 
                        micro.pcoa[,1])
micro.pcoa.ele <- data.frame(micro.pcoa.ele)

colnames(micro.pcoa.ele) = c("Elevation", "Micro.PCoA1")

Micro.PCoA1 = ddply(micro.pcoa.ele, c("Elevation"), summarise, 
                mean = mean(Micro.PCoA1))


micro.pcoa.dacay = cbind(micro.table.sig$slope, Micro.PCoA1)

colnames(micro.pcoa.dacay) = c("Micro.decay", "Elevation", "Micro.pcoa1")

write.csv(micro.pcoa.dacay, "../2.DX2013_145sites/data/micro.pcoa.dacay.csv")

lm.micro = lm(Micro.pcoa1 ~ Micro.decay, micro.pcoa.dacay)
summary(lm.micro)

qu.micro = lm(Micro.pcoa1 ~ Micro.decay + I(Micro.decay^2), micro.pcoa.dacay)
summary(qu.micro)
#Adjusted R-squared:  0.1549 p-value: 0.2302

AIC(lm.micro, qu.micro)

micro.pcoa.dacay$Elevation = factor(micro.pcoa.dacay$Elevation)

AI = data.frame(c("0.99", "0.94", "0.86", "0.78", "0.69", "0.58", "0.45",
                  "0.3", "0.15", "0.002"))

micro.pcoa.dacay = cbind(AI, micro.pcoa.dacay)
colnames(micro.pcoa.dacay)[1] = c("AI")
micro.pcoa.dacay$AI = factor(micro.pcoa.dacay$AI)

p2 = ggplot(micro.pcoa.dacay, aes(x=Micro.decay, y=Micro.pcoa1, color = AI))+ 
  
  geom_point(size = 2)+
  stat_smooth(method="lm", formula=y~x+I(x^2),  color = "black", lty = 2)+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+  
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(strip.background = element_blank())+
  theme(axis.text.x  = element_text(vjust=0.5))+
  theme(text = element_text(size = 12))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))+
  xlab("Micro.decay")+ylab("Micro.pcoa1")

p2 

library(cowplot)
plot_grid(p1,p2)
ggsave("../2.DX2013_145sites/figures/aridity.plant.cov.decay.micro.pcoa.decay.pdf",
       width = 10, height = 5)


ggsave("../2.DX2013_145sites/figures/AI.plant.cov.decay.micro.pcoa.decay1.pdf",
       width = 7.5, height = 3)

# plant.micro = read.csv("../2.DX2013_145sites/data/plant.micro.pcoa.decay.AI.csv")
# 
# colnames(plant.micro)[1] = c("Group")
# 
# colnames(plant.micro)
# 
# plant.micro$AI = factor(plant.micro$AI)
# 
# ggplot(plant.micro, aes(x=Decay, y=Beta, color = AI))+ 
#   
#   geom_point(size = 2)+
#   stat_smooth(method="lm", color = "black")+
#   scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
#   facet_wrap( ~ Group, ncol = 2, scales = "free")
#   theme_bw()+  
#   theme(legend.background=element_rect(colour="Black",size=0.5))+
#   theme(strip.background = element_blank())+
#   theme(axis.text.x  = element_text(vjust=0.5))+
#   theme(text = element_text(size = 12))+
#   theme( panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          axis.line = element_line(colour = "black"))+
#   theme(plot.title = element_text(hjust = 0.5,size = 12),
#         axis.text=element_text(size=12),
#         axis.title.x=element_text(size=12),
#         axis.title.y=element_text(size=12))+
#   xlab("Plant.cov")+ylab("Plant.decay")
# 
# 
# 
# ### plant 2013 ###
# veg13.sp = read.csv("../2.DX2013_145sites/data/plant2013.comm45species.csv",
#                     row.names = 1)
# 
# 
# 
# envi13.tmp = envi[rownames(veg13.sp), ]
# 
# 
# i=1
# for (i in 1:length(ele)){
#   ele.tmp <- ele[i]
#   
#   otu.tmp <- veg13.sp[which(envi13.tmp$Elevation==ele.tmp),]
#   otu.tmp1 = na.omit(otu.tmp)
#   otu.dist.tmp = as.dist(vegdist(otu.tmp, method = "bray")) 
#   
#   geo.tmp <- geo[which(envi13.tmp$Elevation==ele.tmp),]
#   
#   spatial.transform = latlong2grid(geo.tmp[,c("Longtitude","Latitude")])
#   spatial.dist.tmp = as.dist(vegdist(spatial.transform,method = "euclid"))
#   
#   spa.df = data.frame(otu.d = as.vector(otu.dist.tmp),
#                       spa.d = as.vector(spatial.dist.tmp))
#   
#   # all sites
#   geo.lm.tmp = summary(lm(otu.d ~ ln(spa.d), spa.df))
#   geo.all.lm = data.frame(ele = ele[i],
#                           slope = geo.lm.tmp$coefficients[2,1],
#                           Ine = geo.lm.tmp$coefficients[1,1],
#                           p = geo.lm.tmp$coefficients[2,4],
#                           r2 = geo.lm.tmp$adj.r.squared)
#   
#   table.sig.tmp = rbind(geo.all.lm) 
#   
#   if(i==1){
#     plant13.table.sig = table.sig.tmp
#   }
#   else{
#     plant13.table.sig = rbind(plant13.table.sig,table.sig.tmp)
#   }
# }
# 
# plot(plant13.table.sig$ele, plant13.table.sig$slope)


