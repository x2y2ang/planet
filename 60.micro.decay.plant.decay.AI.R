
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
library(RColorBrewer)

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



### plant comm ###
veg14 = read.csv("../2.DX2013_145sites/data/plant.2014.species.number.csv",
                 row.names = 1)

#### standard with 2020 0.25 m2 ########
veg14 = round((veg14/4))


# na value equal 0
veg14[is.na(veg14)] <- 0

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


AI = data.frame(c("0.99", "0.94", "0.86", "0.78", "0.69", "0.58", "0.45",
                  "0.3", "0.15", "0.002"))
AI = as.data.frame(lapply(AI, as.numeric))

plant.decay = cbind(AI, plant.table.sig$ele, plant.table.sig$slope, micro.table.sig$slope)
colnames(plant.decay) = c("AI", "Elevation", "Plant.decay", "Micro.decay")


lm.plant = lm(Plant.decay ~ AI, plant.decay)
qu.plant = lm(Plant.decay ~ AI + I(AI^2), plant.decay)
AIC(lm.plant, qu.plant)
summary(qu.plant)
# Adjusted R-squared:  0.3948  p-value: 0.07156


lm.Micro = lm(Micro.decay ~ AI, plant.decay)
qu.Micro  = lm(Micro.decay ~ AI + I(AI^2), plant.decay)
AIC(lm.Micro, qu.Micro)
summary(lm.Micro)
# Adjusted R-squared: 0.005652  p-value: 0.3352

plant.decay.me = melt(plant.decay, id = c("AI", "Elevation"))

plant.decay 

envi.qu.F1 = plant.decay.me$variable %in% c("Plant.decay")
envi.lm.F1 = plant.decay.me$variable %in% c("Micro.decay")

ggplot(plant.decay.me, aes(x=AI, y=value, color = AI))+ 
  geom_point(size = 2)+
  stat_smooth(method="lm",formula=y~x,size=0.8, lty=2, color="Black",
              se = FALSE, data = plant.decay.me[envi.lm.F1,])+
  # quadratic lm sig==T
  stat_smooth(method="lm",formula=y~x+I(x^2),size=0.8, lty=2, color="Black",
              se = FALSE, data = plant.decay.me[envi.qu.F1,])+
  # scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  facet_wrap(. ~ variable)+
  theme_bw()+  
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(strip.background = element_blank())+
  theme(axis.text.x  = element_text(angle=45,vjust=0.5))+
  theme(text = element_text(size = 12))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))+
  xlab("AI")+ylab("Decay rate")+
  scale_x_continuous(breaks=c(0.002,0.15,0.3,0.45,0.58,0.69,0.78,
                              0.86,0.94,0.99))

ggsave("../2.DX2013_145sites/figures/AI.lm.qu.plant.micro.decay.pdf",
       width = 6.5, height = 3.5)


plant.decay.me$AI = factor(plant.decay.me$AI)

ggplot(plant.decay.me, aes(x=AI, y=value, color = AI))+ 
  geom_point(size = 2)+
  stat_smooth(method="lm",formula=y~x,size=0.8, lty=2, color="Black",
              se = FALSE, data = plant.decay.me[envi.lm.F1,])+
  # quadratic lm sig==T
  stat_smooth(method="lm",formula=y~x+I(x^2),size=0.8, lty=2, color="Black",
              se = FALSE, data = plant.decay.me[envi.qu.F1,])+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  facet_wrap(. ~ variable)+
  theme_bw()+  
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(strip.background = element_blank())+
  theme(axis.text.x  = element_text(angle=45,vjust=0.5))+
  theme(text = element_text(size = 12))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))+
  xlab("Aridity")+ylab("Decay rate")

  
ggsave("../2.DX2013_145sites/figures/AI.plant.micro.decay.pdf",
         width = 6.5, height = 3.5)
