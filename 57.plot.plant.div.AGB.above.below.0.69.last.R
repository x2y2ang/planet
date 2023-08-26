rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
getwd()

# library
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


all.plant = read.csv("../2.DX2013_145sites/data/plant1314.cov.network.mass.rich.for.SEM.mean.csv",
                     row.names = 1)



all.plant$pre.AI = as.factor(all.plant$pre.AI)


### A.above0.69 ###
lm.all.AGB = lm(Veg.mass ~ Plant.rich, all.plant)
summary(lm.all.AGB)
# Adjusted R-squared:  0.1538; p-value: 1.566e-08, slope = 13.9

A.above0.69 = subset(all.plant, Group == "Above")

lm.rich.A.div = lm(Veg.mass ~ Plant.rich, A.above0.69)
summary(lm.rich.A.div)
# Adjusted R-squared:  0.1736, p-value: 2.369e-05, slope = 10.1

### B.below0.69 ###
B.below0.69 = subset(all.plant, Group == "Below")

lm.B.plant.div = lm(Veg.mass ~ Plant.rich, B.below0.69)
summary(lm.B.plant.div)
# Adjusted R-squared:  0.06486, p-value: 0.006792, slope = 8.929


envi.lm.T1 = all.plant$Group %in% c("Above", "Below")


plant.div.AGB = 
  ggplot(all.plant, aes(x=Plant.rich, y=Veg.mass, color = Group))+ 
  
  geom_point(size = 2)+
  # (method="lm", color = "black")+
  # scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  
  stat_smooth(method="lm",formula=y~x, size=1,
              se = T, data = all.plant[envi.lm.T1,])+

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
  xlab("Plant.rich")+ylab("AGB")

plant.div.AGB


##########  plant.net and AGB ###########
lm.all.net = lm(Veg.mass ~ Plant.net, all.plant)
summary(lm.all.net)
# Adjusted R-squared:  0.0834; p-value: 3.457e-05, slope = 254.4

### A.above0.69 ###
lm.plant.A.net = lm(Veg.mass ~ Plant.net, A.above0.69)
summary(lm.plant.A.net)
# Adjusted R-squared:  0.01301, p-value: 0.1427, slope = 60.6

### B.below0.69 ###
lm.plant.B.net = lm(Veg.mass ~ Plant.net, B.below0.69)
summary(lm.plant.B.net)
# Adjusted R-squared:  0.1436, p-value: 7.654e-05, slope = 382.3

net.lm.T1 = all.plant$Group %in% c("Below")
net.lm.F1 = all.plant$Group %in% c("Above")


plant.net.AGB = 
  ggplot(all.plant, aes(x=Plant.net, y=Veg.mass, color = Group))+ 
  
  geom_point(size = 2)+
  
  stat_smooth(method="lm",formula=y~x,size=1, 
              se = T, data = all.plant[net.lm.T1,])+
  stat_smooth(method="lm",formula=y~x,size=1,lty=2, 
              se = T, data = all.plant[net.lm.F1,])+
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
  xlab("Plant.net")+ylab("AGB")

plant.net.AGB

library(cowplot)

plot_grid(plant.div.AGB, plant.net.AGB)


ggsave("../2.DX2013_145sites/final.figures2023/Figure.last.lm.plant.rich.net.AGB.above.below0.69.pdf",
       width = 7.5, height = 3)


#### whole graph ####
Rich.AGB = 
  ggplot(all.plant, aes(x=Plant.rich, y=Veg.mass))+ 
  geom_point(size = 2)+
  stat_smooth(method="lm",formula=y~x, size=1)+
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
  xlab("Plant.rich")+ylab("AGB")

Rich.AGB



Rich.net = 
  ggplot(all.plant, aes(x=Plant.net, y=Veg.mass))+ 
  geom_point(size = 2)+
  stat_smooth(method="lm",formula=y~x, size=1)+
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
  xlab("Plant.rich")+ylab("AGB")

Rich.net


plot_grid(Rich.AGB, Rich.net)

ggsave("../2.DX2013_145sites/final.figures2023/Figure.last.lm.plant.rich.net.AGB.all.pdf",
       width = 5, height = 3)
