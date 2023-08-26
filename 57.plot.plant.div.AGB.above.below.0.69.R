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
A.above0.69 = subset(all.plant, Group == "Above")

lm.plant.div = lm(Veg.mass ~ Plant.rich, A.above0.69)
qu.plant.div = lm(Veg.mass ~ Plant.rich + I(Plant.rich^2), A.above0.69)
AIC(lm.plant.div, qu.plant.div)
# df      AIC
# lm.plant.div  3 927.1890
# qu.plant.div  4 926.5222

summary(qu.plant.div)
# Adjusted R-squared:  0.1884, p-value: 3.825e-05


### B.below0.69 ###
B.below0.69 = subset(all.plant, Group == "Below")

lm.plant.div = lm(Veg.mass ~ Plant.rich, B.below0.69)
qu.plant.div = lm(Veg.mass ~ Plant.rich + I(Plant.rich^2), B.below0.69)
AIC(lm.plant.div, qu.plant.div)
# df      AIC
# lm.plant.div  3 1153.961
# qu.plant.div  4 1155.903

summary(lm.plant.div)
# Adjusted R-squared:  0.06486, p-value: 0.006792

envi.qu.T1 = all.plant$Group %in% c("Above")
envi.lm.T1 = all.plant$Group %in% c("Below")



plant.div.AGB = 
  ggplot(all.plant, aes(x=Plant.rich, y=Veg.mass, color = Group))+ 
  
  geom_point(size = 2)+
  # (method="lm", color = "black")+
  # scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  
  stat_smooth(method="lm",formula=y~x, size=1,
              se = T, data = all.plant[envi.lm.T1,])+
  
  stat_smooth(method="lm",formula=y~x+I(x^2),size=1, 
              se = T, data = all.plant[envi.qu.T1,])+
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

### A.above0.69 ###
lm.plant.net = lm(Veg.mass ~ Plant.net, A.above0.69)
qu.plant.net = lm(Veg.mass ~ Plant.net + I(Plant.net^2), A.above0.69)
AIC(lm.plant.net, qu.plant.net)
# df      AIC
# lm.plant.net  3 943.3495
# qu.plant.net  4 931.8896

summary(qu.plant.net)
# Adjusted R-squared:  0.139, p-value: 0.0005126


### B.below0.69 ###
lm.plant.net = lm(Veg.mass ~ Plant.net, B.below0.69)
qu.plant.net = lm(Veg.mass ~ Plant.net + I(Plant.net^2), B.below0.69)
AIC(lm.plant.net, qu.plant.net)
# df      AIC
# lm.plant.net  3 1145.428
# qu.plant.net  4 1143.483

summary(qu.plant.net)
# Adjusted R-squared:  0.169, p-value: 6.189e-05

net.qu.T1 = all.plant$Group %in% c("Above", "Below")


plant.net.AGB = 
  ggplot(all.plant, aes(x=Plant.net, y=Veg.mass, color = Group))+ 
  
  geom_point(size = 2)+

  stat_smooth(method="lm",formula=y~x+I(x^2),size=1, 
              se = T, data = all.plant[net.qu.T1,])+
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


ggsave("../2.DX2013_145sites/final.figures2023/Figure.plant.rich.net.AGB.above.below0.69.pdf",
       width = 7.5, height = 3)
