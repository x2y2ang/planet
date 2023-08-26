
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# library
library(vegan)
library(lavaan)
library(semPlot)
library(AICcmodavg)
library(data.table)
library(MASS)
library(ggplot2)

# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")

load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")

all.micro = read.csv("../2.DX2013_145sites/data/Micro.rich.inter.biomass.combine.AI.group0.69.csv", 
                     row.names = 1)

micro.coh = read.csv("../2.DX2013_145sites/data/micro.cohesion.last1314.csv")


envi.mass.not.gram = read.csv("../2.DX2013_145sites/data/envi.1314.add.mass.27.7.2022.correct.csv",
                              row.names = 1)

cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

micro.mass = data.frame(envi.mass.not.gram$Elevation, cli1314.last$pre.AI, 
                        envi.mass.not.gram$PLFA, micro.coh$total.coh) 

colnames(micro.mass) = c("Elevation", "pre.AI", "PLFA", "Micro.net")

rownames(micro.mass)= rownames(cli1314.last)

source("../2.DX2013_145sites/function/diversity.functions.R")

dim(comm1314.samp)

micro1314 = quick.diversity(comm1314.samp)

net.plfa = cbind(micro.mass, micro1314$comm.richness)
colnames(net.plfa)[5] = c("Micro.rich")


PLFA = read.csv("../2.DX2013_145sites/data/PLFA.origin.csv", row.names = 1)


net.plfa2 = net.plfa[rownames(PLFA),]

net.plfa2$Group = ifelse(net.plfa2$pre.AI > 0.69, "Above","Below")


### A.above0.69 ###
A.above0.69 = subset(net.plfa2, Group == "Above")

lm.micro.div = lm(PLFA ~ Micro.rich, A.above0.69)
qu.micro.div = lm(PLFA ~ Micro.rich + I(Micro.rich^2), A.above0.69)
AIC(lm.micro.div, qu.micro.div)
# df     AIC
# lm.micro.div  3 133.077
# qu.micro.div  4 135.023

summary(lm.micro.div)
# Adjusted R-squared:  0.1682, p-value: 0.07213


### B.below0.69 ###
B.below0.69 = subset(net.plfa2, Group == "Below")

lm.micro.div = lm(PLFA ~ Micro.rich, B.below0.69)
qu.micro.div = lm(PLFA ~ Micro.rich + I(Micro.rich^2), B.below0.69)
AIC(lm.micro.div, qu.micro.div)
# df      AIC
# lm.micro.div  3 148.1758
# qu.micro.div  4 149.9380

summary(lm.micro.div)
# Adjusted R-squared:  -0.05937, p-value: 0.6502


envi.lm.F1 = net.plfa2$Group %in% c("Above", "Below")


micro.div.PLFA = 
  ggplot(net.plfa2, aes(x=Micro.rich, y=PLFA, color = Group))+ 
  
  geom_point(size = 2)+

  stat_smooth(method="lm",formula=y~x, size=1, lty =2,
              se = T, data = net.plfa2[envi.lm.F1,])+
  
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
  xlab("Micro.rich")+ylab("PLFA")

micro.div.PLFA


#####################################################
### A.above0.69 ###
lm.micro.net = lm(PLFA ~ Micro.net, A.above0.69)
qu.micro.net = lm(PLFA ~ Micro.net + I(Micro.net^2), A.above0.69)
AIC(lm.micro.net, qu.micro.net)
# df      AIC
# lm.micro.net  3 135.1679
# qu.micro.net  4 136.2992

summary(lm.micro.net)
# Adjusted R-squared:  0.0438, p-value: 0.2225


### B.below0.69 ###
lm.micro.net = lm(PLFA ~ Micro.net, B.below0.69)
qu.micro.net = lm(PLFA ~ Micro.net + I(Micro.net^2), B.below0.69)
AIC(lm.micro.net, qu.micro.net)
# df      AIC
# lm.micro.net  3 146.9158
# qu.micro.net  4 145.7585

summary(lm.micro.net)
# Adjusted R-squared:  0.02599, p-value: 0.2622


envi.lm.F1 = net.plfa2$Group %in% c("Above", "Below")



micro.net.PLFA = 
  ggplot(net.plfa2, aes(x=Micro.net, y=PLFA, color = Group))+ 
  
  geom_point(size = 2)+
  
  stat_smooth(method="lm",formula=y~x, size=1, lty =2,
              se = T, data = net.plfa2[envi.lm.F1,])+
  # stat_smooth(method="lm",formula=y~x+I(x^2), size=1, lty =2,
  #             se = T, data = net.plfa2[envi.qu.F1,])+
  
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
  xlab("Micro.net")+ylab("PLFA")

micro.net.PLFA

library(cowplot)

plot_grid(micro.div.PLFA, micro.net.PLFA)


ggsave("../2.DX2013_145sites/final.figures2023/Figure.Micro.rich.net.PLFA.30samples.above.below0.69.pdf",
       width = 7.5, height = 3)



