
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
library(cowplot)



micro.soil.all = read.csv("../2.DX2013_145sites/data/soil.micro.decay.last.add.CV.Asy.pre.AI.csv")

plant.all = read.csv("../2.DX2013_145sites/data/plant.coverage.last.add.CV.Asy.pre.AI.csv")

colnames(micro.soil.all)

micro.soil.all$pre.AI = factor(micro.soil.all$pre.AI)
micro.soil.all$Elevation = factor(micro.soil.all$Elevation)

lm.micro.cv = lm(micro.soil.all$Micro.Rich.CV ~ micro.soil.all$Micro.decay, micro.soil.all)
summary(lm.micro.cv)
# Adjusted R-squared: 0.04545  p-value: 0.1052

qu.micro.cv = lm(micro.soil.all$Micro.Rich.CV ~ micro.soil.all$Micro.decay + 
                   I((micro.soil.all$Micro.decay)^2), micro.soil.all)

summary(qu.micro.cv)


AIC(lm.micro.cv, qu.micro.cv)
# df       AIC
# lm.micro.cv  3 -195.4009
# qu.micro.cv  4 -193.4009


p11 = ggplot(micro.soil.all, aes(Micro.decay, Micro.Rich.CV, color = pre.AI))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y~x, lty=2, color = "black")+
  scale_colour_manual(values = c("#1F78B4", "#43919B", "#6CC4A1", "#B2DF8A", "#DAE2B6",
                                 "#E5BA73", "#FDBF6F", "#FF7F00", "#FB9A99", "#E31A1C"))+
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
  xlab("Micro.decay")+ylab("Micro.Rich.CV")

p11

ggsave("../2.DX2013_145sites/final.figures2023/Figure.S6.Micro.decay.Micro.Rich.CV.legend.pdf",
       width = 6.5, height = 5)

p1 = ggplot(micro.soil.all, aes(Micro.decay, Micro.Rich.CV, color = pre.AI))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y~x, lty=2, color = "black")+
  scale_colour_manual(values = c("#1F78B4", "#43919B", "#6CC4A1", "#B2DF8A", "#DAE2B6",
                                 "#E5BA73", "#FDBF6F", "#FF7F00", "#FB9A99", "#E31A1C"))+
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
  xlab("Micro.decay")+ylab("Micro.Rich.CV")+
  theme(legend.position = 'none')

p1


### micro.Asy
colnames(micro.soil.all)

lm.micro.Asy = lm(micro.soil.all$Micro.Rich.Asy ~ micro.soil.all$Micro.decay, micro.soil.all)
summary(lm.micro.Asy)
# Adjusted R-squared:  0.1868 p-value: 0.00393

qu.micro.Asy = lm(micro.soil.all$Micro.Rich.Asy ~ micro.soil.all$Micro.decay + 
                    I((micro.soil.all$Micro.decay)^2), micro.soil.all)

summary(qu.micro.Asy)


AIC(lm.micro.Asy, qu.micro.Asy)
# df       AIC
# lm.micro.Asy  3 -8.471883
# qu.micro.Asy  4 -6.824824


p2 = ggplot(micro.soil.all, aes(Micro.decay, Micro.Rich.Asy, color = pre.AI))+
  geom_point(size = 2)+
  stat_smooth(method="lm", color = "black")+
  scale_colour_manual(values = c("#1F78B4", "#43919B", "#6CC4A1", "#B2DF8A", "#DAE2B6",
                                 "#E5BA73", "#FDBF6F", "#FF7F00", "#FB9A99", "#E31A1C"))+
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
  xlab("Micro.decay")+ylab("Micro.Rich.Asy")+
  theme(legend.position = 'none')

p2


### Soil.CV
colnames(micro.soil.all)

lm.Soil.CV = lm(micro.soil.all$Soil.PCA1.CV ~ micro.soil.all$Soil.decay, micro.soil.all)
summary(lm.Soil.CV)


qu.Soil.CV = lm(micro.soil.all$Soil.PCA1.CV ~ micro.soil.all$Soil.decay + 
                  I((micro.soil.all$Soil.decay)^2), micro.soil.all)

summary(qu.Soil.CV)
# Adjusted R-squared:  0.03857 p-value: 0.19

AIC(lm.Soil.CV, qu.Soil.CV)
# df      AIC
# lm.Soil.CV  3 86.36424
# qu.Soil.CV  4 85.82081


p3 = ggplot(micro.soil.all, aes(Soil.decay, Soil.PCA1.CV, color = pre.AI))+
  geom_point(size = 2)+
  stat_smooth(method="lm",formula = y~x+I(x^2), lty = 2, color = "black")+
  scale_colour_manual(values = c("#1F78B4", "#43919B", "#6CC4A1", "#B2DF8A", "#DAE2B6",
                                 "#E5BA73", "#FDBF6F", "#FF7F00", "#FB9A99", "#E31A1C"))+
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
  xlab("Soil.decay")+ylab("Soil.PCA1.CV")+
  theme(legend.position = 'none')

p3



### Soil.Asy
colnames(micro.soil.all)

lm.Soil.Asy = lm(micro.soil.all$Soil.PCA1.Asy ~ micro.soil.all$Soil.decay, micro.soil.all)
summary(lm.Soil.Asy)
# Adjusted R-squared:  0.306 p-value: 0.000188

qu.Soil.Asy = lm(micro.soil.all$Soil.PCA1.Asy ~ micro.soil.all$Soil.decay + 
                   I((micro.soil.all$Soil.decay)^2), micro.soil.all)

summary(qu.Soil.Asy)


AIC(lm.Soil.Asy, qu.Soil.Asy)
# df      AIC
# lm.Soil.Asy  3 21.63284
# qu.Soil.Asy  4 22.88815



p4 = ggplot(micro.soil.all, aes(Soil.decay, Soil.PCA1.Asy, color = pre.AI))+
  geom_point(size = 2)+
  stat_smooth(method="lm",formula = y~x, color = "black")+
  scale_colour_manual(values = c("#1F78B4", "#43919B", "#6CC4A1", "#B2DF8A", "#DAE2B6",
                                 "#E5BA73", "#FDBF6F", "#FF7F00", "#FB9A99", "#E31A1C"))+
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
  xlab("Soil.decay")+ylab("Soil.PCA1.Asy")+
  theme(legend.position = 'none')

p4



### Plant.CV
colnames(plant.all)

lm.Plant.CV = lm(plant.all$Plant.rich.CV ~ plant.all$Plant.cov, plant.all)
summary(lm.Plant.CV)


qu.Plant.CV = lm(plant.all$Plant.rich.CV ~ plant.all$Plant.cov + 
                   I((plant.all$Plant.cov)^2), plant.all)

summary(qu.Plant.CV)
# Adjusted R-squared: 0.03863  p-value: 0.1977

AIC(lm.Plant.CV, qu.Plant.CV)
# df       AIC
# lm.Plant.CV  3 -67.33058
# qu.Plant.CV  4 -67.47476

plant.all$pre.AI = as.factor(plant.all$pre.AI)

p5 = ggplot(plant.all, aes(Plant.cov, Plant.rich.CV, color = pre.AI))+
  geom_point(size = 2)+
  stat_smooth(method="lm",formula = y~x+I(x^2), lty = 2, color = "black")+
  scale_colour_manual(values = c("#1F78B4", "#43919B", "#6CC4A1", "#B2DF8A", "#DAE2B6",
                                 "#E5BA73", "#FDBF6F", "#FF7F00", "#FB9A99", "#E31A1C"))+
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
  xlab("Plant.cov")+ylab("Plant.rich.CV")+
  theme(legend.position = 'none')

p5


### Plant.Asy
colnames(plant.all)

lm.Plant.Asy = lm(plant.all$Plant.rich.Asy ~ plant.all$Plant.cov, plant.all)
summary(lm.Plant.Asy)
# Adjusted R-squared: -0.007954 p-value: 0.4008

qu.Plant.Asy = lm(plant.all$Plant.rich.Asy ~ plant.all$Plant.cov + 
                    I((plant.all$Plant.cov)^2), plant.all)

summary(qu.Plant.Asy)


AIC(lm.Plant.Asy, qu.Plant.Asy)
# df      AIC
# lm.Plant.Asy  3 38.17470
# qu.Plant.Asy  4 38.44224

plant.all$pre.AI = as.factor(plant.all$pre.AI)
plant.all$Ele = as.factor(plant.all$Ele)

p6 = ggplot(plant.all, aes(Plant.cov, Plant.rich.Asy, 
                           color = pre.AI))+
  geom_point(size = 2)+
  stat_smooth(method="lm",formula = y~x, lty = 2, color = "black")+
  scale_colour_manual(values = c("#1F78B4", "#43919B", "#6CC4A1", "#B2DF8A", "#DAE2B6",
                                 "#E5BA73", "#FDBF6F", "#FF7F00", "#FB9A99", "#E31A1C"))+
  scale_shape_manual(values = c(1:10))+
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
  xlab("Plant.cov")+ylab("Plant.rich.Asy")+
  theme(legend.position = 'none')

p6


#### no legend ###
plot_grid(p5,p1,p3,p6,p2,p4)

ggsave("../2.DX2013_145sites/final.figures2023/Figure.S6.1314year.plant.cover.micro.soil.decey.cor.CV.Asy.no.legend.pdf",
       width = 6.5, height = 5)

