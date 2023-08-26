
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
cli188 = cli1314.last[rownames(veg.div1314), ]

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

me.micro.soil = melt(vec.tmp, id = c("Elevation"))

prec.tmp = c(round(cli1314.last[, c("pre.bio12")],0))

prec.tmp1 = rep(prec.tmp, 2)

me.micro.soil1 = cbind(prec.tmp1, me.micro.soil)

me.micro.soil1$Elevation = factor(me.micro.soil1$Elevation)

ggplot(me.micro.soil1, aes(x=prec.tmp1, y=value, color = Elevation))+ 
  geom_point(size = 3)+
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
  scale_colour_manual(values = brewer.pal(12,"Paired"))+
  xlab("Precipitation (mm)") + ylab("")

?geom_boxplot
ggplot(me.micro.soil1, aes(prec.tmp1, value, color = Elevation))+
  stat_boxplot(geom="errorbar", 
               width = 0.3, size=0.1, color="#7F7C82",
               position=position_dodge(0.6))+ # aes(fill=Group), lty =2,  color=Group
  geom_boxplot(position=position_dodge(0.6), 
               # lty size
               size=0.2, 
               width=0.6,
               #lty = 2,
               outlier.shape = 19,
               outlier.size = 0.2,
               outlier.stroke = 0.2,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.1)+
  #scale_fill_manual(values = c("#FF1818","#11468F", "#FABB51"))+
  facet_wrap( ~ variable, scales ="free_y", ncol = 3)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))+
  xlab("")+ylab("")
