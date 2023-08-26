
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

envi1314 = data.frame(envi1314)

cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

veg14 = read.csv("../2.DX2013_145sites/data/plant.2014.species.number.csv",
                 row.names = 1)

envi.mass = read.csv("../2.DX2013_145sites/data/envi.1314.add.mass.csv",
                     row.names = 1)

envi.mass.not.gram = read.csv("../2.DX2013_145sites/data/envi.1314.add.mass.csv",
                              row.names = 1)


#### 2014 plant quadrat areas 1 m2 #####

#### standard with 2020 0.25 m2 ########
veg14 = round((veg14/4))

# na value equal 0
veg14[is.na(veg14)] <- 0

source("../2.DX2013_145sites/function/diversity.functions.R")

veg14.div = quick.diversity(veg14)
veg13.div = quick.diversity(Veg.sp)

# combine 1314
veg.div1314 = rbind(veg13.div, veg14.div)

envi188 = envi1314[rownames(veg.div1314), ]
cli188 = cli1314.last[rownames(veg.div1314), ]


vec.tmp = cbind(cli188, veg.div1314)
plot(vec.tmp$Act.ele, vec.tmp$comm.richness)
colnames(vec.tmp)
vec.tmp$Elevation = factor(vec.tmp$Elevation)

veg.rich = subset(vec.tmp[, c("Elevation", "pre.AI", "comm.richness")], )

colnames(veg.rich) = c("Elevation", "pre.AI", "Veg.rich")

veg.rich.me <- melt(veg.rich, id=c("Elevation", "pre.AI"))

p2 = ggplot(veg.rich.me, aes(pre.AI, value, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", color = "black")+
  # scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1",  "#43919B", "#1F78B4", 
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4"))+
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
  xlab("Aridity")+ylab("Plant richness")

p2



### qu lm ###

### veg.rich.me ###
lm.plant = lm(value ~ pre.AI, veg.rich.me)
qu.plant = lm(value ~ pre.AI + I(pre.AI^2), veg.rich.me)
AIC(lm.plant, qu.plant)
# df      AIC
# lm.plant  3 860.3806
# qu.plant  4 862.1163

summary(qu.plant)

# Adjusted R-squared:  0.1137, p-value: 5.261e-06

ggplot(veg.rich.me, aes(pre.AI, value, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula= y ~ x + I(x^2), color = "black")+
  # scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1",  "#43919B", "#1F78B4", 
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4"))+
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
  xlab("Aridity")+ylab("Plant richness")

ggsave("../2.DX2013_145sites/final.figures2023/Figure.plant.rich.quadradic.fit.pdf",
       width = 3.5, height = 3)


# veg mass   
mass2014 = read.csv("../2.DX2013_145sites/data/last.DX2014.plant.biomass.total.species.number.csv")

mass2014 = mass2014[, -1]

colnames(mass2014)[1] = c("")

rownames(mass2014) = mass2014[,1]

mass1314 = c(Veg.div[, c("AGB")], mass2014[, c("Biomass")])

site = c(rownames(Veg.div), rownames(mass2014))

mass.all = cbind(site, mass1314)

colnames(mass.all)[1] = c("")

rownames(mass.all) = mass.all[,1]

ele.tmp = cli1314.last[rownames(mass.all),]

mass.all = data.frame(mass.all)
mass.num = round(as.numeric(mass.all$mass1314))


mass.all.ele = cbind(ele.tmp, mass.num)

mass.all.ele$Elevation = factor(mass.all.ele$Elevation)

veg.mass = subset(mass.all.ele[, c("Elevation", "pre.AI", "mass.num")], )

colnames(veg.mass) = c("Elevation", "pre.AI", "Veg.mass")

# write.csv(veg.mass, "../2.DX2013_145sites/data/veg1314.biomass.last.csv")

veg.mass.me <- melt(veg.mass, id=c("Elevation", "pre.AI"))

lm.veg.mass = lm(veg.mass.me$value ~ veg.mass.me$pre.AI, veg.mass.me)
summary(lm.veg.mass)

qu.veg.mass = lm(veg.mass.me$value ~ veg.mass.me$pre.AI + I((veg.mass.me$pre.AI)^2), veg.mass.me)
summary(qu.veg.mass)
#Adjusted R-squared: 0.2663 p-value: 1.008e-12

AIC(lm.veg.mass, qu.veg.mass)
# qu.veg.mass  4 2033.061

p1 = ggplot(veg.mass.me, aes(pre.AI, value, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula=y~x+I(x^2), 
              color = "black")+
  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1",  "#43919B", "#1F78B4", 
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4"))+
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
  xlab("Aridity")+ylab("Plant biomass")

p1


# micro richness 
micro.rich = quick.diversity(comm1314.samp)

micro.rich.AI = data.frame(cbind(cli1314.last$Elevation, cli1314.last$pre.AI, 
                                 micro.rich$comm.richness))  

colnames(micro.rich.AI) = c("Elevation", "pre.AI", "Micro.rich")

micro.rich.me <- melt(micro.rich.AI, id=c("Elevation", "pre.AI"))

lm.micro.rich = lm(micro.rich.me$value ~ micro.rich.me$pre.AI, micro.rich.me)
summary(lm.micro.rich)

qu.micro.rich = lm(micro.rich.me$value ~ micro.rich.me$pre.AI + I((micro.rich.me$pre.AI)^2), veg.mass.me)
summary(qu.micro.rich)
# Adjusted R-squared:  0.003688 p-value: 0.1785

AIC(lm.micro.rich, qu.micro.rich)
# lm.micro.rich  3 3078.418

micro.rich.me$Elevation = factor(micro.rich.me$Elevation)

gp = c(rep("A", 60), rep("B", 163))

micro.rich.me1 = cbind(gp, micro.rich.me)

A = subset(micro.rich.me1[1:60, ])

summary(lm(A$value ~ A$pre.AI, A))
# Adjusted R-squared:  0.08094, p-value: 0.01569

B = subset(micro.rich.me1[61:223, ])

summary(lm(B$value ~ B$pre.AI, B))
# p-value: 0.6091

lm.T1 = micro.rich.me1$gp %in% c("A")
lm.F1 = micro.rich.me1$gp %in% c("B")
# 
micro.rich.me$Elevation = factor(micro.rich.me$Elevation)

p3 = ggplot(micro.rich.me, aes(pre.AI, value, color = Elevation))+
  geom_point(size = 2)+
  # linear lm sig==T
  stat_smooth(aes(colour=gp),method="lm",formula=y~x,size=0.8,
              se = T,data = micro.rich.me1[lm.T1,])+
  # linear lm sig==F
  stat_smooth(aes(colour=gp),method="lm",formula=y~x,size=0.8, lty = 2,
              se = FALSE,data = micro.rich.me1[lm.F1,])+
  
  # stat_smooth(method="lm", formula=y~x, se = FALSE,  lty = 2,
  #             color = "black")+
  
  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1",  "#43919B", "#1F78B4", 
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4"))+
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
  xlab("Aridity")+ylab("Microbial richness")

p3


############# Micro PLFA ######################

PLFA.la = read.csv("../2.DX2013_145sites/data/micro.rich.inter.biomass.combine.csv", row.names = 1)

# cli13 = cli1314.last[rownames(envi), ]

micro.mass = data.frame(cli1314.last$Elevation, cli1314.last$pre.AI, PLFA.la$Micro.mass) 

rownames(micro.mass) = rownames(cli1314.last)

# write.csv(micro.mass, "../2.DX2013_145sites/data/micro.PLFA145sites.csv")

colnames(micro.mass) = c("Elevation", "pre.AI", "Micro.mass")

micro.mass.me <- melt(micro.mass, id=c("Elevation", "pre.AI"))


lm.micro.mass = lm(micro.mass.me$value ~ micro.mass.me$pre.AI, micro.mass.me)
summary(lm.micro.mass)

qu.micro.mass = lm(micro.mass.me$value ~ micro.mass.me$pre.AI + 
                     I((micro.mass.me$pre.AI)^2), micro.mass.me)
summary(qu.micro.mass)
#Adjusted R-squared:  0.5959 p-value: < 2.2e-16

AIC(lm.micro.mass, qu.micro.mass)
# qu.micro.mass  4 1367.305

micro.mass.me$Elevation = factor(micro.mass.me$Elevation)

p4 = ggplot(micro.mass.me, aes(pre.AI, value, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula=y~x+I(x^2), 
              color = "black")+
  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1",  "#43919B", "#1F78B4", 
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4"))+
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
  xlab("Aridity")+ylab("Microbial PLFA")

p4


############## EMF ################

EMF = cbind(envi1314[, c("TOC", "TN", "DOC", "DON", "TP")], 
            envi.mass[,c("Veg.mass")], PLFA.la$Micro.mass)

# z-score
EMF.tmp = scale(EMF)
# mean
EMF.la = data.frame(apply(EMF.tmp, 1, mean))
colnames(EMF.la) = c("EMF")

EMF.AI = data.frame(envi1314$Elevation, cli1314.last$pre.AI, EMF.la) 


colnames(EMF.AI) = c("Elevation", "pre.AI", "EMF")

EMF.AI.me <- melt(EMF.AI, id=c("Elevation", "pre.AI"))


lm.EMF.AI = lm(EMF.AI.me$value ~ EMF.AI.me$pre.AI, EMF.AI.me)
summary(lm.EMF.AI)

qu.EMF.AI = lm(EMF.AI.me$value ~ EMF.AI.me$pre.AI + I((EMF.AI.me$pre.AI)^2), veg.mass.me)
summary(qu.EMF.AI)
# Adjusted R-squared:  0.7132 p-value: < 2.2e-16

AIC(lm.EMF.AI, qu.EMF.AI)
# qu.EMF.AI  4 356.6067

EMF.AI.me$Elevation = factor(EMF.AI.me$Elevation)

p5 = ggplot(EMF.AI.me, aes(pre.AI, value, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula=y~x+I(x^2), 
              color = "black")+
  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1",  "#43919B", "#1F78B4", 
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4"))+
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
  xlab("Aridity")+ylab("EMF")

p5


############## soil pca1 #################
soil = envi1314[, c("TOC", "TN", "DOC", "DON", "TP")] 

soil.pca = PCA(soil, scale.unit=T, graph=F)
summary(soil.pca)
# Dim.1 66.681
soil.pca.dd <- data.frame(soil.pca$ind$coord)


Soil.pca1.AI = data.frame(envi1314$Elevation, cli1314.last$pre.AI, soil.pca.dd$Dim.1) 


colnames(Soil.pca1.AI) = c("Elevation", "pre.AI", "Soil.pca1")

Soil.pca1.AI.me <- melt(Soil.pca1.AI, id=c("Elevation", "pre.AI"))


lm.Soil.pca1.AI = lm(Soil.pca1.AI.me$value ~ Soil.pca1.AI.me$pre.AI, Soil.pca1.AI.me)
summary(lm.Soil.pca1.AI)

qu.Soil.pca1.AI = lm(Soil.pca1.AI.me$value ~ Soil.pca1.AI.me$pre.AI + I((Soil.pca1.AI.me$pre.AI)^2), veg.mass.me)
summary(qu.Soil.pca1.AI)
# Adjusted R-squared:   0.578  p-value: < 2.2e-16

AIC(lm.Soil.pca1.AI, qu.Soil.pca1.AI)
# qu.Soil.pca1.AI  4 714.9553

Soil.pca1.AI.me$Elevation = factor(Soil.pca1.AI.me$Elevation)

p6 = ggplot(Soil.pca1.AI.me, aes(pre.AI, value, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula=y~x+I(x^2), 
              color = "black")+
  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1",  "#43919B", "#1F78B4", 
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4",
                                 "#FF7F00", "#1F78B4"))+
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
  xlab("Aridity")+ylab("Soil.pca1")

p6


library(cowplot)
plot_grid(p1,p4,p5,p2,p3,p6)

ggsave("../2.DX2013_145sites/final.figures2023/Figure.S2.plot.plant.micro.segment.soil.EMF.AI.patterns.pdf",
       width = 10, height = 5)

