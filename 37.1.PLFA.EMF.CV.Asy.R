

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
library(FactoMineR)


# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")

load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")

cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

veg14 = read.csv("../2.DX2013_145sites/data/plant.2014.species.number.csv",
                 row.names = 1)

envi.mass = read.csv("../2.DX2013_145sites/data/envi.1314.add.mass.csv",
                     row.names = 1)


############# Micro PLFA ######################

micro.mass = data.frame(envi.mass$Elevation, cli1314.last$pre.AI, envi.mass$PLFA) 

rownames(micro.mass) = rownames(envi.mass)

colnames(micro.mass) = c("Elevation", "pre.AI", "Micro.mass")

micro.mass.me <- melt(micro.mass, id=c("Elevation", "pre.AI"))

vec.mean = dcast(micro.mass.me, Elevation ~ variable, mean)

# melt mean 
vec.mean.me = melt(vec.mean, id = c("Elevation"))

# calculate standard
vec.sd = dcast(micro.mass.me, Elevation ~ variable, sd)

# melt sd
vec.sd.me = melt(vec.sd, id = c("Elevation"))


# number of each elevation
vec.len <- aggregate(micro.mass.me$value, by=list(micro.mass.me$Elevation, 
                                                  micro.mass.me$variable), 
                     FUN=length)



# asymmetry as (max - mean)/(mean - min)
vec.max = aggregate(value ~ Elevation, data = micro.mass.me, FUN=max)

vec.min = aggregate(value ~ Elevation, data = micro.mass.me, FUN=min)


# mean, se, CV, asymmetry for vec
vec.me = melt(vec.tmp, id = c("Elevation"))

colnames(vec.me)

vec.mean = dcast(vec.me, Elevation ~ variable, mean)

# melt mean 
vec.mean.me = melt(vec.mean, id = c("Elevation"))

# calculate standard
vec.sd = dcast(vec.me, Elevation ~ variable, sd)

# melt sd
vec.sd.me = melt(vec.sd, id = c("Elevation"))


# number of each elevation
vec.len <- aggregate(vec.me$value, by=list(vec.me$Elevation, vec.me$variable), 
                     FUN=length)



# asymmetry as (max - mean)/(mean - min)
vec.max = aggregate(value ~ Elevation, data = vec.me, FUN=max)

vec.min = aggregate(value ~ Elevation, data = vec.me, FUN=min)


# mean, se, CV, asymmetry for vec
vec.all <- data.frame(vec.mean.me, se = (vec.sd.me$value)/sqrt(vec.len$x),
                      CV = ((vec.sd.me$value)/(vec.mean.me$value)),
                      Asy = ((vec.max$value) - (vec.mean.me$value))/((vec.mean.me$value) - (vec.min$value)))

write.csv(vec.all, "../2.DX2013_145sites/data/Micro.PLFA.CV.Asy.csv")


############## EMF ################

EMF = cbind(envi1314[, c("TOC", "TN", "DOC", "DON", "TP")], 
            envi.mass[,c("Veg.mass", "PLFA")])
# z-score
EMF.tmp = scale(EMF)
# mean
EMF.la = data.frame(apply(EMF.tmp, 1, mean))
colnames(EMF.la) = c("EMF")

EMF.AI = data.frame(envi1314$Elevation,  EMF.la) 

colnames(EMF.AI) = c("Elevation", "EMF")


EMF.AI.me <- melt(EMF.AI, id=c("Elevation"))


vec.mean1 = dcast(EMF.AI.me, Elevation ~ variable, mean)

# melt mean 
vec.mean.me1 = melt(vec.mean1, id = c("Elevation"))

# calculate standard
vec.sd1 = dcast(vec.mean.me1, Elevation ~ variable, sd)

# melt sd
vec.sd.me1 = melt(vec.sd1, id = c("Elevation"))


# number of each elevation
vec.len1 <- aggregate(EMF.AI.me$value, by=list(EMF.AI.me$Elevation, 
                                              EMF.AI.me$variable), 
                     FUN=length)



# asymmetry as (max - mean)/(mean - min)
vec.max1 = aggregate(value ~ Elevation, data = EMF.AI.me, FUN=max)

vec.min1 = aggregate(value ~ Elevation, data = EMF.AI.me, FUN=min)


# mean, se, CV, asymmetry for vec
vec.me1 = melt(EMF.AI, id = c("Elevation"))

colnames(vec.me1)

vec.mean1 = dcast(vec.me1, Elevation ~ variable, mean)

# melt mean 
vec.mean.me1 = melt(vec.mean1, id = c("Elevation"))

# calculate standard
vec.sd1 = dcast(vec.me1, Elevation ~ variable, sd)

# melt sd
vec.sd.me1 = melt(vec.sd1, id = c("Elevation"))


# number of each elevation
vec.len1 <- aggregate(vec.me1$value, by=list(vec.me1$Elevation, vec.me1$variable), 
                     FUN=length)



# asymmetry as (max - mean)/(mean - min)
vec.max1 = aggregate(value ~ Elevation, data = vec.me1, FUN=max)

vec.min1 = aggregate(value ~ Elevation, data = vec.me1, FUN=min)


# mean, se, CV, asymmetry for vec
vec.all1 <- data.frame(vec.mean.me1, se = (vec.sd.me1$value)/sqrt(vec.len1$x),
                      CV = ((vec.sd.me1$value)/(vec.mean.me1$value)),
                      Asy = ((vec.max1$value) - (vec.mean.me1$value))/((vec.mean.me1$value) - (vec.min1$value)))


write.csv(vec.all1, "../2.DX2013_145sites/data/EMF.CV.Asy.csv")


######### plot  ############ 
EMF.PLFA.cv = rbind(vec.all, vec.all1)


pre.AI = cli1314.last[, c("pre.AI", "Elevation")]

vec.me2 = melt(pre.AI, id = c("Elevation"))

colnames(vec.me2)

vec.mean2 = dcast(vec.me2, Elevation ~ variable, mean)

pre.AI2 = rbind(vec.mean2, vec.mean2)

AI.CV.LA = cbind(pre.AI2, EMF.PLFA.cv)
AI.CV.LA1 = AI.CV.LA[,-1]

AI.CV.LA1$Elevation = factor(AI.CV.LA1$Elevation)

ggplot(AI.CV.LA1, aes(x=pre.AI, y=CV, color = Elevation))+
  geom_point(size = 3)+
  facet_wrap( ~ variable, scales = "free_y", ncol=2)+
  scale_colour_manual(values = brewer.pal(12,"Paired"))+
  xlab("Aridity") + ylab("")+
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
        axis.title.y=element_text(size=12))

ggsave("../2.DX2013_145sites/figures/EMF.PLFA.CV.Asymmetic.ele.pdf",
       height = 3,width = 6)


ggplot(AI.CV.LA1, aes(x=pre.AI, y=Asy, color = Elevation))+
  geom_point(size = 3)+
  facet_wrap( ~ variable, scales = "free_y", ncol=2)+
  scale_colour_manual(values = brewer.pal(12,"Paired"))+
  xlab("Aridity") + ylab("Asy")+
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
  

ggsave("../2.DX2013_145sites/figures/EMF.PLFA.Asymmetic.ele.pdf",
       height = 3,width = 6)
