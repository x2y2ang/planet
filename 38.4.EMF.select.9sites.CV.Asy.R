

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


############# EMF ################

EMF = cbind(envi1314[, c("TOC", "TN", "DOC", "DON", "TP")], 
            envi.mass[,c("Veg.mass", "PLFA")])
# z-score
EMF.tmp = scale(EMF)
# mean
EMF.la = data.frame(apply(EMF.tmp, 1, mean))
colnames(EMF.la) = c("EMF")

EMF.AI = data.frame(envi1314$Elevation,  EMF.la) 

colnames(EMF.AI) = c("Elevation", "EMF")


vec.tmp = cbind(cli1314.last, EMF.la, envi1314)

EMF.tmp = data.frame(vec.tmp[, "EMF"])
Act.ele = data.frame(vec.tmp[, c("Act.ele")])
pre.bio12 = data.frame(vec.tmp[, c("pre.bio12")])
pre.AI = data.frame(vec.tmp[, c("pre.AI")])


ele = unique(vec.tmp[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  plant.t <- EMF.tmp[which(vec.tmp$Elevation == ele.tmp),]
  
  plant.com <- t(combn(plant.t,9))
  plant.com <- data.frame(plant.com)
  rownames(plant.com)
  
  plant.mean = apply(plant.com, 1, mean)
  plant.sd = apply(plant.com, 1, sd)
  
  CV.tmp= plant.sd/plant.mean
  
  # Asy strong functions for row
  plant.max = apply(plant.com, 1, max)
  plant.min = apply(plant.com, 1, min)
  
  Asy = ((plant.max - plant.mean)/(plant.mean - plant.min))
  
  
  ######## for mean ele #######
  act.ele.t <- Act.ele[which(vec.tmp$Elevation == ele.tmp),]
  
  act.ele.com <- t(combn(act.ele.t,9))
  act.ele.com <- data.frame(act.ele.com)
  rownames(act.ele.com)
  
  act.ele.mean = apply(act.ele.com, 1, mean)
  
  
  ####### for mean pre.bio12 ########
  pre.bio12.t <- pre.bio12[which(vec.tmp$Elevation == ele.tmp),]
  
  pre.bio12.com <- t(combn(pre.bio12.t,9))
  pre.bio12.com <- data.frame(pre.bio12.com)
  rownames(pre.bio12.com)
  
  pre.bio12.mean = apply(pre.bio12.com, 1, mean)
  
  
  ####### for mean pre.AI ########
  pre.AI.t <- pre.AI[which(vec.tmp$Elevation == ele.tmp),]
  
  pre.AI.com <- t(combn(pre.AI.t,9))
  pre.AI.com <- data.frame(pre.AI.com)
  rownames(pre.AI.com)
  
  pre.AI.mean = apply(pre.AI.com, 1, mean)
  
  
  ##########################################
  act.ele.CV.tmp = cbind(act.ele.mean, pre.bio12.mean, pre.AI.mean, 
                         CV.tmp, Asy)
  
  colnames(act.ele.CV.tmp) = c("Ele", "pre.bio12", "pre.AI", 
                               "plant.CV", "Asy")
  
  if (i==1)
  {
    act.ele.CV = act.ele.CV.tmp
  }
  else{
    act.ele.CV = rbind(act.ele.CV, act.ele.CV.tmp)
  }
}

act.ele.CV = data.frame(act.ele.CV)

ele.r = round(act.ele.CV$Ele, 0)

colnames(act.ele.CV)
CV.Asy = cbind(round(act.ele.CV$Ele, 0), round(act.ele.CV$pre.bio12, 0),
               round(act.ele.CV$pre.AI, 3), act.ele.CV[, c("plant.CV", "Asy")])

CV.Asy[1:50, 1:5]

colnames(CV.Asy) = c("Ele", "pre.bio12", "pre.AI", 
                     "EMF.CV", "EMF.Asy")

library(reshape2)
CV.Asy.md <- melt(CV.Asy, id=c("Ele", "pre.bio12", "pre.AI"))
CV.Asy.md[1:10, 1:5]


CV.me = dcast(CV.Asy.md, Ele  ~ variable, mean)


CV.me.sd = dcast(CV.Asy.md, Ele  ~ variable, sd)


# number of each elevation
hori.len <- aggregate(CV.Asy.md$value, by=list(CV.Asy.md$Ele,  CV.Asy.md$variable), 
                      FUN=length)


CV.me1 = melt(CV.me, id = c("Ele"))


# mean and se for hori
hori.all <- data.frame(CV.me1, se = (CV.me1$value)/sqrt(hori.len$x))

hori.all1 = na.omit(hori.all)

write.csv(hori.all, "../2.DX2013_145sites/data/EMF.1314.CV.Asy.elevation.csv")


ggplot(hori.all, aes(x=Ele, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.4)+
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))

ggsave("../2.DX2013_145sites/figures/EMF1314.CV.Asy.9site.pdf",
       width = 6.5, height = 3)


################# prec  ############### 

CV.Asy.md <- melt(CV.Asy, id=c("Ele", "pre.bio12", "pre.AI"))

prec.CV.me = dcast(CV.Asy.md, pre.bio12  ~ variable, mean)

prec.CV.me.sd = dcast(CV.Asy.md, pre.bio12  ~ variable, sd)


# number of each elevation
memory.limit(size = 35000)

prec.hori.len <- aggregate(CV.Asy.md$value, by=list(CV.Asy.md$pre.bio12,  CV.Asy.md$variable), 
                           FUN=length)


prec.CV.me1 = melt(prec.CV.me, id = c("pre.bio12"))



# mean and se for hori
prec.hori.all <- data.frame(prec.CV.me1, se = (prec.CV.me1$value)/sqrt(prec.hori.len$x))


write.csv(prec.hori.all, "../2.DX2013_145sites/data/EMF.1314.CV.Asy.precipitation.csv")


ggplot(prec.hori.all, aes(x=pre.bio12, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.4)+
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))

ggsave("../2.DX2013_145sites/figures/prec.EMF1314.rich.CV.9site.pdf",
       width = 6.5, height = 3)


colnames(CV.Asy)
AI.CV.Asy2 = CV.Asy[, c("pre.AI", "EMF.CV", "EMF.Asy")]

AI.CV.Asy3 = melt(AI.CV.Asy2, id = ("pre.AI")) 

memory.limit(size = 35000)
AI.CV.me = dcast(AI.CV.Asy3, pre.AI ~ variable, mean)

AI.CV.me.sd = dcast(AI.CV.Asy3, pre.AI ~ variable, sd)


# number of each elevation
AI.hori.len <- aggregate(AI.CV.Asy3$value, by=list(AI.CV.Asy3$pre.AI, 
                                                   AI.CV.Asy3$variable), FUN=length)

AI.CV.me1 = melt(AI.CV.me, id = c("pre.AI"))

# mean and se for hori
AI.hori.all <- data.frame(AI.CV.me1, se = (AI.CV.me1$value)/sqrt(AI.hori.len$x))

write.csv(AI.hori.all, "../2.DX2013_145sites/data/EMF1314.CV.Asy.aridity.csv")


ggplot(AI.hori.all, aes(x=pre.AI, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.5)+
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))

ggsave("../2.DX2013_145sites/figures/Aridity.EMF1314.CV.Asy.9site.pdf",
       width = 6.5, height = 3)

