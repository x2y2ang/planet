
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


# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")

load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")


cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

coh.all = read.csv("../2.DX2013_145sites/data/plant1314.188sites.coverage.cohesion.all.csv",
                   row.names = 1)


act = cli1314.last[rownames(coh.all), ]
coh.all = cbind(act[,1], coh.all)
colnames(coh.all)[1] = c("Act.ele")

# delete 4300 4400
vec.tmp.de = subset(coh.all, Elevation != 4300)

vec43 = subset(coh.all, Elevation == 4300)

vec43.CV = sd(vec43$Total.coh)/mean(vec43$Total.coh)
# 0.4146605

vec43.max = max(vec43$Total.coh)
vec43.min = min(vec43$Total.coh)
vec43.mean = mean(vec43$Total.coh)

vec43.Asy = ((vec43.max - vec43.mean)/(vec43.mean - vec43.min))
# 1.215776

#### loop 
net.tmp = data.frame(vec.tmp.de[, "Total.coh"])
Act.ele = data.frame(vec.tmp.de[, c("Act.ele")])
pre.AI = data.frame(vec.tmp.de[, c("pre.AI")])


ele = unique(vec.tmp.de[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  net.t <- net.tmp[which(vec.tmp.de$Elevation == ele.tmp),]
  
  net.com <- t(combn(net.t,9))
  net.com <- data.frame(net.com)
  rownames(net.com)
  
  plant.mean = apply(net.com, 1, mean)
  plant.sd = apply(net.com, 1, sd)
  
  CV.tmp= plant.sd/plant.mean
  
  # Asy strong functions for row
  plant.max = apply(net.com, 1, max)
  plant.min = apply(net.com, 1, min)
  
  Asy = ((plant.max - plant.mean)/(plant.mean - plant.min))
  
  
  ######## for mean ele #######
  act.ele.t <- Act.ele[which(vec.tmp.de$Elevation == ele.tmp),]
  
  act.ele.com <- t(combn(act.ele.t,9))
  act.ele.com <- data.frame(act.ele.com)
  rownames(act.ele.com)
  
  act.ele.mean = apply(act.ele.com, 1, mean)
  
  
  ####### for mean pre.AI ########
  pre.AI.t <- pre.AI[which(vec.tmp.de$Elevation == ele.tmp),]
  
  pre.AI.com <- t(combn(pre.AI.t,9))
  pre.AI.com <- data.frame(pre.AI.com)
  rownames(pre.AI.com)
  
  pre.AI.mean = apply(pre.AI.com, 1, mean)
  
  
  ##########################################
  act.ele.CV.tmp = cbind(act.ele.mean,  pre.AI.mean, 
                         CV.tmp, Asy)
  
  colnames(act.ele.CV.tmp) = c("Ele", "pre.AI", 
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

act.ele.CV[1:50, 1:4]

act.ele.CV = data.frame(act.ele.CV)

ele.r = round(act.ele.CV$Ele, 0)

colnames(act.ele.CV)
CV.Asy = cbind(round(act.ele.CV$Ele, 0),
               round(act.ele.CV$pre.AI, 3), act.ele.CV[, c("plant.CV", "Asy")])

CV.Asy[1:50, 1:4]

colnames(CV.Asy) = c("Ele",  "pre.AI", 
                     "plant.net.CV", "plant.net.Asy")

library(reshape2)
CV.Asy.md <- melt(CV.Asy, id=c("Ele", "pre.AI"))
CV.Asy.md[1:10, 1:4]


CV.me = dcast(CV.Asy.md, Ele  ~ variable, mean)


CV.me.sd = dcast(CV.Asy.md, Ele  ~ variable, sd)


# number of each elevation
hori.len <- aggregate(CV.Asy.md$value, by=list(CV.Asy.md$Ele,  
                                               CV.Asy.md$variable), 
                      FUN=length)


CV.me1 = melt(CV.me, id = c("Ele"))


# mean and se for hori
hori.all <- data.frame(CV.me1, se = (CV.me1$value)/sqrt(hori.len$x))

hori.all1 = na.omit(hori.all)

write.csv(hori.all, "../2.DX2013_145sites/data/Total.cohesion.plant.1314.CV.Asy.elevation.csv")

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

ggsave("../2.DX2013_145sites/figures/Total.cohesion.plant.1314.CV.Asy.pdf",
       width = 6.5, height = 3)


AI.CV.Asy2 = CV.Asy[, c("pre.AI", "plant.net.CV", "plant.net.Asy")]

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

write.csv(AI.hori.all, "../2.DX2013_145sites/data/Total.cohesion.plant.1314.CV.Asy.aridity.csv")


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

ggsave("../2.DX2013_145sites/figures/Total.cohesion.plant.1314.CV.Asy.aridity.pdf",
       width = 6.5, height = 3)

