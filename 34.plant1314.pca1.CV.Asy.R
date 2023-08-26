
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


# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")


load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")

cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

veg14 = read.csv("../2.DX2013_145sites/data/plant.2014.species.number.csv",
                 row.names = 1)

#### 2014 plant quadrat areas 1 m2 #####
#### standard with 2020 0.25 m2 ########
veg14 = round((veg14/4))

# na value equal 0
veg14[is.na(veg14)] <- 0

# to match the data 2014
sp36.45 = matrix(data=0, nrow = 110, ncol = 10, byrow = FALSE, dimnames = NULL)
colnames(sp36.45) = c("sp36", "sp37", "sp38", "sp39","sp40",
                      "sp41", "sp42", "sp43", "sp44","sp45")

veg13 = cbind(Veg.sp, sp36.45)


# combine veg data 2013 and 2014
veg1314 = rbind(veg13, veg14)


comm.hel.B = decostand(veg1314, "hellinger")
# bray
euc <- vegdist(comm.hel.B, method="bray")
# cmdscale
b.pcoa1 <- cmdscale(euc,k=(nrow(comm.hel.B)-1),eig=TRUE)
b.pcoa1
summary(b.pcoa1)
# The first two principal coordinates of each point
plant.pcoa <- as.data.frame(b.pcoa1$points[,1:2])
# change names
names(plant.pcoa) <- c("PCoA1","PCoA2")


############ 188 sites for plant data ##########
cli1314.tmp = cli1314.last[rownames(plant.pcoa),]
envi1314.tmp = envi1314[rownames(plant.pcoa),]


################    loop for 10 eles    ######################
vec.tmp = cbind(plant.pcoa, cli1314.tmp, envi1314.tmp)


# delete 4300  
vec.tmp.de = subset(vec.tmp, Elevation != 4300)

plant.pcoa1.43 = subset(vec.tmp, Elevation == 4300)

write.csv(plant.pcoa1.43, "../2.DX2013_145sites/data/plant.pcoa1.4300.csv")


plant.PCoA1.tmp = data.frame(vec.tmp.de[, "PCoA1"])
Act.ele = data.frame(vec.tmp.de[, c("Act.ele")])
pre.bio12 = data.frame(vec.tmp.de[, c("pre.bio12")])
pre.AI = data.frame(vec.tmp.de[, c("pre.AI")])


ele = unique(vec.tmp.de[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  PCoA1.t <- plant.PCoA1.tmp[which(vec.tmp.de$Elevation == ele.tmp),]
  
  PCoA1.com <- t(combn(PCoA1.t,9))
  PCoA1.com <- data.frame(PCoA1.com)
  rownames(PCoA1.com)
  
  PCoA1.mean = apply(PCoA1.com, 1, mean)
  PCoA1.sd = apply(PCoA1.com, 1, sd)
  
  CV.tmp= PCoA1.sd/PCoA1.mean
  
  # Asy strong functions for row
  PCoA1.max = apply(PCoA1.com, 1, max)
  PCoA1.min = apply(PCoA1.com, 1, min)
  
  Asy = ((PCoA1.max - PCoA1.mean)/(PCoA1.mean - PCoA1.min))
  
  
  ######## for mean ele #######
  act.ele.t <- Act.ele[which(vec.tmp.de$Elevation == ele.tmp),]
  
  act.ele.com <- t(combn(act.ele.t,9))
  act.ele.com <- data.frame(act.ele.com)
  rownames(act.ele.com)
  
  act.ele.mean = apply(act.ele.com, 1, mean)
  
  
  ####### for mean pre.bio12 ########
  pre.bio12.t <- pre.bio12[which(vec.tmp.de$Elevation == ele.tmp),]
  
  pre.bio12.com <- t(combn(pre.bio12.t,9))
  pre.bio12.com <- data.frame(pre.bio12.com)
  rownames(pre.bio12.com)
  
  pre.bio12.mean = apply(pre.bio12.com, 1, mean)
  
  
  ####### for mean pre.AI ########
  pre.AI.t <- pre.AI[which(vec.tmp.de$Elevation == ele.tmp),]
  
  pre.AI.com <- t(combn(pre.AI.t,9))
  pre.AI.com <- data.frame(pre.AI.com)
  rownames(pre.AI.com)
  
  pre.AI.mean = apply(pre.AI.com, 1, mean)
  
  
  ##########################################
  act.ele.CV.tmp = cbind(act.ele.mean, pre.bio12.mean, pre.AI.mean, 
                         CV.tmp, Asy)
  
  colnames(act.ele.CV.tmp) = c("Ele", "pre.bio12", "pre.AI", 
                               "PCoA1.CV", "Asy")
  
  if (i==1)
  {
    act.ele.CV = act.ele.CV.tmp
  }
  else{
    act.ele.CV = rbind(act.ele.CV, act.ele.CV.tmp)
  }
}

act.ele.CV = data.frame(act.ele.CV)

act.ele.CV[1:50, 1:5]

ele.r = round(act.ele.CV$Ele, 0)

CV.Asy = cbind(round(act.ele.CV$Ele, 0), round(act.ele.CV$pre.bio12, 0),
               round(act.ele.CV$pre.AI, 3), act.ele.CV[, c("PCoA1.CV", "Asy")])

CV.Asy[1:50, 1:5]

colnames(CV.Asy) = c("Ele", "pre.bio12", "pre.AI", 
                     "PCoA1.CV", "Asy")


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

write.csv(hori.all, "../2.DX2013_145sites/data/plant.PCoA1.1314.CV.Asy.elevation.csv")


ggplot(hori.all, aes(x=Ele, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.4)+
  geom_smooth(method = "glm", formula = y ~ poly(x, 2),
              se = T, color = "red")+
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))

ggsave("../2.DX2013_145sites/figures/plant1314.PCoA1.CV.9site.AI.pdf",
       width = 6.5, height = 3)


################# prec  ############### 
CV.Asy.md <- melt(CV.Asy, id=c("Ele", "pre.bio12", "pre.AI"))

prec.CV.me = dcast(CV.Asy.md, pre.bio12  ~ variable, mean)


prec.CV.me.sd = dcast(CV.Asy.md, pre.bio12  ~ variable, sd)


# number of each elevation
prec.hori.len <- aggregate(CV.Asy.md$value, by=list(CV.Asy.md$pre.bio12,  CV.Asy.md$variable), 
                           FUN=length)


prec.CV.me1 = melt(prec.CV.me, id = c("pre.bio12"))


# mean and se for hori
prec.hori.all <- data.frame(prec.CV.me1, se = (prec.CV.me1$value)/sqrt(prec.hori.len$x))


write.csv(prec.hori.all, "../2.DX2013_145sites/data/plant.PCoA1.1314.CV.Asy.precipitation.csv")


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

ggsave("../2.DX2013_145sites/figures/prec.plant1314.PCoA1.CV.9site.AI.pdf",
       width = 6.5, height = 3)




CV.all = subset(prec.hori.all, prec.hori.all$variable == "PCoA1.CV")
CV.all = CV.all[,-4]

library(segmented)
x1=CV.all$pre.bio12
y1=CV.all$value
Turb1<-glm(y1~x1,data=CV.all)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
summary(Turb1.seg)


pdf("../2.DX2013_145sites/figures/prec.plant1314.PCoA1.CV.9site.all.tipping.pdf",
    width = 5, height = 3)
# edge distance 
par(mar=c(4.5,4.5,4.5,4.5))
plot(x1,y1,xlab='Ele', ylab='PCoA1.CV',
     pch = 16, cex = 2, col = "#64C9CF")
plot(Turb1.seg, add=TRUE, link=FALSE, lwd=2, col="red",
     lty=1, conf.level=0.95, shade=F,)
lines(Turb1.seg, col="red", pch = 19, lwd = 2)
title(font.lab = 2, cex.lab = 2)
dev.off()


#### Asy ####
Asy.all = subset(prec.hori.all, prec.hori.all$variable == "Asy")
Asy.all = Asy.all[,-4]


x1=Asy.all$pre.bio12
y1=Asy.all$value
Turb1<-glm(y1~x1,data=Asy.all)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
summary(Turb1.seg)


pdf("../2.DX2013_145sites/figures/prec.plant1314.PCoA1.Asy.9site.all.tipping.pdf",
    width = 6.5, height = 3)
# edge distance 
par(mar=c(4.5,4.5,4.5,4.5))
plot(x1,y1,xlab='Ele', ylab='Asy',
     pch = 16, cex = 2, col = "#64C9CF")
plot(Turb1.seg, add=TRUE, link=FALSE, lwd=2, col="red",
     lty=1, conf.level=0.95, shade=F,)
lines(Turb1.seg, col="red", pch = 19, lwd = 2)
title(font.lab = 2, cex.lab = 2)
dev.off()



#########  pre.AI, CV and asy ########### 
CV.Asy[1:5, 1:5]

AI.CV.Asy2 = CV.Asy[, c("pre.AI", "PCoA1.CV", "Asy")]

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

write.csv(AI.hori.all, "../2.DX2013_145sites/data/plant.PCoA1.1314.CV.Asy.aridity.csv")

# AI.hori.all = read.csv("../2.DX2013_145sites/data/plant.PCoA1.1314.CV.Asy.aridity.csv")

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

pdf("../2.DX2013_145sites/figures/Aridity.plant1314.PCoA1.Asy.9site.pdf",
    width = 6.5, height = 3)