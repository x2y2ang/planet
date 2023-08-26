
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(vegan)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
# display.brewer.all() 
# brewer.pal(12,"Paired")


load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")

colnames(envi)

source('../1.DX2013_145sites/functions/diversity.functions.R')

soil = envi[, c(21:53)]

# except for SWC
soil.la = soil[,-3]
colnames(soil.la)

# soil richness
soil.alp <- quick.diversity(soil.la)
colnames(soil.alp) = c("Soil.rich", "Soil.shan", "Soil.even")

vec.tmp = cbind(envi[, c("Act.ele", "Elevation")], 
                soil.alp[, c("Soil.shan")])

colnames(vec.tmp)[3] = c("Soil.shan")

# for sd
rowVars <- function(x, na.rm=F) {
  # Vectorised version of variance filter
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}


Soil.shan.tmp = data.frame(vec.tmp[, "Soil.shan"])
Act.ele = data.frame(vec.tmp[, c("Act.ele")])

ele = unique(vec.tmp[,"Elevation"])

i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- Soil.shan.tmp[which(vec.tmp$Elevation == ele.tmp),]
  
  rich.com <- t(combn(rich.t,3))
  rich.com <- data.frame(rich.com)
  rownames(rich.com)
  
  # a = data.frame(1:455)
  rich.com1 = data.frame(rownames(rich.com), rich.com)
  colnames(rich.com1) = c("order", "rich1", "rich2", "rich3")
  
  # for mean
  mean.tmp = data.frame(rowMeans(rich.com1[, c("rich1", "rich2", "rich3")]))
  
  
  sd.tmp = data.frame(sqrt(rowVars(rich.com1[, c("rich1", "rich2", "rich3")])))
  
  CV.tmp= sd.tmp/mean.tmp
  
  # Asy strong functions for row
  rich.t = data.frame(rich.t)
  rich.max = apply(rich.com, 1, max)
  rich.min = apply(rich.com, 1, min)
  rich.mean = apply(rich.com, 1, mean)
  
  Asy = ((rich.max - rich.mean)/(rich.mean - rich.min))
  
  
  # for mean ele
  act.ele.t <- Act.ele[which(vec.tmp$Elevation == ele.tmp),]
  
  act.ele.com <- t(combn(act.ele.t,3))
  act.ele.com <- data.frame(act.ele.com)
  rownames(act.ele.com)
  
  act.ele.com1 = data.frame(rownames(act.ele.com), act.ele.com)
  colnames(act.ele.com1) = c("order", "ele1", "ele2", "ele3")
  
  # for mean
  act.mean.tmp = data.frame(rowMeans(act.ele.com1[, c( "ele1", "ele2", "ele3")]))
  
  act.ele.CV.tmp = cbind(act.mean.tmp, CV.tmp, Asy)
  
  colnames(act.ele.CV.tmp) = c("Ele", "Rich.CV", "Asy")
  
  if (i==1)
  {
    act.ele.CV = act.ele.CV.tmp
  }
  else{
    act.ele.CV = rbind(act.ele.CV, act.ele.CV.tmp)
  }
}


p1 = ggplot(act.ele.CV, aes(x=Ele, y=Rich.CV))+ 
  geom_point(size = 0.5, color = "gray")+
  stat_smooth(method="loess", color = "red")+
  xlab("Elevation") + ylab("soil shannon CV") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))
p1

ggsave("../2.DX2013_145sites/figures/Fig.55.soil.shan.CV.3site.10eles.pdf",
       width = 3.5, height = 3) 


act.ele.asy = act.ele.CV[, c(1,3)]
act.ele.asy1 = na.omit(act.ele.asy)

p2 = ggplot(act.ele.asy1, aes(x=Ele, y=Asy))+ 
  geom_point(size = 0.5, color = "gray")+
  stat_smooth(method="loess", color = "red")+
  xlab("Elevation") + ylab("Soil Asy") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))
p2


ggsave("../2.DX2013_145sites/figures/Fig.56.soil.rich.Asy.3site.10eles.pdf",
       width = 3.5, height = 3) 



library(segmented)
x1=act.ele.CV$Ele
y1=act.ele.CV$Rich.CV
Turb1<-glm(y1~x1,data=act.ele.CV)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)


plot(x1,y1,xlab='Ele', ylab='Rich.CV',
     pch = 16, cex = 2, col = "#64C9CF")
plot(Turb1.seg, add=TRUE, link=FALSE, lwd=2, col="#64C9CF",
     lty=1, conf.level=0.95, shade=F,)
lines(Turb1.seg, col="#64C9CF", pch = 19, lwd = 2)
title(font.lab = 2, cex.lab = 2)


summary(Turb1.seg)




