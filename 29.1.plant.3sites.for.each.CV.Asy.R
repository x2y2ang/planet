
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


vec.tmp = envi[, c("Act.ele", "Elevation", "Veg.rich")] 

# for sd
rowVars <- function(x, na.rm=F) {
  # Vectorised version of variance filter
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}

################     loop for 10 ele    ######################

envi110 = envi[rownames(Veg.div),]

vec.tmp = envi110[, c("Act.ele", "Elevation", "Veg.rich")]

Veg.rich.tmp = data.frame(vec.tmp[, "Veg.rich"])
Act.ele = data.frame(vec.tmp[, c("Act.ele")])

ele = unique(vec.tmp[,"Elevation"])

i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- Veg.rich.tmp[which(vec.tmp$Elevation == ele.tmp),]
  
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
  xlab("Elevation") + ylab("Rich.CV") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))
p1

ggsave("../2.DX2013_145sites/figures/Fig.51.veg.rich.CV.3site.10eles.pdf",
       width = 3.5, height = 3) 


act.ele.asy = act.ele.CV[, c(1,3)]
act.ele.asy1 = na.omit(act.ele.asy)

p2 = ggplot(act.ele.asy1, aes(x=Ele, y=Asy))+ 
  geom_point(size = 0.5, color = "gray")+
  stat_smooth(method="loess", color = "red")+
  xlab("Elevation") + ylab("Asy") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))
p2


ggsave("../2.DX2013_145sites/figures/Fig.52.veg.rich.Asy.3site.10eles.pdf",
       width = 3.5, height = 3) 



library(segmented)
x1=act.ele.CV$Ele
y1=act.ele.CV$Asy
Turb1<-glm(y1~x1,data=act.ele.CV)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)

plot(x1,y1,xlab='Ele', ylab='Asy',
     pch = 16, cex = 2, col = "#64C9CF")
plot(Turb1.seg, add=TRUE, link=FALSE, lwd=2, col="#64C9CF",
     lty=1, conf.level=0.95, shade=F,)
lines(Turb1.seg, col="#64C9CF", pch = 19, lwd = 2)
title(font.lab = 2, cex.lab = 2)


summary(Turb1.seg)



#  delect 4300 4400

act.plant = subset(act.ele.CV, Ele > 4400)

p3 = ggplot(act.plant, aes(x=Ele, y=Rich.CV))+ 
  geom_point(size = 0.5, color = "gray")+
  stat_smooth(method="loess", color = "red")+
  xlab("Elevation") + ylab("Plant richness CV") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))
p3

ggsave("../2.DX2013_145sites/figures/Fig.53.veg.delete.43.44.rich.CV.3site.10eles.pdf",
       width = 3.5, height = 3) 


act.plant.asy = act.plant[, c(1,3)]
act.plant.asy1 = na.omit(act.plant.asy)

p4 = ggplot(act.plant.asy1, aes(x=Ele, y=Asy))+ 
  geom_point(size = 0.5, color = "gray")+
  stat_smooth(method="loess", color = "red")+
  xlab("Elevation") + ylab("Asy") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))
p4


ggsave("../2.DX2013_145sites/figures/Fig.52.veg.rich.Asy.de.43.44.3site.10eles.pdf",
       width = 3.5, height = 3) 



library(segmented)
x1=act.plant.asy1$Ele
y1=act.plant.asy1$Asy
Turb1<-glm(y1~x1,data=act.plant.asy1)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)

plot(x1,y1,xlab='Ele', ylab='Asy',
     pch = 16, cex = 2, col = "#64C9CF")
plot(Turb1.seg, add=TRUE, link=FALSE, lwd=2, col="#64C9CF",
     lty=1, conf.level=0.95, shade=F,)
lines(Turb1.seg, col="#64C9CF", pch = 19, lwd = 2)
title(font.lab = 2, cex.lab = 2)


summary(Turb1.seg)



