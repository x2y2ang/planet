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

micro.vec.hori = cbind(envi[, c("Act.ele", "Elevation", "Distance_between_origin",
                                "Distance_origin", "Direction_origin",       
                                "X" ,"Y")], Micro.alp)

vec.tmp = micro.vec.hori[, c("Act.ele", "Elevation", "Micro.rich")] 

# for sd
rowVars <- function(x, na.rm=F) {
  # Vectorised version of variance filter
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}


################     loop for 10 ele    ######################

micro.rich.tmp = data.frame(vec.tmp[, "Micro.rich"])
Act.ele = data.frame(vec.tmp[, c("Act.ele")])

ele = unique(vec.tmp[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- micro.rich.tmp[which(vec.tmp$Elevation == ele.tmp),]
  
  rich.com <- t(combn(rich.t,9))
  rich.com <- data.frame(rich.com)
  rownames(rich.com)
  
  rich.mean = apply(rich.com, 1, mean)
  rich.sd = apply(rich.com, 1, sd)
  
  CV.tmp= rich.sd/rich.mean
  
  # Asy strong functions for row
  rich.max = apply(rich.com, 1, max)
  rich.min = apply(rich.com, 1, min)
  
  Asy = ((rich.max - rich.mean)/(rich.mean - rich.min))
  
  
  # for mean ele
  act.ele.t <- Act.ele[which(vec.tmp$Elevation == ele.tmp),]
  
  act.ele.com <- t(combn(act.ele.t,9))
  act.ele.com <- data.frame(act.ele.com)
  rownames(act.ele.com)
  
  # for mean
  act.ele.mean = apply(act.ele.com, 1, mean)
  
  act.ele.CV.tmp = cbind(act.ele.mean, CV.tmp, Asy)
  
  colnames(act.ele.CV.tmp) = c("Ele", "Rich.CV", "Asy")
  
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

CV.Asy = cbind(ele.r, act.ele.CV)

CV.Asy2 = CV.Asy[, -2]

CV.Asy3 = melt(CV.Asy2, id = ("ele.r"))

# CV.Asy3$ele.r = factor(CV.Asy3$ele.r)

ggplot(CV.Asy3, aes(ele.r, value, fill = variable))+
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
  scale_fill_manual(values = c("#FF1818", "#11468F"))+ 
  facet_wrap( ~ variable, scales ="free_y", ncol = 2)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))+
  xlab("")+ylab("")

# ggsave("../2.DX2013_145sites/figures/Fig.57.micro.rich.CV.9site.AI.pdf",
#        width = 6.5, height = 3) 


CV.me = dcast(CV.Asy3, ele.r ~ variable, mean)

CV.me.sd = dcast(CV.Asy3, ele.r ~ variable, sd)

# number of each elevation
hori.len <- aggregate(CV.Asy3$value, by=list(CV.Asy3$ele.r, CV.Asy3$variable), 
                      FUN=length)

CV.me1 = melt(CV.me, id = c("ele.r"))
# mean and se for hori
hori.all <- data.frame(CV.me1, se = (CV.me1$value)/sqrt(hori.len$x))


# hori.all$ele.r = factor(hori.all$ele.r)

ggplot(hori.all, aes(x=ele.r, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.4)+
  # geom_smooth(method = "glm", formula = y ~ poly(x, 7), 
  #             se = T, color = "red")+ 
  facet_wrap( ~ variable, scales = "free", ncol = 1)+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))


ggsave("../2.DX2013_145sites/figures/Fig.58.micro.rich.CV.9site.all.pdf",
       width = 3.5, height = 5)


CV.all = subset(hori.all, hori.all$variable == "Rich.CV")
CV.all = CV.all[,-4]
  
library(segmented)
x1=CV.all$ele.r
y1=CV.all$value
Turb1<-glm(y1~x1,data=CV.all)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
summary(Turb1.seg)


pdf("../2.DX2013_145sites/figures/Fig.59.micro.rich.CV.9site.all.tipping.pdf",
    width = 4, height = 4.5)
# edge distance 
par(mar=c(4.5,4.5,4.5,4.5))
plot(x1,y1,xlab='Ele', ylab='Rich.CV',
     pch = 16, cex = 2, col = "#64C9CF")
plot(Turb1.seg, add=TRUE, link=FALSE, lwd=2, col="red",
     lty=1, conf.level=0.95, shade=F,)
lines(Turb1.seg, col="red", pch = 19, lwd = 2)
title(font.lab = 2, cex.lab = 2)
dev.off()


################################
Asy.all = subset(hori.all, hori.all$variable == "Asy")
Asy.all = Asy.all[,-4]

x2=Asy.all$ele.r
y2=Asy.all$value
Turb2<-glm(y2~x2,data=Asy.all)
Turb2.seg<-segmented(Turb2,seg.Z=~x1)


pdf("../2.DX2013_145sites/figures/Fig.60.micro.Asymmetry.9sites.all.tipping.pdf",
    width = 4, height = 4.5)
par(mar=c(4.5,4.5,4.5,4.5))

plot(x2,y2,xlab='Ele', ylab='Asymmetry',
     pch = 16, cex = 2, col = "#64C9CF")
plot(Turb2.seg, add=TRUE, link=FALSE, lwd=2, col="red",
     lty=1, conf.level=0.95, shade=F,)
lines(Turb2.seg, col="red", pch = 19, lwd = 2)
title(font.lab = 2, cex.lab = 2)
dev.off()

summary(Turb2.seg)

