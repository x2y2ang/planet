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

micro.vec.hori = cbind(envi[, c("Aridity.index","Act.ele", "Elevation", "Distance_between_origin",
                                "Distance_origin", "Direction_origin",       
                                "X" ,"Y")], Micro.alp)

vec.tmp = micro.vec.hori[, c("Aridity.index", "Elevation", "Micro.rich")] 

# for sd
rowVars <- function(x, na.rm=F) {
  # Vectorised version of variance filter
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}


################     loop for 10 ele    ######################
(3096+3187+3106+3155+3940+3006+3195+2979+3228)/9

micro.rich.tmp = data.frame(vec.tmp[, "Micro.rich"])
Act.ele = data.frame(vec.tmp[, c("Aridity.index")])

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

ele.r = round(act.ele.CV$Ele, 2)

CV.Asy = cbind(ele.r, act.ele.CV)

CV.Asy2 = CV.Asy[, -2]

CV.Asy3 = melt(CV.Asy2, id = ("ele.r"))

CV.Asy3$ele.r = factor(CV.Asy3$ele.r)

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
  xlab("Aridity")+ylab("")
ggsave("../2.DX2013_145sites/figures/Fig.57.micro.rich.CV.9site.AI.pdf",
       width = 6.5, height = 3) 


CV.me = dcast(CV.Asy3, ele.r ~ variable, mean)

CV.me.sd = dcast(CV.Asy3, ele.r ~ variable, sd)




CV.me1 = melt(CV.me, id = c("ele.r"))


ggplot(CV.me1, aes(x=ele.r, y=value))+ 
  geom_point(size = 3, color = "gray")+
  #stat_smooth(method="loess", color = "red")+
  facet_wrap( ~ variable, scales = "free")+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))



ggsave("../2.DX2013_145sites/figures/Fig.57.micro.rich.CV.3site.10eles.pdf",
       width = 3.5, height = 3) 

