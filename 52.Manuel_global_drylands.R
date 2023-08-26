rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(vegan)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(openxlsx)

global = read.csv("../2.DX2013_145sites/data/Manuel_global_drylands1.csv")
colnames(global)[1] = c("Plot")

glo.tmp = global[, c(7, 10:17)]

glo.melt = melt(glo.tmp, id = c("Aridity"))

ggplot(glo.melt, aes(x=Aridity, y=value, color = variable))+ 
  geom_point(size = 2)+
  stat_smooth(method="loess",
              formula=y~x)+
  facet_wrap( ~ variable, scales = "free_y", ncol=2)+
  xlab("Aridity")+ylab("")+
  theme_bw()+
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(axis.text.x  = element_text(vjust=0.5))+ # angle=45, 
  theme(strip.background = element_blank())+
  theme(text = element_text())+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(),
        axis.title.x=element_text(),
        axis.title.y=element_text())+
  theme(axis.text.x  = element_text(vjust=0.5))+
  theme( panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5))+
  theme(legend.position = 'none')


lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


colnames(global)
gp <- unique(global[,"gp2"])
length(gp)

EF <- as.data.frame(global[,c("Mutilfunction")])
colnames(EF) <- c("EF")

AI <- as.data.frame(global[,c("Aridity")])

i=1
for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  EF.t <- EF[which(global$gp2 == gp.tmp),]
  
  EF.mean = mean(EF.t)
  EF.sd = sd(EF.t)
  
  CV.tmp= EF.sd/EF.mean
  
  # Asy strong functions for row
  EF.max = max(EF.t)
  EF.min = min(EF.t)
  
  Asy = ((EF.max - EF.mean)/(EF.mean - EF.min))
  
  ####### for mean pre.AI ########
  pre.AI.t <- AI[which(global$gp2 == gp.tmp),]
  
  pre.AI.mean = mean(pre.AI.t)
  
  
  ##########################################
  act.ele.CV.tmp = cbind(gp.tmp, pre.AI.mean, CV.tmp, Asy)
  
  colnames(act.ele.CV.tmp) = c("Group", "AI", "EF.CV", "EF.Asy")
  
  if (i==1)
  {
    act.ele.CV = act.ele.CV.tmp
  }
  else{
    act.ele.CV = rbind(act.ele.CV, act.ele.CV.tmp)
  }
}

act.ele.CV = data.frame(act.ele.CV)

act.ele.CV = act.ele.CV[order(act.ele.CV$AI), ]

# ## AL3 outlier
act.ele.CV = subset(act.ele.CV, Group != "27" & Group != "36" & Group != "40")


EF.CV.Asy = act.ele.CV[, 2:4]

AI.num = round(as.numeric(EF.CV.Asy$AI), 3)

cv = abs(as.numeric(EF.CV.Asy$EF.CV, 6))

asy = abs(as.numeric(EF.CV.Asy$EF.Asy, 6))

EF.la = cbind(AI.num, cv, asy)

EF.la = data.frame(EF.la)

colnames(EF.la) = c("AI", "EF.CV", "EF.Asy")


library(segmented)
x1=EF.la$AI
y1=EF.la$EF.CV
Turb1<-glm(y1~x1,data=EF.la)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='EF.CV')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.739  0.059



x1=EF.la$AI
y1=EF.la$EF.Asy
Turb1<-glm(y1~x1,data=EF.la)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='EF.Asy')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.855  0.043

