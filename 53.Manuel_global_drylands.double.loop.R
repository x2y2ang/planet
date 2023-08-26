rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(vegan)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(segmented)

global = read.csv("../2.DX2013_145sites/data/Manuel_global_drylands2.csv")
colnames(global)[1] = c("Plot")

glo.tmp = global[, c(7, 10:17)]


pc <- prcomp(global[, c("SOC", "TN", "TP", "Clay", "AP")],
             center = TRUE,
             scale. = TRUE)
attributes(pc)

print(pc)

pc$x

summary(pc)

# Importance of components:
#   PC1    PC2    PC3    PC4     PC5
# Standard deviation     1.5769 0.9851 0.9376 0.7607 0.29193
# Proportion of Variance 0.4973 0.1941 0.1758 0.1157 0.01704
# Cumulative Proportion  0.4973 0.6914 0.8672 0.9830 1.00000

nut.pca1.dd <- data.frame(pc$x)

colnames(global)
alp.ai.gp = cbind(global[, c("Aridity", "Group", "gp2", "Plant.cover", "Plant.rich", "Mutilfunction")], 
                  nut.pca1.dd$PC1)


colnames(alp.ai.gp)[7] = c("Soil.pca1")

global.pca1 = cbind(global, nut.pca1.dd$PC1)

write.csv(global.pca1, "../2.DX2013_145sites/data/Manuel_global_drylands2_add_soil.pca1.csv")

glo.melt = melt(alp.ai.gp, id = c("Aridity", "Group", "gp2"))

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


colnames(alp.ai.gp)
gp <- unique(alp.ai.gp[,"gp2"])
length(gp)

var <- as.data.frame(alp.ai.gp[,c("Plant.cover", "Plant.rich", "Mutilfunction",
                                  "Soil.pca1" )])

class(var)

attri <- unique(colnames(var))
length(attri)


AI <- as.data.frame(alp.ai.gp[,c("Aridity")])

i=1
for (i in 1:length(gp)){
  gp.tmp <- gp[i]

  j=1
  for(j in 1:length(attri)){
    att.tmp = attri[j]
    var1 <- var[which(alp.ai.gp$gp2 == gp.tmp),]
    var.t <- var1[, j]
  
    var.mean = mean(var.t)
    var.sd = sd(var.t)
  
    CV.tmp= var.sd/var.mean
  
    # Asy strong functions for row
    var.max = max(var.t)
    var.min = min(var.t)
  
    Asy = ((var.max - var.mean)/(var.mean - var.min))
  
    ####### for mean pre.AI ########
    pre.AI.t <- AI[which(alp.ai.gp$gp2 == gp.tmp),]
  
    pre.AI.mean = mean(pre.AI.t)
  
  
    ##########################################
    act.ele.CV.tmp = cbind(gp.tmp, att.tmp, pre.AI.mean, CV.tmp, Asy)
    
    colnames(act.ele.CV.tmp) = c("Group", "Attri", "AI", "var.CV", "var.Asy")
    
    if(j==1){
     act.ele.CV = act.ele.CV.tmp
   }
   else{
     act.ele.CV = rbind(act.ele.CV, act.ele.CV.tmp)
       } # end for scale
  }
  
  
  if (i==1)
  {
    act.CV = act.ele.CV
  }
  else{
    act.CV = rbind(act.CV, act.ele.CV)
  }
}

act.CV = data.frame(act.CV)

act.CV1 = act.CV[order(act.CV$Attri), ]

## verified it's correct ####

EF.CV.Asy = act.CV1[, 3:5]

AI.num = round(as.numeric(EF.CV.Asy$AI), 3)

cv = abs(as.numeric(EF.CV.Asy$var.CV, 6))

asy = abs(as.numeric(EF.CV.Asy$var.Asy, 6))

EF.la = cbind(act.CV1[, 1:2], AI.num, cv, asy)

EF.la = data.frame(EF.la)

colnames(EF.la)

EF.melt = melt(EF.la, id = c("Group",  "Attri",  "AI.num"))

colnames(EF.melt)
ggplot(EF.melt, aes(x=AI.num, y=value))+ 
  geom_point(size = 3)+
  stat_smooth(method="loess")+
  facet_wrap(Attri ~ variable, scales = "free_y", ncol=2)+
  # scale_colour_manual(values = brewer.pal(12,"Paired"))+
  xlab("") + ylab("")+
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 9))+
  theme( panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 9),
        axis.text=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9))+
  theme(axis.text.x = element_text( hjust = 0.5, vjust = 0.5))


ggsave("../2.DX2013_145sites/figures/Maunel.global.drylands.CV.Asy.all.pdf",
       height = 10,width = 6.5)


write.csv(EF.melt, "../2.DX2013_145sites/data/Manuel_global_drylands2.CV.Asy.csv")


### EF for tipping ###
EF.tmp = subset(EF.la, Attri == "Mutilfunction")


x1=EF.tmp$AI
y1=EF.tmp$cv
Turb1<-glm(y1~x1,data=EF.tmp)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='EF.CV')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.739  0.059



x1=EF.tmp$AI
y1=EF.tmp$asy
Turb1<-glm(y1~x1,data=EF.tmp)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='EF.Asy')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.855  0.043


### Plant.cover for tipping ###
plant.cov.tmp = subset(EF.la, Attri == "Plant.cover")


x1=plant.cov.tmp$AI
y1=plant.cov.tmp$cv
Turb1<-glm(y1~x1,data=plant.cov.tmp)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='EF.CV')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.406  0.177



x1=plant.cov.tmp$AI
y1=plant.cov.tmp$asy
Turb1<-glm(y1~x1,data=plant.cov.tmp)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='EF.Asy')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.739  0.086


### Plant.rich for tipping ###
Soil.pca1.tmp = subset(EF.la, Attri == "Plant.rich")


x1=Soil.pca1.tmp$AI
y1=Soil.pca1.tmp$cv
Turb1<-glm(y1~x1,data=Soil.pca1.tmp)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='EF.CV')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.701  0.046



x1=Soil.pca1.tmp$AI
y1=Soil.pca1.tmp$asy
Turb1<-glm(y1~x1,data=Soil.pca1.tmp)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='EF.Asy')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.701  0.153



### Plant.rich for tipping ###
Soil.pca1.tmp = subset(EF.la, Attri == "Soil.pca1")


x1=Soil.pca1.tmp$AI
y1=Soil.pca1.tmp$cv
Turb1<-glm(y1~x1,data=Soil.pca1.tmp)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='EF.CV')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.739   0.07



x1=Soil.pca1.tmp$AI
y1=Soil.pca1.tmp$asy
Turb1<-glm(y1~x1,data=Soil.pca1.tmp)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='EF.Asy')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.427  0.226


################## plot ###############

cv.asy.all = read.csv("../2.DX2013_145sites/data/Manuel_global_drylands2.CV.Asy.tipping.point.csv.csv",
                      row.names = 1)

gp <- unique(cv.asy.all[,"gp1"])
length(gp)

resil <- as.data.frame(cv.asy.all[,c("value")])
colnames(resil) <- c("resil")

AI <- as.data.frame(cv.asy.all[,c("AI.num")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- resil[which(cv.asy.all$gp1 == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(cv.asy.all$gp1 == gp.tmp),]
  AI.t = data.frame(AI.t)
  
  
  gp.AI.tmp = data.frame(gp.t = as.vector(gp.t),
                         AI.t = as.vector(AI.t))
  
  lm.tmp = summary(lm(gp.t ~ AI.t, gp.AI.tmp))
  
  table.AI.tmp = data.frame(gp = gp[i],
                            slope = lm.tmp$coefficients[2,1],
                            Ine = lm.tmp$coefficients[1,1],
                            p = lm.tmp$coefficients[2,4],
                            r2 = lm.tmp$adj.r.squared)
  if (i==1)
  {
    env.gp = table.AI.tmp
  }
  else{
    env.gp = rbind(env.gp,table.AI.tmp)
  }
}


lm.cv.asy.all = env.gp
# p.value > 0.05
lm.cv.asy.all$sig = ifelse(lm.cv.asy.all$p > 0.05, "F","T")


envi.lm.T <- lm.cv.asy.all[lm.cv.asy.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = cv.asy.all$gp1 %in% c(lm.t)

envi.lm.F <- lm.cv.asy.all[lm.cv.asy.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = cv.asy.all$gp1 %in% c(lm.f)

class(cv.asy.all)
colnames(cv.asy.all)
cv.asy.all$gp1 = factor(cv.asy.all$gp1)

ggplot(cv.asy.all, aes(x=AI.num, y=value, color = gp1))+ 
  geom_point(size = 2)+
  stat_smooth(aes(colour=gp1),method="lm",
              formula=y~x,size=1, se = TRUE,
              data = cv.asy.all[envi.lm.T1,])+
  # linear lm sig==F
  stat_smooth(aes(colour=gp1),method="lm",
              formula=y~x,size=1, lty = 2,  se = FALSE,
              data = cv.asy.all[envi.lm.F1,])+
  
  facet_wrap(Attri ~ variable, scales="free_y", ncol=2)+
  xlab("Aridity")+ylab("")+
  scale_colour_manual(values = rep(c("#1F78B4", "#33A02C"), 8))+
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


ggsave("../2.DX2013_145sites/figures/Maunel.global.drylands.CV.Asy.all.tipping.point.pdf",
       width = 6, height = 10)
