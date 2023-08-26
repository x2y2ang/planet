
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# library
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(MASS)

CV.Asy.all = read.csv("../2.DX2013_145sites/data/AI.micro.plant.soil.1314.CV.Asy.mean.se.all.last.csv")

colnames(CV.Asy.all)[1] = c("Elevation")

unique(CV.Asy.all$variable)
Asy.all = subset(CV.Asy.all, variable == " Plant.PCoA1.CV" | 
                   variable =="Plant.PCoA1.Asy " |
                   variable == "Micro.PCoA1.CV" |
                   variable =="Micro.PCoA1.Asy") 

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(Asy.all[,"gp"])
length(gp)

resil <- as.data.frame(Asy.all[,c("value")])
colnames(resil) <- c("resil")

AI <- as.data.frame(Asy.all[,c("AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- resil[which(Asy.all$gp == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(Asy.all$gp == gp.tmp),]
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


lm.rich.all = env.gp
# p.value > 0.05
lm.rich.all$sig = ifelse(lm.rich.all$p > 0.05, "F","T")

write.csv(lm.rich.all, "../2.DX2013_145sites/data/lm.p.value.plant.micro.PCoA1.CV.Asy.csv")

unique(Asy.all$variable)
var.list = c(" Plant.PCoA1.CV",  "Micro.PCoA1.CV",
             "Plant.PCoA1.Asy ", "Micro.PCoA1.Asy") 

Asy.all$variable=factor(Asy.all$variable,levels = var.list)
Asy.all=Asy.all[order(Asy.all$variable),]


envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = Asy.all$gp %in% c(lm.t)

envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = Asy.all$gp %in% c(lm.f)


Asy.all$Elevation = factor(Asy.all$Elevation)

colnames(Asy.all)

se = as.numeric(Asy.all$se)

ggplot(Asy.all, aes(x=AI, y=value, color = gp))+ 
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se))+
  
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, se = TRUE,
              data = Asy.all[envi.lm.T1,])+
  # linear lm sig==F
  
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, lty = 2,  se = FALSE,
              data = Asy.all[envi.lm.F1,])+
  
  facet_wrap( ~ variable, scales="free_y", ncol = 2)+
  
  scale_colour_manual(values = rep(c("#1F78B4", "#33A02C"), 6))+
  
  theme_bw()+ 
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
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
  xlab("Aridity")+ylab("")+
  theme(legend.position = 'none')

ggsave("../2.DX2013_145sites/figures/plot.plant.micro.pcoa1.cv.Asy.pdf",
       width = 4.5, height = 5)


########## Plant.PCoA1 ################
library(segmented)
library(ggplot2)

Plant.PCoA1 = subset(Asy.all, variable == " Plant.PCoA1.CV")

windowsFonts(RMN=windowsFont("Times New Roman"))
x1=Plant.PCoA1$AI
y1=Plant.PCoA1$value
Turb1<-glm(y1~x1,data=Plant.PCoA1)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='Plant.PCoA1.CV')

# plot(x1,y1,ylim=c(0,1))

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.329  0.354


#### Micro.PCoA1.CV ###
unique(Asy.all$variable)
Micro.PCoA1.CV = subset(Asy.all, variable == "Micro.PCoA1.CV")

windowsFonts(RMN=windowsFont("Times New Roman"))
x1=Micro.PCoA1.CV$AI
y1=Micro.PCoA1.CV$value
Turb1<-glm(y1~x1,data=Micro.PCoA1.CV)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='Micro.PCoA1.CV')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.864  0.031



#### Plant.PCoA1.Asy ###
unique(Asy.all$variable)
Plant.PCoA1.Asy = subset(Asy.all, variable == "Plant.PCoA1.Asy")

windowsFonts(RMN=windowsFont("Times New Roman"))
x1=Plant.PCoA1.Asy$AI
y1=Plant.PCoA1.Asy$value
Turb1<-glm(y1~x1,data=Plant.PCoA1.Asy)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='Plant.PCoA1.Asy')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.781  0.117



#### Plant.PCoA1.Asy  ###
unique(Asy.all$variable)
Plant.PCoA1.Asy = subset(Asy.all, variable == "Plant.PCoA1.Asy ")

windowsFonts(RMN=windowsFont("Times New Roman"))
x1=Plant.PCoA1.Asy$AI
y1=Plant.PCoA1.Asy$value
Turb1<-glm(y1~x1,data=Plant.PCoA1.Asy)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='Plant.PCoA1.Asy')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)




#### Micro.PCoA1.Asy ###
unique(Asy.all$variable)
Micro.PCoA1.Asy = subset(Asy.all, variable == "Micro.PCoA1.Asy")

windowsFonts(RMN=windowsFont("Times New Roman"))
x1=Micro.PCoA1.Asy$AI
y1=Micro.PCoA1.Asy$value
Turb1<-glm(y1~x1,data=Micro.PCoA1.Asy)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='Micro.PCoA1.Asy')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.342  0.219


#### Soil.PCA1.Asy  ###
unique(Asy.all$variable)
Soil.PCA1.Asy  = subset(Asy.all, variable == " Soil.PCA1.Asy")

windowsFonts(RMN=windowsFont("Times New Roman"))
x1=Soil.PCA1.Asy $AI
y1=Soil.PCA1.Asy $value
Turb1<-glm(y1~x1,data=Soil.PCA1.Asy )
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
plot(x1,y1, xlab='AI', ylab='Soil.PCA1.Asy')

plot(Turb1.seg,add=TRUE,link=FALSE,lwd=2,col=1,  
     lty=1, conf.level=0.95,shade=F,)
lines(Turb1.seg,col=2,pch=19,lwd=2)
summary(Turb1.seg)

# Estimated Break-Point(s):
#   Est. St.Err
# psi1.x1 0.782  0.122