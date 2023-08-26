
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

CV.Asy.all = read.csv("../2.DX2013_145sites/data/Aridity.micro.plant.soil.rich.1314.add.plant.micro.pca1.CV.Asy.csv")

colnames(CV.Asy.all)[1] = c("Elevation")


### AGB PLFA EF ###
Mass.EF = CV.Asy.all[1:110, ]

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(Mass.EF[,"Group"])
length(gp)

resil <- as.data.frame(Mass.EF[,c("value")])
colnames(resil) <- c("resil")

AI <- as.data.frame(Mass.EF[,c("AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- resil[which(Mass.EF$Group == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(Mass.EF$Group == gp.tmp),]
  AI.t = data.frame(AI.t)
  
  
  gp.AI.tmp = data.frame(gp.t = as.vector(gp.t),
                         AI.t = as.vector(AI.t))
  
  lm.tmp = summary(lm(gp.t ~ AI.t, gp.AI.tmp))
  
  table.AI.tmp = data.frame(group = gp[i],
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

unique(Mass.EF$variable)
var.list = c("Plant.Mass.CV", "Micro.PLFA.CV", "EF.CV", 
             "Plant.Mass.Asy", "Micro.PLFA.Asy", "EF.Asy")

Mass.EF$variable=factor(Mass.EF$variable,levels = var.list)
Mass.EF=Mass.EF[order(Mass.EF$variable),]


envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = Mass.EF$Group %in% c(lm.t)

envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = Mass.EF$Group %in% c(lm.f)


Mass.EF$Elevation = factor(Mass.EF$Elevation)


ggplot(Mass.EF, aes(x=AI, y=value, color = Elevation))+
  geom_point(size = 2)+
  # linear lm sig==T
  stat_smooth(aes(colour=Group),method="lm",formula=y~x,size=0.8, se = TRUE,
              data = Mass.EF[envi.lm.T1,])+
  # linear lm sig==F
  stat_smooth(aes(colour=Group),method="lm",formula=y~x,size=0.8, lty = 2,  se = FALSE,
             data = Mass.EF[envi.lm.F1,])+
  
  xlab("")+ylab("")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired"), 
                                 "black", "black", "black",
                                 "black", "black", "black",
                                 "black", "black", "black",
                                 "black", "black", "black"))+
  facet_wrap( ~ variable, scales = "free", ncol = 3)+ 
  xlab("Aridity")+
  theme_bw()+  
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
  theme(legend.position = 'none')


ggsave("../2.DX2013_145sites/figures/AI.last.DX1314.AGB.PLFA.EF.cv.Asy.nonlegend.pdf",
       width = 6.5, height = 5.5)

