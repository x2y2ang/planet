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


### Plan micor Rich Soil.PCoA1 ###
Rich.Soil = CV.Asy.all[111:334, ]

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(Rich.Soil[,"Group"])
length(gp)

resil <- as.data.frame(Rich.Soil[,c("value")])
colnames(resil) <- c("resil")

AI <- as.data.frame(Rich.Soil[,c("AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- resil[which(Rich.Soil$Group == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(Rich.Soil$Group == gp.tmp),]
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

write.csv(lm.rich.all, "../2.DX2013_145sites/data/lm.aridity.plant.micro.rich.soil.PCOA1.csv")

unique(Rich.Soil$variable)
var.list = c("Plant.Rich.CV", "Micro.Rich.CV", " Soil.PCA1.CV ",
             "Plant.Rich.Asy", "Micro.Rich.Asy", " Soil.PCA1.Asy")

Rich.Soil$variable=factor(Rich.Soil$variable,levels = var.list)
Rich.Soil=Rich.Soil[order(Rich.Soil$variable),]

envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = Rich.Soil$Group %in% c(lm.t)

envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = Rich.Soil$Group %in% c(lm.f)


Rich.Soil$Elevation = factor(Rich.Soil$Elevation)


ggplot(Rich.Soil, aes(x=AI, y=value, color = Elevation))+
  geom_point(size = 2)+
  # linear lm sig==T
  stat_smooth(aes(colour=Group),method="lm",formula=y~x,size=0.8, se = TRUE,
              data = Rich.Soil[envi.lm.T1,])+
  # linear lm sig==F
  stat_smooth(aes(colour=Group),method="lm",formula=y~x,size=0.8, lty = 2,  se = FALSE,
              data = Rich.Soil[envi.lm.F1,])+
  
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


ggsave("../2.DX2013_145sites/figures/AI.last.DX1314.plant.micro.rich.soil.cv.Asy.nonlegend.pdf",
       width = 6.5, height = 5.5)

