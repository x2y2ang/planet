
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


lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(CV.Asy.all[,"gp"])
length(gp)

resil <- as.data.frame(CV.Asy.all[,c("value")])
colnames(resil) <- c("resil")

AI <- as.data.frame(CV.Asy.all[,c("AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- resil[which(CV.Asy.all$gp == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(CV.Asy.all$gp == gp.tmp),]
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

unique(CV.Asy.all$variable)
var.list = c( "Plant.Mass.CV", "Micro.PLFA.CV", 
              "Plant.Mass.Asy",  "Micro.PLFA.Asy",
              "EF.CV", "Plant.Rich.CV",  
              "EF.Asy", "Plant.Rich.Asy",
              "Micro.Rich.CV", " Soil.PCA1.CV ", 
              "Micro.Rich.Asy", " Soil.PCA1.Asy",
              " Plant.PCoA1.CV",  "Micro.PCoA1.CV",
              "Plant.PCoA1.Asy ", "Micro.PCoA1.Asy" )

CV.Asy.all$variable=factor(CV.Asy.all$variable,levels = var.list)
CV.Asy.all=CV.Asy.all[order(CV.Asy.all$variable),]


envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = CV.Asy.all$gp %in% c(lm.t)

envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = CV.Asy.all$gp %in% c(lm.f)


CV.Asy.all$Elevation = factor(CV.Asy.all$Elevation)

colnames(CV.Asy.all)

se = as.numeric(CV.Asy.all$se)


ggplot(CV.Asy.all, aes(x=AI, y=value, color = gp))+ 
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se))+
  
  stat_smooth(aes(colour=gp),method="lm",formula=y~x,size=0.8, se = TRUE,
              data = CV.Asy.all[envi.lm.T1,])+
  # linear lm sig==F
  stat_smooth(aes(colour=gp),method="lm",formula=y~x,size=0.8, lty = 2,  se = FALSE,
              data = CV.Asy.all[envi.lm.F1,])+

  facet_wrap( ~ variable, scales="free_y", ncol=4)+
  
  scale_colour_manual(values = rep(c("#1F78B4", "#33A02C"), 16))+
  # scale_colour_manual(values = rep(c("#A6CEE3", "#B2DF8A"), 16))+
  # scale_colour_manual(values = rep(c("#1F78B4", "#6A3D9A"), 16))+
  
  xlab("Aridity") + ylab("")+
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 9))+
  theme( panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 9),
        axis.text=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9))+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5))+
  theme(legend.position = 'none')

ggsave("../2.DX2013_145sites/figures/AI.last.DX1314.all.cv.Asy.nonlegend.pdf",
       width = 6.5, height = 7)
