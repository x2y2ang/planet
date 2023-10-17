
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
getwd()

# library
library(vegan)
library(SpatialEpi)
library(SciViews)
library(segmented)
library(ggplot2)
library(MASS)
library(reshape2)
library(dplyr)
library(plyr)
library(RColorBrewer)

BI.all = read.csv("../data/GCB/Bipattite.plant.microbes.8replicates.attribute.all.last.csv")
colnames(BI.all)

BI.me = melt(BI.all, id = c("Aridity", "Elevation",
                            "Replicate", "Group"))

BI.all1 = BI.all[, c(1:4, 8:11, 14:18)]
BI.me1 = melt(BI.all1, id = c("Aridity", "Elevation", 
                              "Replicate", "Group"))

Bi.mean = dcast(BI.me1, Aridity + Elevation ~ variable, mean)

Bi.sd = dcast(BI.me1, Aridity + Elevation ~ variable, sd)

Bi.length = dcast(BI.me1, Aridity + Elevation ~ variable, length)

Bi.se.la = Bi.sd/sqrt(Bi.length)

write.csv(Bi.mean, "../data/GCB/Bipattite.plant.microbes.mean.8replicates.attribute.all.last.csv")

write.csv(Bi.se.la, "../data/GCB/Bipattite.plant.microbes.se.last.8replicates.attribute.all.last.csv")

write.csv(BI.me1, "../data/GCB/Bipartite.8replicates.all1.delete.some.attri.all.csv")


BI.me1.la = read.csv("../data/GCB/Bipartite.8replicates.all1.delete.some.attri.all.add.group.csv",
                     row.names = 1)

unique(BI.me1.la$variable)

BI.me1.la = subset(BI.me1.la, variable == "Connectance" | variable == "Cluster.coefficient" |
                   variable == "Robustness.at.the.microbial.level")


lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(BI.me1.la[,"gp"])
length(gp)

net <- as.data.frame(BI.me1.la[,c("value")])
colnames(net) <- c("net")


AI <- as.data.frame(BI.me1.la[,c("Aridity")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- net[which(BI.me1.la$gp == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(BI.me1.la$gp == gp.tmp),]
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

write.csv(lm.rich.all, "../data/GCB/lm.qu.Bipartite.8replicates.plant.microbial.Aridity.csv")


unique(BI.me1.la$variable)
var.list = c("Connectance", "Cluster.coefficient",              
             "Robustness.at.the.microbial.level") 

BI.me1.la$variable=factor(BI.me1.la$variable,levels = var.list)
BI.me1.la=BI.me1.la[order(BI.me1.la$variable),]


envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = BI.me1.la$gp %in% c(lm.t)

envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = BI.me1.la$gp %in% c(lm.f)


BI.me1.la$Elevation = factor(BI.me1.la$Elevation)
BI.me1.la$gp = factor(BI.me1.la$gp)
colnames(BI.me1.la)

# net.lm.T1 = BI.me1$Group %in% c("Below")
# net.lm.F1 = BI.me1$Group %in% c("Above")

ggplot(BI.me1.la, aes(x=Aridity, y=value, color=Elevation))+ 
  facet_wrap( ~ variable, scales="free_y", ncol=3)+
  geom_point(size = 2, alpha = 0.8)+
  # geom_smooth(method="loess", size=1, se = F)+
  stat_smooth(aes(colour=Group),method="lm",
              formula=y~x,size=1, se = TRUE,
              data = BI.me1.la[envi.lm.T1,])+
  # linear lm sig==F
  
  stat_smooth(aes(colour=Group),method="lm",
              formula=y~x,size=1, lty = 2,  se = FALSE,
              data = BI.me1.la[envi.lm.F1,])+
  
  scale_colour_manual(values = c("#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1", "#43919B", "#1F78B4", 
                                 "#F94A29","#1F78B4",
                                 "#F94A29","#1F78B4", 
                                 "#F94A29","#1F78B4",
                                 "#F94A29","#1F78B4",
                                 "#F94A29","#1F78B4",
                                 "#F94A29","#1F78B4",
                                 "#F94A29","#1F78B4",
                                 "#F94A29","#1F78B4",
                                 "#F94A29","#1F78B4"))+ 
  
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
  xlab("Aridity")+ylab("")+
  theme(legend.position = 'none')
  

ggsave("../GCB_Figure/Fig.S11.plot.Bipartite.8replicates.plant.microbial.Aridity.patterns.pdf",
       width = 6.5, height = 2.8)


