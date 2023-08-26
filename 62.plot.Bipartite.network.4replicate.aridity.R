
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

BI.all = read.csv("../2.DX2013_145sites/data/62.Bipattite.plant.microbes.4replicates.attribute.all.last.csv")
colnames(BI.all)

BI.me = melt(BI.all, id = c("Aridity", "Elevation",
                            "Replicate", "Group"))

BI.all1 = BI.all[, c(1:4, 8:11, 14:18)]
BI.me1 = melt(BI.all1, id = c("Aridity", "Elevation", 
                              "Replicate", "Group"))

Bi.mean = dcast(BI.me1, Aridity + Elevation ~ variable, mean)

Bi.sd = dcast(BI.me1, Aridity + Elevation ~ variable, sd)

Bi.length = dcast(BI.me1, Aridity + Elevation ~ variable, length)

Bi.sd.la = Bi.sd/2

write.csv(Bi.mean, "../2.DX2013_145sites/final.figures2023/Bipattite.plant.microbes.mean.4replicates.attribute.all.last.csv")

write.csv(Bi.sd.la, "../2.DX2013_145sites/final.figures2023/Bipattite.plant.microbes.sd.last.4replicates.attribute.all.last.csv")


# BI.me1$Elevation = factor(BI.me1$Elevation)
ggplot(BI.me1, aes(x=Aridity, y=value, color=Elevation))+ 
  facet_wrap( ~ variable, scales="free_y", ncol=3)+

#  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
#                                 "#DAE2B6", "#B2DF8A", "#6CC4A1", "#43919B", "#1F78B4"))+
  geom_point(size = 2)+
  stat_smooth(method="loess")+
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
  xlab("Aridity")+ylab("")

ggsave("../2.DX2013_145sites/final.figures2023/Figure.plot.raw.nonliner.Bipartite.plant.microbial.AI.patterns.pdf",
       width = 8, height = 7.5)

### delete consistent patterns

BI.all1 = BI.all[, c(1:4, 8:11, 14:18)]

BI.me1 = melt(BI.all1, id = c("Aridity", "Elevation", 
                              "Replicate", "Group"))

write.csv(BI.me1, "../2.DX2013_145sites/data/63.Bipartite.all1.delete.some.attri.all.csv")


BI.me1.la = read.csv("../2.DX2013_145sites/data/63.Bipartite.all1.delete.some.attri.all.add.group.csv",
                     row.names = 1)

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

write.csv(lm.rich.all, "../2.DX2013_145sites/data/64.lm.qu.Bipartite.plant.microbial.AI.4replictes.csv")


unique(BI.me1.la$variable)
var.list = c("Total.Links", "Connectance","Links.per.species", "Cluster.coefficient",                     
             "Linkage.density", "Robustness.at.the.plant.level",           
             "Robustness.at.the.microbial.level", "Partner.diversity.at.the.plant.level",    
             "Partner.diversity.at.the.microbial.level") 

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
                                 "#F94A29","#1F78B4"))+ #"#E31A1C", 数量一定要对上

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
  xlab("Aridity")+ylab("")

ggsave("../2.DX2013_145sites/final.figures2023/Figure.plot.Bipartite.plant.microbial.AI.patterns.pdf",
       width = 8, height = 7.5)

### PCA 
library(FactoMineR)
BI.AI.tmp = BI.all1[, 4:12]
BI.pca = PCA(BI.AI.tmp, scale.unit=T, graph=F)
BI.pca$ind$coord
summary(BI.pca)

# Eigenvalues
# Dim.1   Dim.2   Dim.3   Dim.4   Dim.5   Dim.6   Dim.7
# Variance               6.461   1.600   0.535   0.273   0.070   0.030   0.021
# % of var.             71.790  17.774   5.941   3.030   0.782   0.329   0.230
# Cumulative % of var.  71.790  89.564  95.505  98.535  99.317  99.646  99.876

BI.pca.t <- BI.pca$ind$coord[,1:2]
colnames(BI.pca.t) = c("Bi.pca1", "Bi.pca2")

BI.LAST = cbind(BI.all1[, c("Elevation","Aridity", "Group")], 
                BI.pca.t)


BI.LAST$Elevation <- as.factor(BI.LAST$Elevation)
BI.LAST$Aridity <- as.factor(BI.LAST$Aridity)

ggplot(BI.LAST, aes(x=Bi.pca1, y=Bi.pca2, color = Aridity)) + # shape = Elevation 
  scale_shape_manual(values= c(1:10))+
  geom_point(size=3)+
  scale_colour_manual(values = c("#1F78B4","#43919B","#6CC4A1", "#B2DF8A","#DAE2B6", 
                                 "#E5BA73","#FDBF6F","#FF7F00","#FB9A99","#E31A1C"))+
  theme_bw()+ 
  guides(col = guide_legend(nrow = 5))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"))+
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
  xlim(-7,9)+ylim(-4,4)

ggsave("../2.DX2013_145sites/final.figures2023/Figure.plot.PCA.Bipartite.plant.microbial.AI.patterns.pdf",
       width = 6.5, height = 5)


ggplot(BI.LAST, aes(x=Bi.pca1, y=Bi.pca2, color = Group)) +
  scale_shape_manual(values= c(1:10))+
  geom_point(size=3)+
  theme_bw()+ 
  guides(col = guide_legend(nrow = 5))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"))+
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
  stat_ellipse(aes(fill=BI.LAST$Group))

ggsave("../2.DX2013_145sites/final.figures2023/Figure.plot.PCA.ellipse.Bipartite.plant.microbial.AI.patterns.pdf",
       width = 6.5, height = 5)

