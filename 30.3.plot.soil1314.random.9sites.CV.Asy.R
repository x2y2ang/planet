
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# library
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

prec.hori.all = read.csv("../2.DX2013_145sites/data/soil.rich.1314.CV.Asy.precipitation.delete.5200.csv")

colnames(prec.hori.all)

prec.hori.all$Elevation = factor(prec.hori.all$Elevation)


ggplot(prec.hori.all, aes(x=pre.prec, y=value, color = Elevation,
                        fill = Elevation))+ 
  geom_point(size = 3)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se))+
  geom_errorbar(aes(xmin=value-se, xmax=value+se))+
  # scale_colour_manual(values = brewer.pal(10, 'RdYlBu'))
  scale_colour_manual(values = brewer.pal(12,"Paired"))+
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
  xlab("") + ylab("")+
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 9))+
  theme( panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 9),
        axis.text=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9))+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5))+
  xlim(300, 510)


##################################
library(vegan)
library(MASS)
lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(prec.hori.all[,"Group"])
length(gp)

resil <- as.data.frame(prec.hori.all[,c("value")])
colnames(resil) <- c("resil")

prec <- as.data.frame(prec.hori.all[,c("pre.prec")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- resil[which(prec.hori.all$Group == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  prec.t <- prec[which(prec.hori.all$Group == gp.tmp),]
  prec.t = data.frame(prec.t)
  
  
  gp.prec.tmp = data.frame(gp.t = as.vector(gp.t),
                           prec.t = as.vector(prec.t))
  
  lm.tmp = summary(lm(gp.t ~ prec.t, gp.prec.tmp))
  
  table.prec.tmp = data.frame(group = gp[i],
                              slope = lm.tmp$coefficients[2,1],
                              Ine = lm.tmp$coefficients[1,1],
                              p = lm.tmp$coefficients[2,4],
                              r2 = lm.tmp$adj.r.squared)
  if (i==1)
  {
    env.gp = table.prec.tmp
  }
  else{
    env.gp = rbind(env.gp, table.prec.tmp)
  }
}


lm.rich.all = env.gp
# p.value > 0.05
lm.rich.all$sig = ifelse(lm.rich.all$p > 0.05, "F","T")

lm.T1 = prec.hori.all$Group %in% c("A", "E")
lm.F1 = prec.hori.all$Group %in% c("B", "D", "F")

D = subset(prec.hori.all, Group == "D")
summary(lm(D$value ~ D$pre.prec))
# p-value: 0.07257

ggplot(prec.hori.all, aes(x=pre.prec, y=value,  color = Elevation))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se))+
  geom_errorbar(aes(xmin=value-se, xmax=value+se))+
  # linear lm sig==T
  stat_smooth(aes(colour=Group),method="lm",formula=y~x,size=0.8,
              se = FALSE,data = prec.hori.all[lm.T1,])+
  # linear lm sig==F
  stat_smooth(aes(colour=Group),method="lm",formula=y~x,size=0.8, lty = 2, 
              se = FALSE,data = prec.hori.all[lm.F1,])+

  xlab("Precipitation (mm)")+ylab("")+
  scale_colour_manual(values = c(brewer.pal(9,"Paired"), "black", "black",
                                 "black", "black", "black"))+
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
  ggtitle("")+
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
  xlim(300, 510)


pdf("../2.DX2013_145sites/figures/Fig.72.prec.soil1314.pca1.Asy.9site.pdf",
    width = 6.5, height = 3)
