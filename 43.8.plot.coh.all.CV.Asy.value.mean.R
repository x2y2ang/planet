rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(vegan)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


coh.last = read.csv("../2.DX2013_145sites/data/Aridity-cohesion-CV-Asy-all-mean.csv")
colnames(coh.last)[1] = c("pre.AI")

a = dcast(coh.last, variable ~ Elevation, mean)
b = dcast(coh.last, variable ~ Elevation, sd)
b[is.na(b)] <- 0
c = dcast(coh.last, variable ~ Elevation, length)
c1 = as.matrix(as.data.frame(lapply(c[, 2:11], as.numeric)))

std_error <-  b[, 2:11] / sqrt(length(c1))

mean.coh = t(a)

se.coh = t(std_error)

coh.ris = rbind(mean.coh, se.coh)

write.csv(coh.ris, "../2.DX2013_145sites/data/Cohesion.all.CV.Asy.mean.se.csv")


se.mean.coh = read.csv("../2.DX2013_145sites/data/Cohesion.all.CV.Asy.mean.se.order.csv")

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(se.mean.coh[,"Group"])
length(gp)

net <- as.data.frame(se.mean.coh[,c("mean")])
colnames(net) <- c("net")

AI <- as.data.frame(se.mean.coh[,c("AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- net[which(se.mean.coh$Group == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(se.mean.coh$Group == gp.tmp),]
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

write.csv(lm.rich.all, "../2.DX2013_145sites/data/lm.p.value.AI.se.coh.CV.Asy.csv")

unique(se.mean.coh$variable)
var.list = c("plant.net.CV", "micro.net.CV", "pms.net.CV",
             "plant.net.Asy", "micro.net.Asy", "pms.net.Asy") 

se.mean.coh$variable=factor(se.mean.coh$variable,levels = var.list)
se.mean.coh=se.mean.coh[order(se.mean.coh$variable),]


envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = se.mean.coh$Group %in% c(lm.t)

envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = se.mean.coh$Group %in% c(lm.f)


se.mean.coh$Elevation = factor(se.mean.coh$Elevation)

colnames(se.mean.coh)

# se = as.numeric(se.mean.coh$se)


ggplot(se.mean.coh, aes(x=AI, y=mean, color = Group))+ 
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width = 0.00)+
  
  stat_smooth(aes(colour=Group),method="lm",
              formula=y~x,size=1, se = TRUE,
              data = se.mean.coh[envi.lm.T1,])+
  # linear lm sig==F
  
  stat_smooth(aes(colour=Group),method="lm",
              formula=y~x,size=1, lty = 2,  se = FALSE,
              data = se.mean.coh[envi.lm.F1,])+
  
  facet_wrap( ~ variable, scales="free_y", ncol = 3)+
  
  scale_colour_manual(values = rep(c("#1F78B4", "#33A02C"), 6))+
  
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


ggsave("../2.DX2013_145sites/figures/AI.mean.se.cohesion.CV.Asy.last.pdf",
       width = 6.5, height = 5)

