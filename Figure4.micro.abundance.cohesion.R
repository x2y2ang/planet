rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(vegan)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

net.all = read.csv("../2.DX2013_145sites/data/micro.cohesion.last1314.group.csv")

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(net.all[,"gp"])
length(gp)

net <- as.data.frame(net.all[,c("value")])
colnames(net) <- c("net")

AI <- as.data.frame(net.all[,c("pre.AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- net[which(net.all$gp == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(net.all$gp == gp.tmp),]
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

write.csv(lm.rich.all, "../2.DX2013_145sites/data/lm.p.value.micro1314.coh.csv")

unique(net.all$variable)
var.list = c("Total.coh", "Neg.coh", 
             "Pos.coh", "Neg.Pos" ) 

net.all$variable=factor(net.all$variable,levels = var.list)
net.all=net.all[order(net.all$variable),]


envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = net.all$gp %in% c(lm.t)

envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = net.all$gp %in% c(lm.f)


net.all$Elevation = factor(net.all$Elevation)


ggplot(net.all, aes(x=pre.AI, y=value, color = Elevation))+ 
  geom_point(size = 2, alpha = 0.8)+
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, se = TRUE,
              data = net.all[envi.lm.T1,])+
  # linear lm sig==F
  
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, lty = 2,  se = FALSE,
              data = net.all[envi.lm.F1,])+
  
  facet_wrap( ~ variable, scales="free_y", ncol = 2)+
  
  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1",  "#43919B", "#1F78B4", 
                                 "#F94A29", "#1F78B4",
                                 "#F94A29", "#1F78B4",
                                 "#F94A29", "#1F78B4",
                                 "#F94A29", "#1F78B4"))+
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


ggsave("../2.DX2013_145sites/final.figures2023/Figure4.micro1314.cohesion.AI.patterns.pdf",
       width = 4.5, height = 5)
