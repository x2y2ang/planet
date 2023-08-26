
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

ele.hori.all = read.csv("../2.DX2013_145sites/data/aridity.micro.plant.soil.rich.1314.CV.Asy.csv")

colnames(ele.hori.all)

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(ele.hori.all[,"Group"])
length(gp)

resil <- as.data.frame(ele.hori.all[,c("value")])
colnames(resil) <- c("resil")

ele <- as.data.frame(ele.hori.all[,c("pre.AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- resil[which(ele.hori.all$Group == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  ele.t <- ele[which(ele.hori.all$Group == gp.tmp),]
  ele.t = data.frame(ele.t)
  
  
  gp.ele.tmp = data.frame(gp.t = as.vector(gp.t),
                          ele.t = as.vector(ele.t))
  
  lm.tmp = summary(lm(gp.t ~ ele.t, gp.ele.tmp))
  
  table.ele.tmp = data.frame(group = gp[i],
                             slope = lm.tmp$coefficients[2,1],
                             Ine = lm.tmp$coefficients[1,1],
                             p = lm.tmp$coefficients[2,4],
                             r2 = lm.tmp$adj.r.squared)
  if (i==1)
  {
    env.gp = table.ele.tmp
  }
  else{
    env.gp = rbind(env.gp,table.ele.tmp)
  }
}


lm.rich.all = env.gp
# p.value > 0.05
lm.rich.all$sig = ifelse(lm.rich.all$p > 0.05, "F","T")

lm.T1 = ele.hori.all$Group %in% c("A", "D", "E",  "G", "H", "I", "J", "K", "L")
lm.F1 = ele.hori.all$Group %in% c("F","B", "C")

ele.hori.all$Elevation = factor(ele.hori.all$Elevation)

ggplot(ele.hori.all, aes(x=pre.AI, y=value, color = Elevation))+
  geom_point(size = 3)+
  # linear lm sig==T
  stat_smooth(aes(colour=Group),method="lm",formula=y~x,size=0.8,
              se = FALSE,data = ele.hori.all[lm.T1,])+
  # linear lm sig==F
  stat_smooth(aes(colour=Group),method="lm",formula=y~x,size=0.8, lty = 2,
              se = FALSE,data = ele.hori.all[lm.F1,])+

  xlab("")+ylab("")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired"), 
                                 "black", "black", "black",
                                 "black", "black", "black",
                                 "black", "black", "black",
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
        axis.title.y=element_text(size=12)) 
  xlim(4300, 5200)

ggsave("../2.DX2013_145sites/figures/Aridity.micro.plant.soil.1314.rich.Asy.9site.pdf",
    width = 6.5, height = 5)
