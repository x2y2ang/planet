rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(vegan)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

coh.all = read.csv("../2.DX2013_145sites/data/envi.plant.micro.1314.cohesion.last1.csv")

colnames(coh.all)[1] = c("Group")

coh.all$Elevation = factor(coh.all$Elevation)

colnames(coh.all)

micro.net = subset(coh.all, Group == "Micro")

micro.net1 = micro.net[,c(2:5, 10:11)]

coh.all = melt(micro.net1, id = c("Elevation", "pre.AI"))

write.csv(micro.net1, "../2.DX2013_145sites/data/micro.cohesion.last1314.csv")

unique(coh.all$variable)
var.list = c("total.coh",  "abs.neg.coh", "pos.coh", "neg.pos")
coh.all$variable=factor(coh.all$variable,levels = var.list)
coh.all=coh.all[order(coh.all$variable),]

ggplot(coh.all, aes(pre.AI, value, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y ~ x + I(x^2), color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
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

ggsave("../2.DX2013_145sites/figures/elevation.aridity.micro.net.cohesion.1314year.pdf",
       width = 6.5, height = 7.5)

# use SPSS add se error 
net.se = read.csv("../2.DX2013_145sites/data/micro.cohesion.last1314.add.se.csv")

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(net.se[,"gp"])
length(gp)

net <- as.data.frame(net.se[,c("value")])
colnames(net) <- c("net")

AI <- as.data.frame(net.se[,c("AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- net[which(net.se$gp == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(net.se$gp == gp.tmp),]
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

# write.csv(lm.rich.all, "../2.DX2013_145sites/data/lm.p.value.plant.micro.PCoA1.CV.Asy.csv")

unique(net.se$variable)
var.list = c("Total.cohesion", "Negtitve.cohesion", 
             "Positive.cohesion", "Neg.pos.ratio" ) 

net.se$variable=factor(net.se$variable,levels = var.list)
net.se=net.se[order(net.se$variable),]


envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = net.se$gp %in% c(lm.t)

envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = net.se$gp %in% c(lm.f)


net.se$Elevation = factor(net.se$Ele)

colnames(net.se)

se = as.numeric(net.se$se)



p1 = ggplot(net.se, aes(x=AI, y=value, color = gp))+ 
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width = 0.00)+
  
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, se = TRUE,
              data = net.se[envi.lm.T1,])+
  # linear lm sig==F
  
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, lty = 2,  se = FALSE,
              data = net.se[envi.lm.F1,])+
  
  facet_wrap( ~ variable, scales="free_y", ncol = 2)+
  
  scale_colour_manual(values = rep(c("#33A02C", "#1F78B4"), 6))+
  
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

  
p1

ggsave("../2.DX2013_145sites/figures/AI.micro.net.cohesion.1314year.last.pdf",
       width = 4.5, height = 5)


################# each site ##########
colnames(micro.net1)
micro.coh.me = melt(micro.net1, id = c("Elevation", "pre.AI"))
write.csv(micro.coh.me, "../2.DX2013_145sites/data/micro1314.coh.last.csv")

coh.gp = read.csv("../2.DX2013_145sites/data/micro1314.coh.last.gp.csv",
                  row.names = 1)


lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(coh.gp[,"gp"])
length(gp)

net <- as.data.frame(coh.gp[,c("value")])
colnames(net) <- c("net")

AI <- as.data.frame(coh.gp[,c("pre.AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- net[which(coh.gp$gp == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(coh.gp$gp == gp.tmp),]
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

write.csv(lm.rich.all, "../2.DX2013_145sites/data/lm.p.value.micro.coh.all.csv")

unique(coh.gp$variable)
var.list = c("total.coh", "abs.neg.coh", "pos.coh", "neg.pos") 

coh.gp$variable=factor(coh.gp$variable,levels = var.list)
coh.gp=coh.gp[order(coh.gp$variable),]


envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = coh.gp$gp %in% c(lm.t)

envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = coh.gp$gp %in% c(lm.f)


coh.gp$Elevation = factor(coh.gp$Ele)

colnames(coh.gp)

ggplot(coh.gp, aes(x=pre.AI, y=value, color = gp))+ 
  geom_point(size = 2)+
  # geom_errorbar(aes(ymin=value-se, ymax=value+se), width = 0.00)+
  
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, se = TRUE,
              data = coh.gp[envi.lm.T1,])+
  # linear lm sig==F
  
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, lty = 2,  se = FALSE,
              data = coh.gp[envi.lm.F1,])+
  
  facet_wrap( ~ variable, scales="free_y", ncol = 2)+
  
  scale_colour_manual(values = rep(c("#33A02C", "#1F78B4"), 6))+
  
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


ggsave("../2.DX2013_145sites/figures/AI.micro.net.cohesion.1314year.last.pdf",
       width = 4.5, height = 5)

# sessionInfo()
