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

unique(coh.all$Group)
var.list = c("Micro", "Veg-Micro", "Envi-Veg-Micro")
coh.all$Group=factor(coh.all$Group,levels = var.list)
coh.all=coh.all[order(coh.all$Group),]

p1 = ggplot(coh.all, aes(pre.AI, total.coh, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y ~ x + I(x^2), color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+ 
  facet_wrap( ~ Group, scales = "free", ncol = 3)+ 
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
  xlab("Aridity")+ylab("Total cohesion")+
  theme(legend.position = 'none')+
  ylim(0.29, 0.4)

p1


colnames(coh.all)

p2 = ggplot(coh.all, aes(pre.AI, pos.coh, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y ~ x + I(x^2), color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+ 
  facet_wrap( ~ Group, scales = "free", ncol = 3)+ 
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
  xlab("Aridity")+ylab("Positive cohesion")+
  theme(legend.position = 'none')+
  ylim(0.16, 0.24)

p2


colnames(coh.all)

p3 = ggplot(coh.all, aes(pre.AI, neg.pos, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y ~ x + I(x^2), color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+ 
  facet_wrap( ~ Group, scales = "free", ncol = 3)+ 
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
  xlab("Aridity")+ylab("Negtive/Positive cohesion")+
  theme(legend.position = 'none')+
  ylim(0.65, 0.9)

p3


library(cowplot)
plot_grid(p1,p2,p3, ncol = 1)
ggsave("../2.DX2013_145sites/figures/elevation.aridity.envi.plant.soil.net.cohesion.1314year.pdf",
       width = 6.5, height = 7.5)


#### for standard ####
colnames(coh.all)
p4 = ggplot(coh.all, aes(pre.AI, std.total, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y ~ x + I(x^2), color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+ 
  facet_wrap( ~ Group, scales = "free", ncol = 3)+ 
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(strip.background = element_blank())+
  theme(axis.text.x  = element_text(vjust=0.5))+
  theme(text = element_text(size = 12))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10),
        axis.text=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))+
  xlab("Aridity")+ylab("Std total cohesion")+
  theme(legend.position = 'none')+
  ylim(-0.02, 0.06)

p4


colnames(coh.all)

p5 = ggplot(coh.all, aes(pre.AI, std.pos, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y ~ x + I(x^2), color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+ 
  facet_wrap( ~ Group, scales = "free", ncol = 3)+ 
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(strip.background = element_blank())+
  theme(axis.text.x  = element_text(vjust=0.5))+
  theme(text = element_text(size = 12))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10),
        axis.text=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))+
  xlab("Aridity")+ylab("Std positive cohesion")+
  theme(legend.position = 'none')+
  ylim(-0.02, 0.06)

p5


colnames(coh.all)

p6 = ggplot(coh.all, aes(pre.AI, std.neg.pos.ratio, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y ~ x + I(x^2), color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+ 
  facet_wrap( ~ Group, scales = "free", ncol = 3)+ 
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(strip.background = element_blank())+
  theme(axis.text.x  = element_text(vjust=0.5))+
  theme(text = element_text(size = 12))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10),
        axis.text=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))+
  xlab("Aridity")+ylab("Std negtive/positive cohesion")+
  theme(legend.position = 'none')+
  ylim(-0.15, 0.1)

p6


library(cowplot)
plot_grid(p4,p5,p6, ncol = 1)
ggsave("../2.DX2013_145sites/figures/Elevation.std.mean.aridity.envi.plant.soil.net.cohesion.1314year.pdf",
       width = 6.5, height = 7.5)



lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

coh.tmp = coh.all[,c(1:9, 11)]

coh.me = melt(coh.tmp, id = c("Group", "pre.AI"))

gp <- unique(coh.me[,"variable"])
length(gp)

net <- as.data.frame(coh.me[,c("value")])
colnames(net) <- c("net")

AI <- as.data.frame(coh.me[,c("pre.AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[1]
  gp.t <- net[which(coh.me$variable == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(coh.me$variable == gp.tmp),]
  AI.t = data.frame(AI.t)
  
  
  gp.AI.tmp = data.frame(gp.t = as.vector(gp.t),
                          AI.t = as.vector(AI.t))
  
  step.aic.tmp=stepAIC(lm(gp.t~AI.t+I(AI.t^2), gp.AI.tmp),trace = FALSE)
  tmp.var <- attr(terms(step.aic.tmp), "term.labels")
  
  qua.tmp=lm(gp.t~AI.t+I(AI.t^2), gp.AI.tmp)
  summary(qua.tmp)
  AIC.qua=AIC(lm(gp.t~AI.t+I(AI.t^2), gp.AI.tmp))
  lin.tmp=lm(gp.t~AI.t,gp.AI.tmp)
  summary(lin.tmp)
  AIC.lin=AIC(lm(gp.t~AI.t,gp.AI.tmp))
  
  if(length(tmp.var)==0){
    if(AIC.qua<AIC.lin){
      env.gp.tem=data.frame(
        group = gp[i],
        model="quadratic",
        p.value=lm_pvalue(qua.tmp),
        r2=summary(qua.tmp)$adj.r.squared)
      
    }
    else{
      env.gp.tem=data.frame(
        group = gp[i],
        model="linear",
        p.value=lm_pvalue(lin.tmp),
        r2=summary(lin.tmp)$adj.r.squared)
    }
  }else if(length(tmp.var)==1){
    env.gp.tem=data.frame(
      group = gp[i],
      model="linear",
      p.value=lm_pvalue(lin.tmp),
      r2=summary(lin.tmp)$adj.r.squared)
    
  }else if(length(tmp.var)==2){
    env.gp.tem=data.frame(
      group = gp[i],
      model="quadratic",
      p.value=lm_pvalue(qua.tmp),
      r2=summary(qua.tmp)$adj.r.squared)
  }
  if (i==1)
  {
    env.envi = env.gp.tem
  }
  else{
    env.envi = rbind(env.envi,env.gp.tem)
  }
}

lm.rich.all = env.gp
# p.value > 0.05
lm.rich.all$sig = ifelse(lm.rich.all$p > 0.05, "F","T")

lm.T1 = ele.hori.all$Group %in% c("A","B", "C", "E", "G", "L", "I", "J")
lm.F1 = ele.hori.all$Group %in% c("D", "F", "H", "K")
