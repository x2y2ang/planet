
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

PV.all = read.csv("../2.DX2013_145sites/data/PV.random9sample.all.components.repli8times.deal.with.PLFA.plant.4300m.true.form.csv")
colnames(PV.all)[1] = c("Elevation")

PV.all.tmp = read.csv("../2.DX2013_145sites/data/PV.random9sample.all.components.repli8times.deal.with.PLFA.plant.4300m.csv")
colnames(PV.all.tmp)[1] = c("variable")

PV.all.tmp1 = PV.all.tmp[,c(1:11)]

PV.all = melt(PV.all.tmp1, id = c("variable", "Elevation", "Aridity"))

write.csv(PV.all, "../2.DX2013_145sites/data/melt.PV.random9sample.all.components.repli8times.deal.with.PLFA.plant.4300m.csv")
  
lm_pvalue <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  
gp <- unique(PV.all[,"gp"])
length(gp)
  
resil <- as.data.frame(PV.all[, c("value")])
colnames(resil) <- c("resil")
  
AI <- as.data.frame(PV.all[,c("AI")])
  
  
for (i in 1:length(gp)){
    gp.tmp <- gp[i]
    gp.t <- resil[which(PV.all$gp == gp.tmp),]
    gp.t = data.frame(gp.t)
    
    AI.t <- AI[which(PV.all$gp == gp.tmp),]
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
  
write.csv(lm.rich.all, "../2.DX2013_145sites/data/lm.p.value.AGB.PLFA.EF.rich.soil.PV.csv")
  
unique(PV.all$Group)
var.list = c( "Plant.rich", "Micro.rich", "EMF",
              "AGB", "PLFA", "Soil.pca1")
PV.all$Group=factor(PV.all$Group,levels = var.list)
PV.all=PV.all[order(PV.all$Group),]
  
  
envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = PV.all$gp %in% c(lm.t)
  
envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = PV.all$gp %in% c(lm.f)
  
  
PV.all$Elevation = factor(PV.all$Elevation)
PV.all$gp = factor(PV.all$gp)

colnames(PV.all)
  
se = as.numeric(PV.all$se)
  
ggplot(PV.all, aes(x=AI, y=value, color = Elevation))+ 
    geom_point(size = 2)+
    geom_errorbar(aes(ymin=value-se, ymax=value+se))+
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, se = TRUE,
              data = PV.all[envi.lm.T1,])+
  # linear lm sig==F
  
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, lty = 2,  se = FALSE,
              data = PV.all[envi.lm.F1,])+
  
  facet_wrap( ~ Group, scales="free_y", ncol=3)+
  
  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1",  "#43919B", "#1F78B4", 
                                 "#F94A29", "#1F78B4", 
                                 "#F94A29", "#1F78B4",
                                 "#F94A29", "#1F78B4", 
                                 "#F94A29", "#1F78B4", 
                                 "#F94A29", "#1F78B4", 
                                 "#F94A29", "#1F78B4"))+
  
  # scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+ 
  xlab("Aridity") + ylab("PV")+
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))+
  theme( panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5))+
  theme(legend.position = 'none')
  
  
ggsave("../2.DX2013_145sites/figures/plot.AGB.PLFA.EF.rich.soil.PV.pdf",
         width = 6.5, height = 5)



