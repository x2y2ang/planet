
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# library
library(vegan)
library(readxl)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(permute)
library(lattice)
library(splitstackshape)
library(MASS)
library(igraph)
library(impute)
library(GO.db)
library(preprocessCore)
library(AnnotationDbi)
library(psych)
library(Hmisc)


# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")

load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")

PLFA.all = read_xls("../2.DX2013_145sites/data/PLFA1314_add_ele_aridity_combine_bac.act.xls")

colnames(PLFA.all)[1] = c("")

rownames(PLFA.all) =  as.matrix(PLFA.all[, 1])

# PLFA.all = PLFA.all[, -1]

PLFA.temp = PLFA.all[, c(4, 6:12)]


lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

colnames(PLFA.temp)

# last envi ele patterns
ENVI.LAST = PLFA.temp

# choose ele
ele.c <- ENVI.LAST[,c("pre.AI")]
# cbind ID to envi
ele.id <- cbind(rownames(ele.c),ele.c)
colnames(ele.id) <- c("pre.AI")


colnames(ENVI.LAST)
envi.id <- cbind(rownames(ENVI.LAST),ENVI.LAST[,2:8])
names(envi.id)[1] <- "ID"
envi.melt <- melt(envi.id, id="ID")
envi <- unique(envi.melt[,"variable"])
length(envi)
envi.value <- as.data.frame(envi.melt[,c("value")])
rownames(envi.value) <- rownames(envi.melt$ID)
colnames(envi.value) <- c("envi.var")

i=1
for (i in 1:length(envi)){
  envi.tmp <- envi[i]
  envi.t <- envi.value[which(envi.melt$variable == envi.tmp),]
  ele.t <- ele.id[,1]

  ele.t = as.numeric(ele.t)
  
  #lm
  envi.ele.tmp = data.frame(envi.t = as.vector(envi.t),
                            ele.t = as.vector(ele.t))
  
  step.aic.tmp=stepAIC(lm(envi.t~ele.t+I(ele.t^2), envi.ele.tmp),trace = FALSE)
  tmp.var <- attr(terms(step.aic.tmp), "term.labels")
  
  qua.tmp=lm(envi.t~ele.t+I(ele.t^2), envi.ele.tmp)
  summary(qua.tmp)
  AIC.qua=AIC(lm(envi.t~ele.t+I(ele.t^2), envi.ele.tmp))
  lin.tmp=lm(envi.t~ele.t,envi.ele.tmp)
  summary(lin.tmp)
  AIC.lin=AIC(lm(envi.t~ele.t,envi.ele.tmp))
  
  if(length(tmp.var)==0){
    if(AIC.qua<AIC.lin){
      env.envi.tem=data.frame(
        envi = envi[i],
        model="quadratic",
        p.value=lm_pvalue(qua.tmp),
        r2=summary(qua.tmp)$adj.r.squared)
      
    }
    else{
      env.envi.tem=data.frame(
        envi = envi[i],
        model="linear",
        p.value=lm_pvalue(lin.tmp),
        r2=summary(lin.tmp)$adj.r.squared)
    }
  }else if(length(tmp.var)==1){
    env.envi.tem=data.frame(
      envi = envi[i],
      model="linear",
      p.value=lm_pvalue(lin.tmp),
      r2=summary(lin.tmp)$adj.r.squared)
    
  }else if(length(tmp.var)==2){
    env.envi.tem=data.frame(
      envi = envi[i],
      model="quadratic",
      p.value=lm_pvalue(qua.tmp),
      r2=summary(qua.tmp)$adj.r.squared)
  }
  if (i==1)
  {
    env.envi = env.envi.tem
  }
  else{
    env.envi = rbind(env.envi,env.envi.tem)
  }
}

lm.ele.envi <- env.envi
# p.value > 0.05
lm.ele.envi$sig = ifelse(lm.ele.envi$p.value > 0.05, "F","T")
# lm.ele.envi$sig = ifelse(lm.ele.envi$p.value > 0.01, "F1","T1")

# Melt the data by altitude into a long table
ENVI.LAST1 = cbind(PLFA.all$Elevation, ENVI.LAST)
colnames(ENVI.LAST1)[1] = c("Elevation")

envi.tmp <- melt(ENVI.LAST1, id=c("pre.AI", "Elevation"))
colnames(envi.tmp) <- c("pre.AI", "Elevation", "Variables","Value")

envi.tmp = subset(envi.tmp, Variables != "Total")

# choose linear lm P<0.05 in lm.ele.envi
envi.lm.T <- lm.ele.envi[lm.ele.envi$model=="linear" &
                           lm.ele.envi$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = envi.tmp$Variables %in% c(lm.t)


# choose quadratic  lm P<0.05 in lm.ele.envi
envi.qu.T <- lm.ele.envi[lm.ele.envi$model=="quadratic" &
                           lm.ele.envi$sig=="T",]
qu.t <- as.vector(envi.qu.T[,1])

# choose quadratic  lm P<0.05 in envi.tmp
envi.qu.T1 = envi.tmp$Variables %in% c(qu.t)


# pl.tmep = melt(PLFA.temp, id = c("Elevation", "pre.AI"))

envi.tmp$Elevation = factor(envi.tmp$Elevation)


unique(envi.tmp$Variables)
var.list = c("Fungal", "Bacterial",  "Fungal/bacterial",
             "Gram-positive",   "Gram-negative",              
             "Gram-positive/gram-negative")

envi.tmp$Variables=factor(envi.tmp$Variables,levels = var.list)
envi.tmp=envi.tmp[order(envi.tmp$Variables),]


ggplot(envi.tmp, aes(x = pre.AI,y = Value, color = Elevation))+
  geom_point(size = 2)+
  facet_wrap( ~ Variables, scales ="free_y", ncol=3)+
  stat_smooth(method="lm",formula=y~x,size=0.8, color="Black",
              se = TRUE, data = envi.tmp[envi.lm.T1,])+
  # quadratic lm sig==T
  stat_smooth(method="lm",formula=y~x+I(x^2),size=0.8, color="Black",
              se = TRUE, data = envi.tmp[envi.qu.T1,])+
  xlab("Aridity")+ylab("")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+

  theme_bw()+
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(axis.text.x  = element_text( vjust=0.5))+ # angle=45,
  theme(strip.background = element_blank())+
  theme(text = element_text())+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(),
        axis.title.x=element_text(),
        axis.title.y=element_text())

ggsave("../2.DX2013_145sites/figures/plot.PLFA.bac.fungal.last.pdf",
       width = 6.5, height = 4.3)


###################################

envi.mass.not.gram = read.csv("../2.DX2013_145sites/data/envi.1314.add.mass.27.7.2022.correct.csv",
                              row.names = 1)
cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

micro.mass = data.frame(envi.mass.not.gram$Elevation, cli1314.last$pre.AI, 
                        envi.mass.not.gram$PLFA) 
colnames(micro.mass) = c("Elevation", "pre.AI", "Micro.mass")

micro.mass.me <- reshape2::melt(micro.mass, id=c("Elevation", "pre.AI"))


lm.micro.mass = lm(micro.mass.me$value ~ micro.mass.me$pre.AI, micro.mass.me)
summary(lm.micro.mass)

qu.micro.mass = lm(micro.mass.me$value ~ micro.mass.me$pre.AI + 
                     I((micro.mass.me$pre.AI)^2), micro.mass.me)
summary(qu.micro.mass)
# Adjusted R-squared: 0.5959  p-value: < 2.2e-16

AIC(lm.micro.mass, qu.micro.mass)
# qu.micro.mass  4 1897.805

micro.mass.me$Elevation = factor(micro.mass.me$Elevation)

p4 = ggplot(micro.mass.me, aes(pre.AI, value, color = Elevation))+
  geom_point(size = 1.8)+
  stat_smooth(method="lm", formula=y~x+I(x^2), 
              color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
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
  xlab("Aridity")+ylab("Microbial PLFA")

p4

ggsave("../2.DX2013_145sites/figures/plot.PLFA.bac.fungal.last.pdf",
       width = 3.2, height = 2.3)



source("../2.DX2013_145sites/function/diversity.functions.R")

dim(comm1314.samp)

micro1314 = quick.diversity(comm1314.samp)


micro.coh = read.csv("../2.DX2013_145sites/data/micro.cohesion.last1314.csv")



all.mass.rich.coh = cbind(micro.mass, micro1314$comm.richness,
                            micro.coh$total.coh)

write.csv(all.mass.rich.coh, "micro.rich.inter.biomass.combine.csv")



colnames(all.mass.rich.coh) = c("Elevation", "pre.AI", "Micro.mass",             
                                "Micro.rich", "Micro.net")

all = read.csv("micro.rich.inter.biomass.combine.csv", row.names = 1)

cor.all<-rcorr(as.matrix(as.data.frame(lapply(all, as.numeric))), 
                 type = "pearson")
cor.all.r = cor.all$r 
cor.all.p = cor.all$P 



## aridity <0.69

all.below = subset(all, Elevation > 4700)

cor.below<-rcorr(as.matrix(as.data.frame(lapply(all.below, as.numeric))), 
                 type = "pearson")
cor.below.r = cor.below$r 
cor.below.p = cor.below$P 


## aridity >0.69

all.above = subset(all, Elevation <= 4600)

cor.above<-rcorr(as.matrix(as.data.frame(lapply(all.above, as.numeric))), 
             type = "pearson")
cor.above.r = cor.above$r 
cor.above.p = cor.above$P 



# library
library(vegan)
library(FactoMineR)
library(lavaan)
library(semPlot)
library(AICcmodavg)
library(data.table)
library(MASS)


net.high.data = all.above

net.high.SEM.md = na.omit(net.high.data)
net.high.SEM = scale(net.high.SEM.md,center=T,scale=T) 
net.high.SEM = as.data.frame(net.high.SEM)
net.high.SEM = data.table(net.high.SEM)
colnames(net.high.SEM)


lvmod.1 <- ' 
# Regressions
 Micro.rich ~ pre.AI 
 Micro.net ~ pre.AI + Micro.rich
 Micro.mass ~ pre.AI + Micro.rich + Micro.net
'
lvmod.1.fit <- sem(lvmod.1, data=net.high.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.1.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.2 <- ' 
# Regressions
 Micro.rich ~ pre.AI 
 Micro.net ~ pre.AI + Micro.rich
 Micro.mass ~ pre.AI + Micro.net # + Micro.rich
'
lvmod.2.fit <- sem(lvmod.2, data=net.high.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)


lvmod.3 <- ' 
# Regressions
 Micro.rich ~ pre.AI 
 Micro.net ~ pre.AI # + Micro.rich
 Micro.mass ~ pre.AI + Micro.net #+ Micro.rich 
'
lvmod.3.fit <- sem(lvmod.3, data=net.high.SEM, fixed.x=F,std.lv=TRUE,
                   orthogonal=TRUE )
summary(lvmod.3.fit, rsq=T, standardized=T,fit.measures = TRUE)

source("../2.DX2013_145sites/function/lavaan.modavg.R")


aictab.lavaan(list(lvmod.1.fit,lvmod.2.fit,lvmod.3.fit),
              c("model1","model2","model3")) 


semPaths(lvmod.2.fit, what = 'std', layout = 'tree', residuals = FALSE,
         edge.label.cex = 1)


summary(lvmod.2.fit, rsq=T, standardized=T,fit.measures = TRUE)

fitMeasures(lvmod.2.fit, c("cfi", "rmsea", "srmr"))


