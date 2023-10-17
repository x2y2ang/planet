
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

envi.mass = read.csv("../data/GCB/envi.1314.add.mass.27.7.2022.correct.csv",
                     row.names = 1)

### PLFA only three replicates ####

PLFA.Origin = read.csv("../data/GCB/DX1314.PLFA.orgin.csv",
                       row.names = 1)

PLFA.LA = envi.mass[rownames(PLFA.Origin), ]

pv.fun1 <- function(var1){
  if(var1[1]==var1[2]){return(0)}else{
    return((1-(min(var1[1],var1[2])/max(var1[1],var1[2]))))}
}
pv.fun <- function(var, k=0){
  var<-var[!is.na(var)]
  if(length(var)<3){return(NA)}
  if (min(var)<0){var <- var + abs(min(var)) + 0.01*(max(var)-min(var))}
  var <- var + k
  if(length(which(is.na(var)))>=(length(var)-1)){return(NA)}else{
    var<-as.numeric(na.omit(var))
    if (any(var<0)){var <- var + abs(min(var)) + 1}
    return(mean(as.numeric(combn(var, 2, FUN = pv.fun1, simplify = F))))}}


PLFA.tmp = data.frame(PLFA.LA[, "PLFA"])
ele = unique(PLFA.LA[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  data <- PLFA.tmp[which(PLFA.LA$Elevation == ele.tmp),]
  # Calculate mean and standard deviation
  Dmean <- mean(data)
  Dsd <- sd(data)
  
  # Asy functions
  Dmax <- max(data)
  Dmin <- min(data)
  
  cv_values <- (Dsd / Dmean)
  asy_values <- ((Dmax - Dmean)/(Dmean - Dmin))
  pv_values<- pv.fun(data)
  ele.cv.asy.pv <- data.frame(cbind(ele.tmp, cv_values, asy_values, pv_values))
  
  if (i==1)
  {
    PLFA.cv.asy.pv = ele.cv.asy.pv
  }
  else{
    PLFA.cv.asy.pv = rbind(PLFA.cv.asy.pv, ele.cv.asy.pv)
  }
}

write.csv(PLFA.cv.asy.pv, "../data/GCB/PLFA.cv.asy.pv.calulate.once.csv")

CV.Asy.all = read.csv("../data/GCB/AI.micro.plant.soil.1314.CV.Asy.mean.se.all.last.csv")

colnames(CV.Asy.all)[1] = c("Elevation")

CV.all = subset(CV.Asy.all, variable == "Plant.Mass.CV" | variable =="Micro.PLFA.CV" |
                  variable == "EF.CV" | variable =="Plant.Rich.CV" | 
                  variable =="Micro.Rich.CV" | variable =="Soil.PCA1.CV")

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(CV.all[,"gp"])
length(gp)

resil <- as.data.frame(CV.all[,c("value")])
colnames(resil) <- c("resil")

AI <- as.data.frame(CV.all[,c("AI")])


for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- resil[which(CV.all$gp == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(CV.all$gp == gp.tmp),]
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

write.csv(lm.rich.all, "../data/GCB/lm.p.value.AGB.PLFA.EF.rich.soil.CV.csv")

# order variable
unique(CV.all$variable)
var.list = c( "Plant.Rich.CV", "Micro.Rich.CV",
              "EF.CV", "Plant.Mass.CV", "Micro.PLFA.CV",
              "Soil.PCA1.CV")

CV.all$variable=factor(CV.all$variable,levels = var.list)
CV.all=CV.all[order(CV.all$variable),]


envi.lm.T <- lm.rich.all[lm.rich.all$sig=="T",]
lm.t <- as.vector(envi.lm.T[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.T1 = CV.all$gp %in% c(lm.t)

envi.lm.F <- lm.rich.all[lm.rich.all$sig=="F",]
lm.f <- as.vector(envi.lm.F[,1])
# choose linear lm P<0.05 in envi.tmp
envi.lm.F1 = CV.all$gp %in% c(lm.f)


CV.all$Elevation = factor(CV.all$Elevation)
CV.all$gp = factor(CV.all$gp)

colnames(CV.all)

se = as.numeric(CV.all$se)

ggplot(CV.all, aes(x=AI, y=value, color = Elevation))+ 
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se))+
  
  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, se = TRUE,
              data = CV.all[envi.lm.T1,])+ #
  # linear lm sig==F

  stat_smooth(aes(colour=gp),method="lm",
              formula=y~x,size=1, lty = 2,  se = FALSE,
              data = CV.all[envi.lm.F1,])+ # 

  facet_wrap( ~ variable, scales="free_y", ncol=3)+
  
  scale_colour_manual(values = c("#E31A1C", "#FB9A99","#FF7F00","#FDBF6F","#E5BA73",
                                 "#DAE2B6", "#B2DF8A", "#6CC4A1",  "#43919B", "#1F78B4", 
                                 "#1F78B4", "#F94A29", 
                                 "#1F78B4", "#F94A29", 
                                 "#1F78B4", "#F94A29", 
                                 "#1F78B4", "#F94A29", 
                                 "#1F78B4", "#F94A29", 
                                 "#1F78B4", "#F94A29"))+
  theme_bw()+ 
  xlab("Aridity") + ylab("")+
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


ggsave("../GCB_Figure/Fig.S3.AGB.PLFA.EF.rich.soil.CV.pdf", width = 6.5, height = 5)


