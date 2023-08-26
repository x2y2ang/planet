
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")


library(vegan)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
# display.brewer.all() 
# brewer.pal(12,"Paired")

veg188.cli.tmp = read.csv("../2.DX2013_145sites/data/DX1314.Plant.rich.AI.for.PV.csv",
                          row.names = 1)

# veg188.cli.tmp = subset(veg188.cli.tmp, X < 10000)

Plant.micro.soil.EMF = read.csv("../2.DX2013_145sites/data/DX1314.Plant.micro.soil.EMF.AI.for.PV.csv",
                                row.names = 1)
# Plant.micro.soil.EMF = subset(Plant.micro.soil.EMF, X < 10000)

# load functions
d.fun <- function (var, k=0){
  var<-var[!is.na(var)]
  if(length(var)<3){return(NA)}
  if(sd(var)==0){return(0)}
  if (min(var)<0){var <- var + abs(min(var)) + 0.01*(max(var)-min(var))}
  if (length(var)<2){return(NA)}else{
    var <- var + k
    aux <- numeric(length(var)-1)
    for (i in 1:(length(var)-1)){
      aux [i] <- abs(log(var[i+1]/var[i]))
    }
  }
  return(mean(aux, na.rm=T))}
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
se.fun <- function(var){return(sd(var)/sqrt(length(var)))}
ar1.fun <- function (var){ 
  var<-var[!is.na(var)]
  if(length(var)<3){return(NA)}
  if(sd(var)==0){return(0)} else {return(as.numeric(acf(var, plot = F)[1]$acf))}} # wrong usage


plant.rich.tmp = data.frame(veg188.cli.tmp[, "Plant.rich"])
Act.ele = data.frame(veg188.cli.tmp[, c("pre.AI")])

ele = unique(veg188.cli.tmp[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- plant.rich.tmp[which(veg188.cli.tmp$Elevation == ele.tmp),]
  
  
  N <- 9
  s <- sample(x=rich.t, size=N, replace=TRUE)
  table(s)
  s
  
  Rich.PV <- pv.fun(s)
  
  
  # for mean ele
  act.ele.t <- Act.ele[which(veg188.cli.tmp$Elevation == ele.tmp),]
  
  act.ele.mean = mean(act.ele.t)
  
  act.ele.PV.tmp = cbind(act.ele.mean, Rich.PV)
  
  colnames(act.ele.PV.tmp) = c("Aridity", "Rich.PV")
  
  if (i==1)
  {
    act.ele.PV = act.ele.PV.tmp
  }
  else{
    act.ele.PV = rbind(act.ele.PV, act.ele.PV.tmp)
  }
}


plant.rich.PV = data.frame(act.ele.PV)

plot(plant.rich.PV$Aridity, plant.rich.PV$Rich.PV)

write.csv(plant.rich.PV, "../2.DX2013_145sites/data/69.third.random9sites.2.PV.Test.plant.rich.csv")



############################ micro.rich ################################

Micro.rich.tmp = data.frame(Plant.micro.soil.EMF[, "Micro.rich"])
Act.ele = data.frame(Plant.micro.soil.EMF[, c("pre.AI")])

ele = unique(Plant.micro.soil.EMF[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- Micro.rich.tmp[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  N <- 9
  s <- sample(x=rich.t, size=N, replace=TRUE)
  table(s)
  s
  
  Rich.PV <- pv.fun(s)
  
  
  # for mean ele
  act.ele.t <- Act.ele[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  act.ele.mean = mean(act.ele.t)
  
  act.ele.PV.tmp = cbind(act.ele.mean, Rich.PV)
  
  colnames(act.ele.PV.tmp) = c("Aridity", "Rich.PV")
  
  if (i==1)
  {
    act.ele.PV = act.ele.PV.tmp
  }
  else{
    act.ele.PV = rbind(act.ele.PV, act.ele.PV.tmp)
  }
}


Micro.rich.PV = data.frame(act.ele.PV)

plot(Micro.rich.PV$Aridity, Micro.rich.PV$Rich.PV)

write.csv(Micro.rich.PV, "../2.DX2013_145sites/data/69.third.random9sites.3.PV.Test.Micro.rich.csv")




############################ EMF ################################

EMF.tmp = data.frame(Plant.micro.soil.EMF[, "EMF"])
Act.ele = data.frame(Plant.micro.soil.EMF[, c("pre.AI")])

ele = unique(Plant.micro.soil.EMF[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- EMF.tmp[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  N <- 9
  s <- sample(x=rich.t, size=N, replace=TRUE)
  table(s)
  s
  
  Rich.PV <- pv.fun(s)
  
  # for mean ele
  act.ele.t <- Act.ele[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  act.ele.mean = mean(act.ele.t)
  
  act.ele.PV.tmp = cbind(act.ele.mean, Rich.PV)
  
  colnames(act.ele.PV.tmp) = c("Aridity", "Rich.PV")
  
  if (i==1)
  {
    act.ele.PV = act.ele.PV.tmp
  }
  else{
    act.ele.PV = rbind(act.ele.PV, act.ele.PV.tmp)
  }
}


EMF.PV = data.frame(act.ele.PV)

plot(EMF.PV$Aridity, EMF.PV$Rich.PV)

write.csv(EMF.PV, "../2.DX2013_145sites/data/69.third.random9sites.4.PV.Test.EMF.csv")




############################ AGB ################################

AGB.tmp = data.frame(Plant.micro.soil.EMF[, "AGB"])
Act.ele = data.frame(Plant.micro.soil.EMF[, c("pre.AI")])

ele = unique(Plant.micro.soil.EMF[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- AGB.tmp[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  N <- 9
  s <- sample(x=rich.t, size=N, replace=TRUE)
  table(s)
  s
  
  Rich.PV <- pv.fun(s)
  
  
  # for mean ele
  act.ele.t <- Act.ele[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  act.ele.mean = mean(act.ele.t)
  
  act.ele.PV.tmp = cbind(act.ele.mean, Rich.PV)
  
  colnames(act.ele.PV.tmp) = c("Aridity", "Rich.PV")
  
  if (i==1)
  {
    act.ele.PV = act.ele.PV.tmp
  }
  else{
    act.ele.PV = rbind(act.ele.PV, act.ele.PV.tmp)
  }
}


AGB.PV = data.frame(act.ele.PV)

plot(AGB.PV$Aridity, AGB.PV$Rich.PV)

write.csv(AGB.PV, "../2.DX2013_145sites/data/69.third.random9sites.5.PV.Test.AGB.csv")



############################ PLFA ################################

PLFA.tmp = data.frame(Plant.micro.soil.EMF[, "PLFA"])
Act.ele = data.frame(Plant.micro.soil.EMF[, c("pre.AI")])

ele = unique(Plant.micro.soil.EMF[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- PLFA.tmp[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  N <- 9
  s <- sample(x=rich.t, size=N, replace=TRUE)
  table(s)
  s
  
  Rich.PV <- pv.fun(s)
  
  
  # for mean ele
  act.ele.t <- Act.ele[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  act.ele.mean = mean(act.ele.t)
  
  act.ele.PV.tmp = cbind(act.ele.mean, Rich.PV)
  
  colnames(act.ele.PV.tmp) = c("Aridity", "Rich.PV")
  
  if (i==1)
  {
    act.ele.PV = act.ele.PV.tmp
  }
  else{
    act.ele.PV = rbind(act.ele.PV, act.ele.PV.tmp)
  }
}


PLFA.PV = data.frame(act.ele.PV)

plot(PLFA.PV$Aridity, PLFA.PV$Rich.PV)

write.csv(PLFA.PV, "../2.DX2013_145sites/data/69.third.random9sites.6.PV.Test.PLFA.csv")



############################ Soil.pca1 ################################

Soil.pca1.tmp = data.frame(Plant.micro.soil.EMF[, "Soil.pca1"])
Act.ele = data.frame(Plant.micro.soil.EMF[, c("pre.AI")])

ele = unique(Plant.micro.soil.EMF[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- Soil.pca1.tmp[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  N <- 9
  s <- sample(x=rich.t, size=N, replace=TRUE)
  table(s)
  s
  
  Rich.PV <- pv.fun(s)
  
  # for mean ele
  act.ele.t <- Act.ele[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  act.ele.mean = mean(act.ele.t)
  
  act.ele.PV.tmp = cbind(act.ele.mean, Rich.PV)
  
  colnames(act.ele.PV.tmp) = c("Aridity", "Rich.PV")
  
  if (i==1)
  {
    act.ele.PV = act.ele.PV.tmp
  }
  else{
    act.ele.PV = rbind(act.ele.PV, act.ele.PV.tmp)
  }
}


Soil.pca1.PV = data.frame(act.ele.PV)

plot(Soil.pca1.PV$Aridity, Soil.pca1.PV$Rich.PV)

write.csv(Soil.pca1.PV, "../2.DX2013_145sites/data/69.third.random9sites.7.PV.Test.Soil.pca1.csv")




