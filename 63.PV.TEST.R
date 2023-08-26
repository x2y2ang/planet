
set.seed(123)

a = sample(50:100,3,replace=T)

# 75 56 91

a.mean = mean(a)

a.sd = sd(a)

# sqrt(40.835)


b = sample(50:100,3,replace=T)

# 63 74 75

b.mean = mean(b)

b.sd = sd(b)


c = sample(50:100,3,replace=T)

# 78 84 57 

(78+84+57)/3

c.mean = mean(c)

c.sd = sd(c)


a.CV = sd(a)/mean(a)

b.CV = sd(b)/mean(b)

c.CV = sd(c)/mean(c)


((20/76)+(51/91))/3

((20/77)+(18/76))/3

((33/84)+(21/78))/3


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


load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")

colnames(envi)

micro.vec.hori = cbind(envi[, c("Aridity.index","Act.ele", "Elevation", "Distance_between_origin",
                                "Distance_origin", "Direction_origin",       
                                "X" ,"Y")], Micro.alp)

vec.tmp = micro.vec.hori[, c("Aridity.index", "Elevation", "Micro.rich")] 

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

# pv.fun(a)
# pv.fun(b)
# pv.fun(c)

micro.rich.tmp = data.frame(vec.tmp[, "Micro.rich"])
Act.ele = data.frame(vec.tmp[, c("Aridity.index")])

ele = unique(vec.tmp[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- micro.rich.tmp[which(vec.tmp$Elevation == ele.tmp),]
  
  # rich.com <- t(combn(rich.t,9))
  # rich.com <- data.frame(rich.com)
  # rownames(rich.com)
  
  Rich.PV <- pv.fun(rich.t)
 
  
  # for mean ele
  act.ele.t <- Act.ele[which(vec.tmp$Elevation == ele.tmp),]
  
  # act.ele.com <- t(combn(act.ele.t,9))
  # act.ele.com <- data.frame(act.ele.com)
  # rownames(act.ele.com)
  
  # for mean
  # act.ele.mean = apply(act.ele.t, 1, mean)
  
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


micro.rich.PV = data.frame(act.ele.PV)
plot(micro.rich.PV$Aridity, micro.rich.PV$Rich.PV)

write.csv(micro.rich.PV, "../2.DX2013_145sites/data/65.1.PV.Test.micro.rich.csv")



# 
# 
# ### party starts ###
# ### Add k=1 to both disparity and PV ###
# 
# # generate dataframes for PV, D, AR1 - all variables
# ## cams
# cams.pv <- matrix(ncol=12, nrow=nrow(cams.res))
# cams.d <- matrix(ncol=12, nrow=nrow(cams.res))
# cams.ar1 <- matrix(ncol=12, nrow=nrow(cams.res))
# cams.mean <- matrix(ncol=12, nrow=nrow(cams.res))
# 
# ## Loops to calculate pv, d, mean per pixel
# aux <- numeric()
# for (i in 1:nrow(cams.res)){ 
#   for (y in 1:12){
#     for (x in 1:((length(colnames(cams.res))-2)/12)){
#       a <- cams.res[,(3-12+((x*12):(11+x*12)))] # we choose the 12 months of the X year
#       # 1981 >>> 3:14
#       # 1982 >>> 15:26
#       # 1983 >>> 27:38 # it works!
#       aux [x] <- a [i,y] # we select the month of interest
#     }
#     cams.pv [i,y] <- pv.fun(aux, k=1)
#     cams.d [i,y] <- d.fun(aux, k=1)
#     # cams.ar1 [i,y] <- ar1.fun(aux)
#     cams.mean [i,y] <- mean(aux, na.rm=T)
#   }
# }
