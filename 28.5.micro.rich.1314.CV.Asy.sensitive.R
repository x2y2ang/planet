
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# library
library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(splitstackshape)


# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")


load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")


source("../2.DX2013_145sites/function/diversity.functions.R")
dim(comm1314.samp)
micro1314 = quick.diversity(comm1314.samp)

cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

cli1314.last = cbind(cli1314.last, envi1314[, c("SWC")])

colnames(cli1314.last)[11] = c("SWC")


################    loop for 10 eles    ######################
vec.tmp = cbind(cli1314.last, micro1314, envi1314)

micro.rich.tmp = data.frame(vec.tmp[, "comm.richness"])
Act.ele = data.frame(vec.tmp[, c("Act.ele")])
pre.bio12 = data.frame(vec.tmp[, c("pre.bio12")])
pre.AI = data.frame(vec.tmp[, c("pre.AI")])
SWC = data.frame(vec.tmp[, c("SWC")])


ele = unique(vec.tmp[,"Elevation"])



i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- micro.rich.tmp[which(vec.tmp$Elevation == ele.tmp),]
  
  rich.com <- t(combn(rich.t,9))
  rich.com <- data.frame(rich.com)

  
  ######## for mean ele #######
  act.ele.t <- Act.ele[which(vec.tmp$Elevation == ele.tmp),]
  
  act.ele.com <- t(combn(act.ele.t,9))
  act.ele.com <- data.frame(act.ele.com)

  ######## for prec ele #######
  pre.bio12.t <- pre.bio12[which(vec.tmp$Elevation == ele.tmp),]
  
  pre.bio12.com <- t(combn(pre.bio12.t,9))
  pre.bio12.com <- data.frame(pre.bio12.com)

  j=1
  for (j in 1:length(rich.com$X1)){

    ele.mean = apply(act.ele.com[j,], 1, mean)
    prec.mean = apply(pre.bio12.com[j,], 1, mean)
    
    rich.a = as.numeric(c(rich.com[j,]))
    prec.b = as.numeric(pre.bio12.com[j,])
    
    prec.rich.df = data.frame(rich.a = as.vector(rich.a),
                              prec.b = as.vector(prec.b))
    
    # all sites
    prec.rich.lm.tmp = summary(lm(rich.a ~ prec.b, prec.rich.df))
    table.prec.tmp = data.frame(ele = ele.mean,
                                pre.prec = prec.mean,
                                slope = prec.rich.lm.tmp$coefficients[2,1],
                                Ine = prec.rich.lm.tmp$coefficients[1,1],
                                p = prec.rich.lm.tmp$coefficients[2,4],
                                r2 = prec.rich.lm.tmp$adj.r.squared)
 
   if (j==1)
  {
    table.prec = table.prec.tmp
  }
  else{
    table.prec = rbind(table.prec, table.prec.tmp)
  }
} 
  
  if (i==1)
  {
    rich.prec.sensi = table.prec
  }
  else{
    rich.prec.sensi = rbind(rich.prec.sensi, table.prec)
  }
}






i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[1]
  rich.t <- micro.rich.tmp[which(vec.tmp$Elevation == ele.tmp),]
  
  rich.com <- t(combn(rich.t,9))
  rich.com <- data.frame(rich.com)
  rownames(rich.com)
  
  j=1
  for (j in 1:length(rownames(rich.com))){
    rich.a = as.numeric(c(rich.com[j,]))
  }
      
  
  
  ######## for mean ele #######
  act.ele.t <- Act.ele[which(vec.tmp$Elevation == ele.tmp),]
  
  act.ele.com <- t(combn(act.ele.t,9))
  act.ele.com <- data.frame(act.ele.com)
  rownames(act.ele.com)
  
  j=1
  for (j in 1:length(rownames(act.ele.com))){
    ele.mean = apply(act.ele.com[j,], 1, mean)
  }

  

  ####### for mean pre.bio12 ########
  pre.bio12.t <- pre.bio12[which(vec.tmp$Elevation == ele.tmp),]
  
  pre.bio12.com <- t(combn(pre.bio12.t,9))
  pre.bio12.com <- data.frame(pre.bio12.com)
  rownames(pre.bio12.com)
  
  j=1
  for (j in 1:length(rownames(act.ele.com))){
    prec.mean = apply(pre.bio12.com[j,], 1, mean)
    
    prec.b = as.numeric(pre.bio12.com[j,])
    
    prec.rich.df = data.frame(rich.a = as.vector(rich.a),
                              prec.b = as.vector(prec.b))
    
    # all sites
    prec.rich.lm.tmp = summary(lm(rich.a ~ prec.b, prec.rich.df))
    table.prec.tmp = data.frame(ele = ele.mean,
                                pre.prec = prec.mean,
                                slope = prec.rich.lm.tmp$coefficients[2,1],
                                Ine = prec.rich.lm.tmp$coefficients[1,1],
                                p = prec.rich.lm.tmp$coefficients[2,4],
                                r2 = prec.rich.lm.tmp$adj.r.squared)
  }
  
  

  
  ####### for mean pre.AI ########
  pre.AI.t <- pre.AI[which(vec.tmp$Elevation == ele.tmp),]
  
  pre.AI.com <- t(combn(pre.AI.t,9))
  pre.AI.com <- data.frame(pre.AI.com)
  rownames(pre.AI.com)
  
  AI.mean = apply(pre.AI.com[1,], 1, mean)
  
  AI.b = as.numeric(pre.AI.com[1,])
  
  AI.rich.df = data.frame(rich.a = as.vector(rich.a),
                          AI.b = as.vector(AI.b))
  
  # all sites
  AI.rich.lm.tmp = summary(lm(rich.a ~ AI.b, AI.rich.df))
  table.AI.tmp = data.frame(ele = ele.mean,
                            pre.AI = AI.mean,
                            slope = AI.rich.lm.tmp$coefficients[2,1],
                            Ine = AI.rich.lm.tmp$coefficients[1,1],
                            p = AI.rich.lm.tmp$coefficients[2,4],
                            r2 = AI.rich.lm.tmp$adj.r.squared)
  
  ####### for mean pre.SWC ########
  SWC.t <- SWC[which(vec.tmp$Elevation == ele.tmp),]
  
  SWC.com <- t(combn(SWC.t,9))
  SWC.com <- data.frame(SWC.com)
  rownames(SWC.com)
  
  SWC.mean = apply(SWC.com[1,], 1, mean)
  
  SWC.b = as.numeric(SWC.com[1,])
  
  SWC.rich.df = data.frame(rich.a = as.vector(rich.a),
                          SWC.b = as.vector(SWC.b))
  
  # all sites
  SWC.rich.lm.tmp = summary(lm(rich.a ~ SWC.b, SWC.rich.df))
  table.SWC.tmp = data.frame(ele = ele.mean,
                            pre.SWC = SWC.mean,
                            slope = SWC.rich.lm.tmp$coefficients[2,1],
                            Ine = SWC.rich.lm.tmp$coefficients[1,1],
                            p = SWC.rich.lm.tmp$coefficients[2,4],
                            r2 = SWC.rich.lm.tmp$adj.r.squared)
  
  
  ##########################################
  sensitivity.tmp = cbind(table.prec.tmp, table.AI.tmp, table.SWC.tmp)
  

  if (i==1)
  {
    sensitivity = sensitivity.tmp
  }
  else{
    sensitivity = rbind(sensitivity, sensitivity.tmp)
  }
}


sensitivity = data.frame(sensitivity)


act.ele.CV[1:5, 1:13]
dim(act.ele.CV)

act.ele.CV.correct = act.ele.CV[,c(1:3, 12:13)]

colnames(act.ele.CV.correct) = c("Ele", "pre.bio12", "pre.AI", 
                                 "Rich.CV", "Asy")
act.ele.CV.correct[1:5, 1:5]

ele.r = round(act.ele.CV.correct$Ele, 0)

CV.Asy = cbind(ele.r, act.ele.CV.correct)

CV.Asy[1:5, 1:6]

CV.Asy2 = CV.Asy[, c("ele.r", "Rich.CV", "Asy")]

CV.Asy3 = melt(CV.Asy2, id = ("ele.r")) # , "pre.bio12", "pre.AI"


#########  ele, CV and asy ########### 
CV.me = dcast(CV.Asy3, ele.r ~ variable, mean)
CV.me[1:3, 1:3]

CV.me.sd = dcast(CV.Asy3, ele.r ~ variable, sd)

# number of each elevation
hori.len <- aggregate(CV.Asy3$value, by=list(CV.Asy3$ele.r, CV.Asy3$variable), 
                      FUN=length)

CV.me1 = melt(CV.me, id = c("ele.r"))

# mean and se for hori
hori.all <- data.frame(CV.me1, se = (CV.me1$value)/sqrt(hori.len$x))


# hori.all$ele.r = factor(hori.all$ele.r)

ggplot(hori.all, aes(x=ele.r, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.4)+
  # geom_smooth(method = "glm", formula = y ~ poly(x, 7), 
  #             se = T, color = "red")+ 
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))

ggsave("../2.DX2013_145sites/figures/Fig.61.micro1314.rich.CV.9site.AI.pdf",
       width = 6.5, height = 3)


CV.all = subset(hori.all, hori.all$variable == "Rich.CV")
CV.all = CV.all[,-4]

library(segmented)
x1=CV.all$ele.r
y1=CV.all$value
Turb1<-glm(y1~x1,data=CV.all)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
summary(Turb1.seg)


pdf("../2.DX2013_145sites/figures/Fig.62.micro1314.rich.CV.9site.all.tipping.pdf",
    width = 6.5, height = 3)
# edge distance 
par(mar=c(4.5,4.5,4.5,4.5))
plot(x1,y1,xlab='Ele', ylab='Rich.CV',
     pch = 16, cex = 2, col = "#64C9CF")
plot(Turb1.seg, add=TRUE, link=FALSE, lwd=2, col="red",
     lty=1, conf.level=0.95, shade=F,)
lines(Turb1.seg, col="red", pch = 19, lwd = 2)
title(font.lab = 2, cex.lab = 2)
dev.off()


#### Asy ####
Asy.all = subset(hori.all, hori.all$variable == "Asy")
Asy.all = Asy.all[,-4]


x1=Asy.all$ele.r
y1=Asy.all$value
Turb1<-glm(y1~x1,data=Asy.all)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
summary(Turb1.seg)


pdf("../2.DX2013_145sites/figures/Fig.63.micro1314.rich.Asy.9site.all.tipping.pdf",
    width = 6.5, height = 3)
# edge distance 
par(mar=c(4.5,4.5,4.5,4.5))
plot(x1,y1,xlab='Ele', ylab='Asy',
     pch = 16, cex = 2, col = "#64C9CF")
plot(Turb1.seg, add=TRUE, link=FALSE, lwd=2, col="red",
     lty=1, conf.level=0.95, shade=F,)
lines(Turb1.seg, col="red", pch = 19, lwd = 2)
title(font.lab = 2, cex.lab = 2)
dev.off()


#########  pre.bio12, CV and asy ########### 
CV.Asy[1:5, 1:6]

prec.CV.Asy2 = CV.Asy[, c("pre.bio12", "Rich.CV", "Asy")]

prec.CV.Asy3 = melt(prec.CV.Asy2, id = ("pre.bio12")) 


prec.CV.me = dcast(prec.CV.Asy3, pre.bio12 ~ variable, mean)

prec.CV.me.sd = dcast(prec.CV.Asy3, pre.bio12 ~ variable, sd)

# number of each elevation
prec.hori.len <- aggregate(prec.CV.Asy3$value, by=list(prec.CV.Asy3$pre.bio12, 
                                                       prec.CV.Asy3$variable), 
                           FUN=length)

prec.CV.me1 = melt(prec.CV.me, id = c("pre.bio12"))

# mean and se for hori
prec.hori.all <- data.frame(prec.CV.me1, se = (prec.CV.me1$value)/sqrt(prec.hori.len$x))


ggplot(prec.hori.all, aes(x=pre.bio12, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.4)+
  # geom_smooth(method = "glm", formula = y ~ poly(x, 7), 
  #             se = T, color = "red")+ 
  facet_wrap( ~ variable, scales = "free", ncol = 1)+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))



#########  ele, CV and asy ########### 
CV.me = dcast(CV.Asy3, pre.AI ~ variable, mean)

CV.me.sd = dcast(CV.Asy3, pre.AI ~ variable, sd)

# number of each elevation
hori.len <- aggregate(CV.Asy3$value, by=list(CV.Asy3$pre.AI, CV.Asy3$variable), 
                      FUN=length)

CV.me1 = melt(CV.me, id = c("pre.AI"))
# mean and se for hori
hori.all <- data.frame(CV.me1, se = (CV.me1$value)/sqrt(hori.len$x))


ggplot(hori.all, aes(x=pre.AI, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.4)+
  # geom_smooth(method = "glm", formula = y ~ poly(x, 7), 
  #             se = T, color = "red")+ 
  facet_wrap( ~ variable, scales = "free", ncol = 1)+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))

