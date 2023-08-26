
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
library(reshape2)

comm1314 <- read.csv("../damxung145sites/data/comm_YXQ.csv", row.names = 1)
envi1314 <- read.csv("../damxung145sites/data/ENVI7.13.csv", row.names = 1)

taxa1314 <- data.frame(otu.id = rownames(comm1314),
                   taxonomy = comm1314$taxonomy)
comm1314 <- comm1314[, ! colnames(comm1314) %in% "taxonomy"]; dim(comm1314)
comm1314 <- t(comm1314); dim(comm1314)

# delect HDX111 
comm1314.de <- comm1314[rownames(envi1314),]
dim(comm1314.de)

#split character to mitiple columns
split <- strsplit(as.character(taxa1314$taxonomy), ";")
split

Kingdom <- sapply(split, "[", 1)
Phylum <- sapply(split, "[", 2)
Class <- sapply(split, "[", 3)
Order <- sapply(split, "[", 4)
Family <- sapply(split, "[", 5)
Genus <- sapply(split, "[", 6)

taxa1314 <- cbind(taxa1314, Kingdom, Phylum, Class, Order, Family, Genus)


# summary comm
summary(rowSums(comm1314.de))

# rarefactions use Min 14399 of comm_Bacteria 
comm1314.samp <- rrarefy(comm1314.de, 14399)

# The otu of the column is greater than zero
comm1314.samp <- comm1314.samp[, colSums(comm1314.samp) >0]
# dim(comm.samp)

# Match id to the sample point in comm.samp

rownames(taxa1314) <- taxa1314[,1]
taxa1314.samp <- taxa1314[colnames(comm1314.samp),]
# dim(taxa.Bac.samp)

# save.data
save(envi1314, taxa1314.samp, comm1314.samp,  
     file = "../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")


load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")


source("../2.DX2013_145sites/function/diversity.functions.R")
dim(comm1314.samp)
micro1314 = quick.diversity(comm1314.samp)

vec13.tmp = envi[, c("Act.ele", "Elevation", "bio12", "bio1", "Aridity.index",
                     "GTmin")] 

write.csv(vec13.tmp, "../2.DX2013_145sites/data/envi13.bio12.AI.GT.csv")


colnames(envi1314)
vec1314.tmp = envi1314[, c("Act.ele", "Elevation", "bio12", "bio1")] 

write.csv(vec1314.tmp, "../2.DX2013_145sites/data/envi1314.bio12.bio1.csv")

cli1314 = read.csv("../2.DX2013_145sites/data/envi1314.bio12.bio1.add.AI.GT.csv",
                   row.names = 1)

# predit elevation
quick.precdict <- function(var) {
  for(i in 3:6){
    x <- var[,2]
    y <- var[,i]
    lm.tmp <- lm(y ~ I(x^2) + x)
    Elevation <- var[,1]
    bio.new <- data.frame(x = Elevation)
    bio.tmp <- predict(lm.tmp, bio.new)
    if(i==3)
    {
      bio.pre = bio.tmp
    }
    else{
      bio.pre = data.frame(bio.pre, bio.tmp)
    }
  }
  return(bio.pre)
}

act.cli1314 <- quick.precdict(cli1314)

cli.act.ele = cbind(cli1314, act.cli1314)

colnames(cli.act.ele)[7:10] = c("pre.bio12", "pre.bio1", "pre.AI", "pre.GT")

write.csv(cli.act.ele, "../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.csv")


cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                   row.names = 1)

# for sd
rowVars <- function(x, na.rm=F) {
  # Vectorised version of variance filter
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}


################    loop for 10 eles    ######################
vec.tmp = cbind(cli1314.last, micro1314, envi1314)

write.csv(vec.tmp, "../2.DX2013_145sites/data/micro.diversity.1314.bio12.bio1.AI.GT.csv")



micro.rich.tmp = data.frame(vec.tmp[, "comm.richness"])
Act.ele = data.frame(vec.tmp[, c("Act.ele")])
pre.bio12 = data.frame(vec.tmp[, c("pre.bio12")])
pre.AI = data.frame(vec.tmp[, c("pre.AI")])


ele = unique(vec.tmp[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  rich.t <- micro.rich.tmp[which(vec.tmp$Elevation == ele.tmp),]
  
  rich.com <- t(combn(rich.t,9))
  rich.com <- data.frame(rich.com)
  rownames(rich.com)
  
  rich.mean = apply(rich.com, 1, mean)
  rich.sd = apply(rich.com, 1, sd)
  
  CV.tmp= rich.sd/rich.mean
  
  # Asy strong functions for row
  rich.max = apply(rich.com, 1, max)
  rich.min = apply(rich.com, 1, min)
  
  Asy = ((rich.max - rich.mean)/(rich.mean - rich.min))
  
  
  ######## for mean ele #######
  act.ele.t <- Act.ele[which(vec.tmp$Elevation == ele.tmp),]
  
  act.ele.com <- t(combn(act.ele.t,9))
  act.ele.com <- data.frame(act.ele.com)
  rownames(act.ele.com)

  act.ele.mean = apply(act.ele.com, 1, mean)
  
  
  ####### for mean pre.bio12 ########
  pre.bio12.t <- pre.bio12[which(vec.tmp$Elevation == ele.tmp),]
  
  pre.bio12.com <- t(combn(pre.bio12.t,9))
  pre.bio12.com <- data.frame(pre.bio12.com)
  rownames(pre.bio12.com)

  pre.bio12.mean = apply(pre.bio12.com, 1, mean)
  
  
  ####### for mean pre.AI ########
  pre.AI.t <- pre.AI[which(vec.tmp$Elevation == ele.tmp),]
  
  pre.AI.com <- t(combn(pre.AI.t,9))
  pre.AI.com <- data.frame(pre.AI.com)
  rownames(pre.AI.com)
  
  pre.AI.mean = apply(pre.AI.com, 1, mean)
  
  
  ##########################################
  act.ele.CV.tmp = cbind(act.ele.mean, pre.bio12.mean, pre.AI.mean, 
                         CV.tmp, Asy)
  
  colnames(act.ele.CV.tmp) = c("Ele", "pre.bio12", "pre.AI", 
                               "Rich.CV", "Asy")
  
  if (i==1)
  {
    act.ele.CV = act.ele.CV.tmp
  }
  else{
    act.ele.CV = rbind(act.ele.CV, act.ele.CV.tmp)
  }
}

act.ele.CV = data.frame(act.ele.CV)

act.ele.CV[1:50, 1:5]

ele.r = round(act.ele.CV$Ele, 0)

CV.Asy = cbind(round(act.ele.CV$Ele, 0), round(act.ele.CV$pre.bio12, 0),
               round(act.ele.CV$pre.AI, 3), act.ele.CV[, c("Rich.CV", "Asy")])

CV.Asy[1:50, 1:5]

colnames(CV.Asy) = c("Ele", "pre.bio12", "pre.AI", 
                             "Rich.CV", "Asy")


CV.Asy.md <- melt(CV.Asy, id=c("Ele", "pre.bio12", "pre.AI"))
CV.Asy.md[1:10, 1:5]


CV.me = dcast(CV.Asy.md, Ele  ~ variable, mean)


CV.me.sd = dcast(CV.Asy.md, Ele  ~ variable, sd)


# number of each elevation
hori.len <- aggregate(CV.Asy.md$value, by=list(CV.Asy.md$Ele,  CV.Asy.md$variable), 
                      FUN=length)


CV.me1 = melt(CV.me, id = c("Ele"))


# mean and se for hori
hori.all <- data.frame(CV.me1, se = (CV.me1$value)/sqrt(hori.len$x))

hori.all1 = na.omit(hori.all)

write.csv(hori.all, "../2.DX2013_145sites/data/micro.rich.1314.CV.Asy.elevation.csv")


ggplot(hori.all, aes(x=Ele, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.4)+
  geom_smooth(method = "glm", formula = y ~ poly(x, 2),
              se = T, color = "red")+
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


################# prec  ############### 
CV.Asy.md <- melt(CV.Asy, id=c("Ele", "pre.bio12", "pre.AI"))

prec.CV.me = dcast(CV.Asy.md, pre.bio12  ~ variable, mean)


prec.CV.me.sd = dcast(CV.Asy.md, pre.bio12  ~ variable, sd)


# number of each elevation
prec.hori.len <- aggregate(CV.Asy.md$value, by=list(CV.Asy.md$pre.bio12,  CV.Asy.md$variable), 
                      FUN=length)


prec.CV.me1 = melt(prec.CV.me, id = c("pre.bio12"))


# mean and se for hori
prec.hori.all <- data.frame(prec.CV.me1, se = (prec.CV.me1$value)/sqrt(prec.hori.len$x))


write.csv(prec.hori.all, "../2.DX2013_145sites/data/micro.rich.1314.CV.Asy.precipitation.csv")


ggplot(prec.hori.all, aes(x=pre.bio12, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.4)+
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))

ggsave("../2.DX2013_145sites/figures/Fig.62.prec.micro1314.rich.CV.9site.AI.pdf",
       width = 6.5, height = 3)




CV.all = subset(prec.hori.all, prec.hori.all$variable == "Rich.CV")
CV.all = CV.all[,-4]

library(segmented)
x1=CV.all$pre.bio12
y1=CV.all$value
Turb1<-glm(y1~x1,data=CV.all)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
summary(Turb1.seg)


pdf("../2.DX2013_145sites/figures/Fig.63.prec.micro1314.rich.CV.9site.all.tipping.pdf",
    width = 5, height = 3)
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
Asy.all = subset(prec.hori.all, prec.hori.all$variable == "Asy")
Asy.all = Asy.all[,-4]


x1=Asy.all$pre.bio12
y1=Asy.all$value
Turb1<-glm(y1~x1,data=Asy.all)
Turb1.seg<-segmented(Turb1,seg.Z=~x1)
summary(Turb1.seg)


pdf("../2.DX2013_145sites/figures/Fig.64.prec.micro1314.rich.Asy.9site.all.tipping.pdf",
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



#########  pre.AI, CV and asy ########### 
CV.Asy[1:5, 1:5]

AI.CV.Asy2 = CV.Asy[, c("pre.AI", "Rich.CV", "Asy")]

AI.CV.Asy3 = melt(AI.CV.Asy2, id = ("pre.AI")) 

memory.limit(size = 35000)
AI.CV.me = dcast(AI.CV.Asy3, pre.AI ~ variable, mean)

AI.CV.me.sd = dcast(AI.CV.Asy3, pre.AI ~ variable, sd)


# number of each elevation
AI.hori.len <- aggregate(AI.CV.Asy3$value, by=list(AI.CV.Asy3$pre.AI, 
                          AI.CV.Asy3$variable), FUN=length)

AI.CV.me1 = melt(AI.CV.me, id = c("pre.AI"))

# mean and se for hori
AI.hori.all <- data.frame(AI.CV.me1, se = (AI.CV.me1$value)/sqrt(AI.hori.len$x))

write.csv(AI.hori.all, "../2.DX2013_145sites/data/micro.rich.1314.CV.Asy.aridity.csv")

AI.hori.all = read.csv("../2.DX2013_145sites/data/micro.rich.1314.CV.Asy.aridity.csv")

ggplot(AI.hori.all, aes(x=pre.AI, y=value))+ 
  geom_point(size = 3, color = "#54BAB9")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), color = "#54BAB9")+
  stat_smooth(method="loess", color = "#4B6587", span = 0.5)+
  facet_wrap( ~ variable, scales = "free", ncol = 2)+ 
  xlab("") + ylab("") + # Set axis labels
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid=element_blank()) +
  # theme(strip.text.y = element_text(angle = 180))+
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))

pdf("../2.DX2013_145sites/figures/Fig.66.Aridity.micro1314.rich.Asy.9site.pdf",
      width = 6.5, height = 3)
