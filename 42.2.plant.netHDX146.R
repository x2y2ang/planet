
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(WGCNA)
library(igraph)
library(psych)
library(vegan)
library(FactoMineR)
library(plyr)
library(reshape2)
library(igraph)
library(impute)
library(GO.db)
library(preprocessCore)
library(AnnotationDbi)
library(psych)
library(Hmisc)
library(splitstackshape)


HDX146 = read.csv("../2.DX2013_145sites/data/plant-HDX146.csv")
colnames(HDX146)[1] = c("Group")
HDX146.tmp = melt(HDX146, id = c("Group"))

split <- strsplit(as.character(HDX146.tmp$value), "/")
split

species <- sapply(split, "[", 1)
height <- sapply(split, "[", 2)
Number <- sapply(split, "[", 3)

Number = data.frame(Number)
split2 <- strsplit(as.character(Number$Number), " ")
split2

Count <- sapply(split2, "[", 1)
cover <- sapply(split2, "[", 2)


HDX146taxa <- cbind(species, height, Count, cover)

HDX146last = na.omit(HDX146taxa)

HDX146last = t(HDX146last)

colnames(HDX146last) = HDX146last[1,]

HDX146last1 = t(HDX146last)

HDX146last1 = data.frame(HDX146last1)

HDX146last2 = as.data.frame(lapply(HDX146last1, as.numeric))

order = HDX146last2[order(HDX146last2$species), ]

order1 = order[,c(1,4)]

a = dcast(order1, cover ~ species)

b = dcast(order1, species ~ cover)

b = melt(order1, species ~ cover)

?spread
spread(order1, )

write.csv(order1, "../2.DX2013_145sites/data/order1.HDX146.cover.csv")

# library(tidyr)
# 
# order2<-spread(order1, species, cover)

cover2 = read.csv("../2.DX2013_145sites/data/order1.HDX146.cover2.csv")
corHDX146A = cor(cover2)
corHDX146 = cor((cover2/100))

occor<-rcorr(cover2,type = "spearman")
occor = Hmisc::rcorr(as.matrix(as.data.frame(lapply(cover2, as.numeric))), type = "spearman") 

occor.r = occor$r 
occor.p = occor$P 
# Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
occor.p<-p.adjust(occor.p, method="BH")
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0 
diag(occor.r) <- 0 

sum(occor.r)/21

occor.r = data.frame(occor.r)
upper = occor.r[upper.tri(occor.r)]
mean(((upper-mean(upper))/sd(upper))^4)

# HDX032 1.290421
# HDX145 4.080024


occor.r = data.frame(occor.r)
upper = occor.r[upper.tri(occor.r)]
mean(((upper-mean(upper))/sd(upper))^4)
# 1.290421

igraph = graph_from_adjacency_matrix(occor.r,mode="directed",weighted=TRUE,diag=FALSE)
igraph



igraph.weight = E(igraph)$weight
sum(igraph.weight>0)
sum(igraph.weight<0)

E.color = igraph.weight
E.color = ifelse(E.color>0, "red",ifelse(E.color<0, "blue","grey"))
E(igraph)$color = as.character(E.color)

set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


K = tr(exp(igraph))/tr(exp(abs(igraph)))


K = tr(exp(occor.r))/tr(exp(abs(occor.r)))

4/7
