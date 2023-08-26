
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


# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")


load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")

cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

veg14 = read.csv("../2.DX2013_145sites/data/plant.2014.species.number.csv",
                 row.names = 1)

#### 2014 plant quadrat areas 1 m2 #####
#### standard with 2020 0.25 m2 ########
veg14 = round((veg14/4))

# na value equal 0
veg14[is.na(veg14)] <- 0

# to match the data 2014
sp36.45 = matrix(data=0, nrow = 110, ncol = 10, byrow = FALSE, dimnames = NULL)
colnames(sp36.45) = c("sp36", "sp37", "sp38", "sp39","sp40",
                      "sp41", "sp42", "sp43", "sp44","sp45")

veg13 = cbind(Veg.sp, sp36.45)


# combine veg data 2013 and 2014
veg1314 = rbind(veg13, veg14)

HDX032 = read.csv("../2.DX2013_145sites/data/plant-HDX032.csv")
colnames(HDX032)[1] = c("Group")
HDX032.tmp = melt(HDX032, id = c("Group"))

split <- strsplit(as.character(HDX032.tmp$value), "/")
split

species <- sapply(split, "[", 1)
height <- sapply(split, "[", 2)
Number <- sapply(split, "[", 3)

Number = data.frame(Number)
split2 <- strsplit(as.character(Number$Number), " ")
split2

Count <- sapply(split2, "[", 1)
cover <- sapply(split2, "[", 2)


HDX032taxa <- cbind(species, height, Count, cover)

HDX032last = na.omit(HDX032taxa)

HDX032last = t(HDX032last)

colnames(HDX032last) = HDX032last[1,]

HDX032last1 = t(HDX032last)

HDX032last1 = data.frame(HDX032last1)

HDX032last2 = as.data.frame(lapply(HDX032last1, as.numeric))

order = HDX032last2[order(HDX032last2$species), ]

order1 = order[,c(1,4)]

write.csv(order1, "../2.DX2013_145sites/data/order1.HDX032.cover.csv")

# library(tidyr)
# 
# order2<-spread(order1, species, cover)

cover2 = read.csv("../2.DX2013_145sites/data/order1.HDX032.cover2.csv")

corHDX032 = cor(cover2)

occor<-rcorr(cover2,type = "spearman")
occor = Hmisc::rcorr(as.matrix(as.data.frame(lapply(cover2, as.numeric))), type = "spearman") 

occor.r = occor$r 
occor.p = occor$P 
# Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
occor.p<-p.adjust(occor.p, method="BH")
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0 
diag(occor.r) <- 0 

sum(occor.r)/21

tr(2.13^occor.r)

tr(2.13^abs(occor.r))

tmp = occor.r[-1, -1]

tr(2.13^tmp)

tr(2.13^abs(tmp))

tmp2 = tmp*(-1)

tr(2.13^tmp2)

tr(2.13^abs(tmp2))


trl = tr(2.13^l)/tr(2.13^abs(l))
?graph_from_adjacency_matrix
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


tr = tr(2.13^igraph)/tr(2.13^abs(igraph))
K = tr(exp(igraph))/tr(exp(abs(igraph)))

K = tr(exp(occor.r))/tr(exp(abs(occor.r)))

igraph
tr(igraph)

m <- matrix(1:16,ncol=4)
m
tr(m)


library(igraph)

# create a simple signed two mode network
el <- matrix(c(1,"a",1,"b",1,"c",2,"a",2,"b"),ncol = 2,byrow = TRUE)
g <- graph_from_edgelist(el,directed = FALSE)
g

E(g)


E(g)$sign <- c(1,1,-1,1,-1)
V(g)$type <- c(FALSE,TRUE,TRUE,TRUE,FALSE)

# convert to unsigned two-mode network and project
# install.packages("signnet")
library(signnet)


as_unsigned_2mode(g, primary = TRUE)
l <- as_unsigned_2mode(g, primary = TRUE)

trl = tr(2.13^as.numeric(l))/tr(2.13^abs(l))


tr(2.13^as.numeric(l))

l2 = as.matrix(as.data.frame(lapply(v(l), as.numeric)))

p <- bipartite_projection(l,which="true")

# turn the unsigned projection back to a signed network
as_signed_proj(p)


HDX132me = melt(HDX032last2, id = c("species", "height"))


HDX032last3 = dcast(HDX132me, species ~ variable, sum)

HDX032last4 = cbind(HDX032last3, (HDX032last3$cover)/100)

rbind(HDX032last1, HDX032last1, HDX032last1)

b = cor(HDX032last2)

tmp = cor(veg1314)

cov = sum(as.numeric(HDX032last$cover))

