
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(vegan)
library(reshape2)
library(plyr)
library(dplyr)
library(openxlsx)
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
library(ieggr)
library(RColorBrewer)


plant.RA <- read.csv("../2.DX2013_145sites/data/last.plant1314.188sites.species.number.for.saiz.method.csv",
                     row.names = 1)

plant.RA[is.na(plant.RA)] <- 0


occor = Hmisc::rcorr(as.matrix(as.data.frame(lapply(plant.RA, as.numeric))), type = "spearman") 

occor.r = occor$r 
occor.p = occor$P 
# Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
occor.p<-p.adjust(occor.p, method="BH")
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0 
diag(occor.r) <- 0 


### coverage 
plant.cov <- read.csv("../2.DX2013_145sites/data/last.plant1314.188sites.coverage.for.saiz.method1.csv",
                     row.names = 1)

plant.cov[is.na(plant.cov)] <- 0

occor.cov = Hmisc::rcorr(as.matrix(as.data.frame(lapply(plant.cov, as.numeric))), type = "spearman") 

cov.r = occor.cov$r 
cov.p = occor.cov$P 
# Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
cov.p<-p.adjust(cov.p, method="BH")
cov.r[cov.p>0.05] = 0 # |abs(cov.r)<0.6
diag(cov.r) <- 0 

igraph = graph_from_adjacency_matrix(cov.r,mode="undirected",weighted=TRUE,diag=FALSE)
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


# network property
# 边数量 The size of the graph (number of edges)
num.edges = length(E(igraph)) # length(curve_multiple(igraph))
num.edges
# 顶点数量 Order (number of vertices) of a graph
num.vertices = length(V(igraph))# length(diversity(igraph, weights = NULL, vids = V(igraph)))
num.vertices
# 连接数(connectance) 网络中物种之间实际发生的相互作用数之和（连接数之和）占总的潜在相互作用数（连接数）的比例，可以反映网络的复杂程度
connectance = edge_density(igraph,loops=FALSE)# 同 graph.density;loops如果为TRUE,允许自身环（self loops即A--A或B--B）的存在
connectance
# 平均度(Average degree)
average.degree = mean(igraph::degree(igraph))# 或者为2M/N,其中M 和N 分别表示网络的边数和节点数。
average.degree
# 平均路径长度(Average path length)
average.path.length = average.path.length(igraph) # 同mean_distance(igraph) # mean_distance calculates the average path length in a graph
average.path.length

# 聚集系数(Clustering coefficient)：分局域聚类系数和全局聚集系数，是反映网络中节点的紧密关系的参数，也称为传递性。整个网络的全局聚集系数C表征了整个网络的平均的“成簇性质”。
clustering.coefficient = transitivity(igraph) 
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
# 介数中心性(Betweenness centralization)
centralization.betweenness = centralization.betweenness(igraph)$centralization 
centralization.betweenness
# 度中心性(Degree centralization)
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree



### sunnetwork
cor.otu = cor(plant.cov, method = "pearson")
dat=cor.otu

rownames(dat)

#check and match otu names in above two files
# ??match.name
match.rowname=match.name(both.list = list(dat.a=dat),
                         cn.list = list(comm.a=plant.cov))
comm.b=match.rowname$comm.a
dat.b=match.rowname$dat.a
dim(comm.b); dim(dat.b)
#transfer dat into matrix for igraph
m=as.matrix(dat)

#unique otu names appeared in each sample
dim(plant.cov)
comm=t(plant.cov)
# ?sapply
id=sapply(1:ncol(comm), function(i){rownames(comm[which(comm[,i]!=0),])})

#construct whole network with all samples
memory.limit(size = 35000)
gn=graph.adjacency(m, mode="undirected", weighted=NULL) # this will create an 'igraph object'
gn

#construct sub-networks for each sample and calculating each of the average degree and clustering coefficient
sample.out=matrix(0,ncol(comm),9)
for(i in 1:ncol(comm))
{
  g1 <- induced_subgraph(gn, v=id[[i]])
  knn=(knn(simplify(g1,remove.multiple = FALSE, remove.loops = TRUE))$knn)
  knn=data.matrix(knn)
  avg.dgree=colMeans(knn,na.rm = TRUE)
  cluster=transitivity(g1)
  #cluster=mean(transitivity(g1, type =  "global", vids =NULL, weights = NULL, isolates = "NaN"))
  #modularity=modularity(g1, membership=id[[i]])
  wtc <- cluster_walktrap(g1)
  modularity=modularity(g1, membership(wtc))
  
  num.edges = length(E(g1)) 
  num.vertices = length(V(g1))
  connectance = edge_density(g1,loops=FALSE)
  average.path.length = average.path.length(g1) 
  clustering.coefficient = transitivity(g1, type = c("global") )
  centralization.betweenness = centralization.betweenness(g1)$centralization 
  centralization.degree = centralization.degree(g1)$centralization
  
  sample.out[i,1]=avg.dgree
  sample.out[i,2]=cluster
  sample.out[i,3]=modularity
  sample.out[i,4]=num.edges
  sample.out[i,5]=num.vertices
  sample.out[i,6]=connectance
  sample.out[i,7]=average.path.length
  sample.out[i,8]=centralization.betweenness
  sample.out[i,9]=centralization.degree
  
}

rownames(sample.out)=colnames(comm)
colnames(sample.out)=c("average degree","clustering coefficient",
                       "modularity", "num.edges", "num.vertices",
                       "connectance", "average.path.length", 
                       "centralization.betweenness", "centralization.degree")

sample.out = data.frame(sample.out)



cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

cli188 = cli1314.last[rownames(plant.cov),]


subnet.plant = cbind(sample.out[, c(3:6)], cli188[, c(2:3)])


write.csv(subnet.plant, "../2.DX2013_145sites/data/subnetwork.for.plant1314.188site.csv")


veg.net.me <- melt(subnet.plant, id=c("Elevation", "pre.AI"))
veg.net.me$Elevation = factor(veg.net.me$Elevation)

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gp <- unique(veg.net.me[,"variable"])
length(gp)

resil <- as.data.frame(veg.net.me[,c("value")])
colnames(resil) <- c("resil")

AI <- as.data.frame(veg.net.me[,c("pre.AI")])

library(vegan)
library(MASS)
library(adespatial)
library(OTUtable)
library(reshape2)
library(ggplot2)

for (i in 1:length(gp)){
  gp.tmp <- gp[i]
  gp.t <- resil[which(veg.net.me$variable == gp.tmp),]
  gp.t = data.frame(gp.t)
  
  AI.t <- AI[which(veg.net.me$variable == gp.tmp),]
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
    env.envi = rbind(env.envi, env.gp.tem)
  }
}


lm.rich.all = env.envi
# p.value > 0.05
lm.rich.all$sig = ifelse(lm.rich.all$p > 0.05, "F","T")

write.csv(lm.rich.all, "../2.DX2013_145sites/data/lm.plant1314.188sites.subnetwork.csv")

ggplot(veg.net.me, aes(pre.AI, value, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  facet_wrap( ~ variable,nrow = 2, scales = "free")+
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
  xlab("Aridity")+ylab("")


ggsave("../2.DX2013_145sites/figures/plant1314.coverage.subnetwork.pdf",
       width = 6.5, height = 5)



