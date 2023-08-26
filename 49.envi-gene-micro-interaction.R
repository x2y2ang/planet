rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")


library(Hmisc)
library(reshape2)
library(igraph)


####??ȡ????????
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


# load data
# otu table relative abandance
rowSums(comm1314.samp)
RA.otu = (comm1314.samp/14399)*100
RA.otu[1:10, 1:10]

# genus relative abandance
Sp.ra = data.frame((colSums(RA.otu)/223))

rownames(Sp.ra) = colnames(comm1314.samp)
dim(comm1314.samp)

colnames(Sp.ra)[1] = c("RA")


# reserve RA > 0.02%
Sp.ra.or1 = subset(Sp.ra, RA >= 0.02)

comm.tmp.RA = comm1314.samp[, rownames(Sp.ra.or1)]

# occurring in 60% of all samples 223*0.6 = 134
dim(comm.tmp.RA)
comm.tmp.R1 = comm.tmp.RA
comm.tmp.R1[comm.tmp.R1 >0] <- 1
comm.mv <- comm.tmp.RA[, which(colSums(comm.tmp.R1) >= 134)]  
dim(comm.mv)


# use 6 for plant abundance data and as majority selection criterion
veg.tmp = t(veg1314)
veg.mv <- veg.tmp[rowSums(veg.tmp)>= 6,]

# occurring in 30% of all samples 188*0.3 = 56
dim(veg.tmp)
veg.R1 = veg.tmp
veg.R1[veg.R1 >0] <- 1
veg.mv <- veg.tmp[which(rowSums(veg.R1) >= 56), ]  
dim(veg.mv)


veg.shift = t(veg.mv)
micro188 = comm.mv[rownames(veg.shift), ]
envi188 = envi1314[rownames(veg.shift),]
envi188.select = envi188[, c(11:19)]


#### micro-plant spearman ####
genus_args_corr <- rcorr(as.matrix(micro188), as.matrix(veg.shift), type = 'spearman')

# Matrix of correlation coefficient r-values and significance p-values
r <- genus_args_corr$r
p <- genus_args_corr$P

# Only the correlation coefficient of microbial abundance-plant abundance is retained
# Remove the correlation coefficient between microbes-microbes, plant-plant
r <- r[colnames(micro188),colnames(veg.shift)]
p <- p[colnames(micro188),colnames(veg.shift)]

# Threshold filtering
# Remove the relationship with the spearman correlation coefficient lower than 0.7, that is, r>=0.7
# In this mode, we must pay attention to whether the choice of negative value is appropriate, because in many cases, the negative correlation may be meaningless
r[abs(r) < 0.2] <- 0

#ѡȡ?????? p ֵС?? 0.05 ??????ϵ?????? p<0.05
p <- p.adjust(p, method = 'BH')    #??ѡ p ֵУ????????ʹ?? BH ??У?? p ֵ
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

# Keep the data according to the above filtered r-values and p-values
z1 <- r * p
z1 <- data.frame(z1, check.names = FALSE)

#### microbes and envi188.select spearman correlation
micro188_envi188.select_corr <- rcorr(as.matrix(micro188), as.matrix(envi188.select), type = 'spearman')


r <- micro188_envi188.select_corr$r
p <- micro188_envi188.select_corr$P


r <- r[colnames(micro188),colnames(envi188.select)]
p <- p[colnames(micro188),colnames(envi188.select)]

r[abs(r) < 0.2] <- 0

p <- p.adjust(p, method = 'BH')    
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0


z2 <- r * p
z2 <- data.frame(z2, check.names = FALSE)

#### plant and envi spearman correlation
envi188.select_veg.shift_corr <- rcorr(as.matrix(envi188.select), as.matrix(veg.shift), type = 'spearman')

r <- envi188.select_veg.shift_corr$r
p <- envi188.select_veg.shift_corr$P

r <- r[colnames(envi188.select),colnames(veg.shift)]
p <- p[colnames(envi188.select),colnames(veg.shift)]

r[abs(r) < 0.15] <- 0


p <- p.adjust(p, method = 'BH')  
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0


z3 <- r * p
z3 <- data.frame(z3, check.names = FALSE)

#### Then convert to a symmetric matrix for subsequent conversion to a network in igraph
len <- ncol(micro188) + ncol(veg.shift) + ncol(envi188.select)
z <- data.frame(matrix(data=NA, nrow = len, ncol = len, byrow = FALSE, dimnames = NULL))
rownames(z) <- c(colnames(micro188), colnames(veg.shift), colnames(envi188.select))
colnames(z) <- c(colnames(micro188), colnames(veg.shift), colnames(envi188.select))

z[colnames(micro188), colnames(veg.shift)] <- z1[colnames(micro188), colnames(veg.shift)]
z[colnames(micro188), colnames(envi188.select)] <- z2[colnames(micro188), colnames(envi188.select)]
z[colnames(envi188.select), colnames(veg.shift)] <- z3[colnames(envi188.select), colnames(veg.shift)]

z[colnames(veg.shift), colnames(micro188)] <- t(z[colnames(micro188), colnames(veg.shift)])
z[colnames(envi188.select), colnames(micro188)] <- t(z[colnames(micro188), colnames(envi188.select)])
z[colnames(veg.shift), colnames(envi188.select)] <- t(z[colnames(envi188.select), colnames(veg.shift)])

z[colnames(micro188), colnames(micro188)] <- 0
z[colnames(envi188.select), colnames(envi188.select)] <- 0
z[colnames(veg.shift), colnames(veg.shift)] <- 0

write.table(z, 'micro188_veg.shift_envi188.select_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

write.table(z, '../2.DX2013_145sites/data/micro188_veg.shift_envi188.select_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

write.table(z, '../2.DX2013_145sites/data/micro188_veg.shift_envi188.select_corr.matrix.csv')

# Convert adjacency matrix to adjacency list for igraph network
# Construct a weighted undirected network, and the weights represent the spearman correlation coefficient between microbial-plant-envi factors
g <- graph_from_adjacency_matrix(as.matrix(z), 
        weighted = TRUE, mode = 'undirected')
g

plot(g)


# Orphan node deletion (delete nodes with degree 0)
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#In this mode, the edge weight represents the correlation coefficient
#Because the weight is usually a positive value, it is best to take an absolute value, and the correlation coefficient re-copy a column
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)


plot(g)

##Network file output, omitted, adjacency matrix, edge list, etc. can be obtained by referring to the above
#For example, gml format, which can be opened and edited visually with cytoscape software
write.graph(g, '../2.DX2013_145sites/data/57.microbial-plant-envinetwork.gml', format = 'gml')


memory.limit(size = 35000)

knn=(knn(simplify(g,remove.multiple = FALSE, remove.loops = TRUE))$knn)


#construct sub-networks for each sample and calculating each of the average degree and clustering coefficient
sample.out=matrix(0,ncol(veg.shift),3)
for(i in 1:ncol(veg.shift))
{
  g1 <- induced_subgraph(g, v=id[[i]])
  knn=(knn(simplify(g1,remove.multiple = FALSE, remove.loops = TRUE))$knn)
  knn=data.matrix(knn)
  avg.dgree=colMeans(knn,na.rm = TRUE)
  cluster=transitivity(g1)
  #cluster=mean(transitivity(g1, type =  "global", vids =NULL, weights = NULL, isolates = "NaN"))
  #modularity=modularity(g1, membership=id[[i]])
  wtc <- cluster_walktrap(g1)
  modularity=modularity(g1, membership(wtc))
  sample.out[i,1]=avg.dgree
  sample.out[i,2]=cluster
  sample.out[i,3]=modularity
}

rownames(sample.out)=colnames(comm)
colnames(sample.out)=c("average degree","clustering coefficient", "modularity")
sample.out

write.csv(sample.out,"../Desktop/Net-YXQ/3.1.network property for each sample.csv")


# ģ???? modularity
fc = cluster_fast_greedy(igraph,weights =NULL)# cluster_walktrap cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
modularity = modularity(igraph,membership(fc))
# ????ģ??Ϊ?ڵ???ɫ
comps = membership(fc)
colbar = rainbow(max(comps))
V(igraph)$color = colbar[comps] 

set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))



# ?Ƿ?ȥ????��???㣬?????Լ?ʵ??????
# remove isolated nodes????ȥ????????otu?????????Ե?otu ??ʡ?ԣ?ǰ?ھ????Ѵ?????
bad.vs = V(igraph)[degree(igraph) == 0]
igraph = delete.vertices(igraph, bad.vs)
igraph

# ??igraph weight???Ը?ֵ??igraph.weight
igraph.weight = E(igraph)$weight

# ??ͼǰȥ??igraph??weightȨ?أ???Ϊ??ͼʱĳЩlayout???ܵ???Ӱ??
E(igraph)$weight = NA


# ?򵥳?ͼ
# ?趨??????????????????ͼ????ͬһ??????????????????֤ǰ????ͼ??״????Ӧ
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

# ????????????ʱ??weighted=NULL,?˲??費??ͳ??
sum(igraph.weight>0)# number of postive correlation
sum(igraph.weight<0)# number of negative correlation

# set edge color??postive correlation ?趨Ϊred, negative correlation?趨Ϊblue
E.color = igraph.weight
E.color = ifelse(E.color>0, "red",ifelse(E.color<0, "blue","grey"))
E(igraph)$color = as.character(E.color)

# ?ı?edge??ɫ????ͼ
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


# ?????趨edge?Ŀ? ??set edge width?????罫????ϵ????edge width??��
E(igraph)$width = abs(igraph.weight)*4

# ?ı?edge???Ⱥ???ͼ
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


# ????OTUע????Ϣ???????൥Ԫ?ͷ???
# ????????????vertices size, vertices color��????????ά?ȵ?????
# ע??otu_pro.txt?ļ?Ϊ???????????????ݣ?????????ͼ???ܲ????????ض???ģʽ?????ɡ?
otu_pro = read.table("otu_pro.txt",head=T,row.names=1)
# set vertices size
igraph.size = otu_pro[V(igraph)$name,] # ɸѡ??ӦOTU????
igraph.size1 = log((igraph.size$abundance)*100) # ԭʼ??????ʲô??Ϊʲô*100??ȡe????
V(igraph)$size = igraph.size1

# set vertices color
igraph.col = otu_pro[V(igraph)$name,]
levels(igraph.col$phylum)
levels(igraph.col$phylum) = c("green","deeppink","deepskyblue","yellow","brown","pink","gray","cyan","peachpuff") # ֱ???޸?levles????��ֵȫ????Ӧ?滻
V(igraph)$color = as.character(igraph.col$phylum)

set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
