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

### coverage 
plant.cov <- read.csv("../2.DX2013_145sites/data/last.plant1314.188sites.coverage.for.saiz.method1.csv",
                      row.names = 1)
plant.cov[is.na(plant.cov)] <- 0

# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")

load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")

cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

cli.pre.AI = cli1314.last[rownames(plant.cov), ]

# write.csv(cli.pre.AI, "../2.DX2013_145sites/data/cli.pre.AI188sites.csv")

dim(comm1314.samp)

# otu table relative abandance
rowSums(comm1314.samp)
RA.otu = (comm1314.samp/14399)*100
RA.otu[1:10, 1:10]

# genus relative abandance
Sp.ra = data.frame((colSums(RA.otu)/223))

colnames(Sp.ra)[1] = c("RA")

# reserve RA > 0.02%
Sp.ra.or1 = subset(Sp.ra, RA >= 0.02)

comm.tmp.RA = comm1314.samp[, rownames(Sp.ra.or1)]

# occurring in 60% of all samples 223*0.6 = 134
223*0.6 
dim(comm.tmp.RA)
comm.tmp.R1 = comm.tmp.RA
comm.tmp.R1[comm.tmp.R1 >0] <- 1
comm.mv <- comm.tmp.RA[, which(colSums(comm.tmp.R1) >= 134)]  
dim(comm.mv)


comm.tmp = comm.mv[rownames(plant.cov), ]


envi.tmp = envi1314[rownames(plant.cov), ]
colnames(envi.tmp)

### combine together ###
PMS = cbind(plant.cov, envi.tmp[, c(11:19)], comm.tmp)

PMS[1:8, 1:8]


ele <- unique(envi.tmp[,"Elevation"])
length(ele)

i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  otu.tmp <- comm.mv[which(envi.tmp$Elevation==ele.tmp),]
  occor = corr.test(otu.tmp, use="pairwise", method="spearman",
                    adjust="fdr", alpha=.05)
  occor.r = occor$r 
  occor.p = occor$p 
  occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0 
  igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",
                                       weighted=TRUE,diag=FALSE)
  fc = cluster_fast_greedy(igraph,weights =NULL)
  modularity = modularity(igraph,membership(fc))
  num.edges = length(E(igraph)) 
  num.vertices = length(V(igraph))
  connectance = edge_density(igraph,loops=FALSE)
  average.degree = mean(igraph::degree(igraph))
  average.path.length = average.path.length(igraph)
  clustering.coefficient = transitivity(igraph) 
  no.clusters = no.clusters(igraph)
  centralization.betweenness = centralization.betweenness(igraph)$centralization 
  centralization.degree = centralization.degree(igraph)$centralization
  
  net.attr = cbind(ele = ele[i],
                   modularity, num.edges, num.vertices, connectance, 
                   average.degree, average.path.length,
                   clustering.coefficient, no.clusters, 
                   centralization.betweenness,
                   centralization.degree)
  
  table.sig.tmp = data.frame(net.attr)
  
  if(i==1){
    table.sig = table.sig.tmp
  }
  else{
    table.sig = rbind(table.sig, table.sig.tmp)
  }
}

write.csv(table.sig, "../3.DX1314.Network/result/1.ele10.plant.soil.microbial.network.attribute.csv")

net.sig = read.csv("../3.DX1314.Network/result/1.ele10.plant.soil.microbial.network.attribute.csv")


### plot ###

net.sig1 = net.sig[order(net.sig$ele), ]
net.sig1 = net.sig1[, c(2:12)]


aridity = as.numeric(c("0.99", "0.94", "0.86",
                       "0.78", "0.69", "0.58",
                       "0.45", "0.30", "0.15",
                       "0.002"))

net.ai = cbind(aridity, net.sig1)

colnames(net.ai)[1:2] = c("AI", "Elevation")

net.melt = melt(net.ai, id = c("AI", "Elevation"))

## 4300 only 4 samples delete
net.melt1 = subset(net.melt, Elevation != 4300)
net.melt1$Elevation = factor(net.melt1$Elevation)


library(ggplot2)
library(RColorBrewer)

ggplot(net.melt1, aes(x = AI, y = value, color = Elevation))+
  geom_point(size = 3)+
  facet_wrap( ~ variable, scales ="free_y", ncol = 3)+
  stat_smooth(method="lm", color="black", formula=y~x, se = F)+
  stat_smooth(method="lm", color="red", formula=y~x+I(x^2), se = F)+
  scale_colour_manual(values = brewer.pal(12,"Paired"))+
  theme_bw()+
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))+  
  theme(strip.background = element_blank())+
  theme(text = element_text())+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(),
        axis.title.x=element_text(),
        axis.title.y=element_text())+
  xlab("Aridity index")+ylab("")


ggsave("../2.DX2013_145sites/figures/plant.soil.micro.network.attributes.pdf",
       height = 7,width = 7.5)


## all
net.melt$Elevation = factor(net.melt$Elevation)

ggplot(net.melt, aes(x = AI, y = value, color = Elevation))+
  geom_point(size = 3)+
  facet_wrap( ~ variable, scales ="free_y", ncol = 3)+
  stat_smooth(method="lm", color="black", formula=y~x, se = F)+
  stat_smooth(method="lm", color="red", formula=y~x+I(x^2), se = F)+
  scale_colour_manual(values = brewer.pal(12,"Paired"))+
  theme_bw()+
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))+  
  theme(strip.background = element_blank())+
  theme(text = element_text())+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(),
        axis.title.x=element_text(),
        axis.title.y=element_text())+
  xlab("Aridity index")+ylab("")


ggsave("../2.DX2013_145sites/figures/ELE10.plant.soil.micro.network.attributes.pdf",
       height = 7,width = 7.5)


