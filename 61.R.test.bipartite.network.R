
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("phyloseq")

library(phyloseq)

# remotes::install_github("taowenmicro/ggClusterNet")
library(ggClusterNet)

#---细菌和真菌OTU网络-域网络-二分网络#-------
# 仅仅关注细菌和真菌之间的相关，不关注细菌内部和真菌内部相关

Envnetplot<- paste("./16S_ITS_network",sep = "")
dir.create(Envnetplot)

data(psITS)

#--细菌和真菌ps对象中的map文件要一样
ps16s.otu_table = data.frame(phyloseq::otu_table(ps16s))
ps16s.tax_table = data.frame(phyloseq::tax_table(ps16s))
ps16s.sample_data = data.frame(phyloseq::sample_data(ps16s))
ps16s.refseq = data.frame(phyloseq::refseq(ps16s))
ps16s.phy_tree = phyloseq::phy_tree(ps16s)

psITS.otu_table = data.frame(phyloseq::otu_table(psITS))
psITS.tax_table = data.frame(phyloseq::tax_table(psITS))
psITS.sample_data = data.frame(phyloseq::sample_data(psITS))
psITS.refseq = data.frame(phyloseq::refseq(psITS))
psITS.phy_tree = phyloseq::phy_tree(psITS)

merge.otu_table = data.frame(phyloseq::otu_table(ps.merge))
merge.tax_table = data.frame(phyloseq::tax_table(ps.merge))


ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps16s,
                                       psITS = psITS,
                                       N16s = 500,
                                       NITS = 500
)

ps.merge

map =  phyloseq::sample_data(ps.merge)

head(map)
map$Group = "one"
phyloseq::sample_data(ps.merge) <- map

# data =  data.frame(ID = c("fun_ASV_205","fun_ASV_316","fun_ASV_118"),
#                    c("Coremode","Coremode","Coremode"))

# source("F:/Shared_Folder/Function_local/R_function/my_R_packages/ggClusterNet/R/corBionetwork2.R")

# data1 = env
# envRDA.s = vegan::decostand(envRDA,"hellinger")
# data1[,-1] = envRDA.s
# 
# Gru = data.frame(ID = colnames(env)[-1],group = "env" )
# head(Gru)
# 

#install.packages("WGCNA")
library(WGCNA)
library(igraph)
#--导入所需R包#-------

library(network)
library(sna)
library(tidyverse)
# install.packages("tidyfst")
library(tidyfst)

result <- corBionetwork(ps = ps.merge,
                        N = 0,
                        lab = data,
                        r.threshold = 0.8, # 相关阈值
                        p.threshold = 0.05,
                        group = "Group",
                        # env = data1, # 环境指标表格
                        # envGroup = Gru,# 环境因子分组文件表格
                        # layout = "fruchtermanreingold",
                        path = Envnetplot,# 结果文件存储路径
                        fill = "Phylum", # 出图点填充颜色用什么值
                        size = "igraph.degree", # 出图点大小用什么数据
                        scale = TRUE, # 是否要进行相对丰度标准化
                        bio = TRUE, # 是否做二分网络
                        zipi = F, # 是否计算ZIPI
                        step = 100, # 随机网络抽样的次数
                        width = 12,
                        label = TRUE,
                        height = 10,
                        big = TRUE,
                        select_layout = TRUE,
                        layout_net = "model_maptree2",
                        clu_method = "cluster_fast_greedy"
                        
                        
)

tem <- model_maptree(cor =result[[5]],
                     method = "cluster_fast_greedy",
                     seed = 12
)
node_model = tem[[2]]
head(node_model)

p = result[[1]]
p
# 全部样本网络参数比对
data = result[[2]]

plotname1 = paste(Envnetplot,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 20,height = 19)
tablename <- paste(Envnetplot,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)
tablename <- paste(Envnetplot,"/node_model_imformation",".csv",sep = "")
write.csv(node_model,tablename)

tablename <- paste(Envnetplot,"/nodeG_plot",".csv",sep = "")
write.csv(result[[4]],tablename)
tablename <- paste(Envnetplot,"/edge_plot",".csv",sep = "")
write.csv(result[[3]],tablename)
tablename <- paste(Envnetplot,"/cor_matrix",".csv",sep = "")
write.csv(result[[5]],tablename)




library(phyloseq)
library(ggClusterNet)
library(tidyverse)
library(Biostrings)

metadata = read.delim("https://raw.githubusercontent.com/taowenmicro/R-_function/main/metadata.tsv",row.names = 1)
otutab = read.delim("https://raw.githubusercontent.com/taowenmicro/R-_function/main/otutab.txt", row.names=1)
taxonomy = read.table("https://raw.githubusercontent.com/taowenmicro/R-_function/main/taxonomy.txt", row.names=1)
# tree  = read_tree("https://raw.githubusercontent.com/taowenmicro/R-_function/main/otus.tree")
# rep = readDNAStringSet("https://raw.githubusercontent.com/taowenmicro/R-_function/main/otus.fa")

ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy))#,
              # phy_tree(tree),
              # refseq(rep)
)


##################################################################################

# load data
load("../2.DX2013_145sites/data/resample.envi1314.comm1314.taxa1314.Rdata")

load("../2.DX2013_145sites/data/1.comm.envi.div.all.Rdata")

cli1314.last = read.csv("../2.DX2013_145sites/data/predict.act.ele1314.bio12.bio1.AI.GT.last.csv",
                        row.names = 1)

plant.RA <- read.csv("../2.DX2013_145sites/data/last.plant1314.188sites.species.number.for.saiz.method.csv",
                     row.names = 1)

plant.RA[is.na(plant.RA)] <- 0


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

taxa.comm.mv = taxa1314.samp[colnames(comm.mv), ]

write.csv(taxa.comm.mv, "../2.DX2013_145sites/data/Bipartite20230609/4400/taxa.comm.mv.csv")

### DX data ###
envi.4400 = read.csv("../2.DX2013_145sites/data/Bipartite20230609/4400/envi.4400.csv", row.names = 1)
comm.4400 = read.csv("../2.DX2013_145sites/data/Bipartite20230609/4400/otu4400.csv", row.names = 1)
plant.4400 = read.csv("../2.DX2013_145sites/data/Bipartite20230609/4400/plant.RA4400.csv", row.names = 1)
taxa.comm.mv.la = read.csv("../2.DX2013_145sites/data/Bipartite20230609/4400/taxa.comm.mv.csv", row.names = 1)
taxa.plant.la = read.csv("../2.DX2013_145sites/data/Bipartite20230609/4400/taxa.plant.csv", row.names = 1)



ps.comm4400 = phyloseq(sample_data(envi.4400),
              otu_table(as.matrix(comm.4400), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxa.comm.mv.la))#,
              # phy_tree(tree),
              # refseq(rep)
)

ps.plant4400 = phyloseq(sample_data(envi.4400),
                       otu_table(as.matrix(plant.4400), taxa_are_rows=TRUE),
                       tax_table(as.matrix(taxa.plant.la))#,
                       # phy_tree(tree),
                       # refseq(rep)
)

out.plant = rbind(comm.4400, plant.4400)
out.plant.taxa = rbind(taxa.comm.mv.la, taxa.plant.la)


ps.plant4400 = phyloseq(sample_data(envi.4400),
                        otu_table(as.matrix(plant.4400), taxa_are_rows=TRUE),
                        tax_table(as.matrix(taxa.plant.la))#,
                        # phy_tree(tree),
                        # refseq(rep)
)


#--细菌和真菌ps对象中的map文件要一样

ps.mergeDX <- ggClusterNet::merge16S_ITS(ps16s = ps.comm4400,
                                         psITS = ps.plant4400)

ps.mergeDX

mapDX =  phyloseq::sample_data(ps.mergeDX)

head(mapDX)
mapDX$Group = "4400"
phyloseq::sample_data(ps.mergeDX) <- mapDX

ps.mergeDX.otu_table = data.frame(phyloseq::otu_table(ps.mergeDX))

DXresult <- corBionetwork(ps = ps.mergeDX,
                        N = 0,
                        lab = data,
                        r.threshold = 0.2, # 相关阈值
                        p.threshold = 0.05,
                        group = "Group",
                        # env = data1, # 环境指标表格
                        # envGroup = Gru,# 环境因子分组文件表格
                        # layout = "fruchtermanreingold",
                        path = Envnetplot,# 结果文件存储路径
                        fill = "Phylum", # 出图点填充颜色用什么值
                        size = "igraph.degree", # 出图点大小用什么数据
                        scale = TRUE, # 是否要进行相对丰度标准化
                        bio = TRUE, # 是否做二分网络
                        zipi = F, # 是否计算ZIPI
                        step = 100, # 随机网络抽样的次数
                        width = 12,
                        label = TRUE,
                        height = 10,
                        big = TRUE,
                        select_layout = TRUE,
                        layout_net = "model_maptree2",
                        clu_method = "cluster_fast_greedy"
                        
                        
)

tem <- model_maptree(cor =DXresult[[5]],
                     method = "cluster_fast_greedy",
                     seed = 12
)
node_model = tem[[2]]
head(node_model)

p = DXresult[[1]]
p
# 全部样本网络参数比对
data = DXresult[[2]]
