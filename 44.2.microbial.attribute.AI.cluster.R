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

AI = cli1314.last[, c("pre.AI")]

AI.comm = cor(comm.mv, AI)

AI.test = corr.test(comm.mv, AI, use="pairwise", method="spearman",
                    adjust="fdr", alpha=.05)

occor.r = AI.test$r 
occor.p = AI.test$p 
occor.r[occor.p>0.05|abs(occor.r)<0.5] = 0 

write.csv(occor.r, "../2.DX2013_145sites/data/60.AI.otu.comm.cor.test.p.r.csv")


GTmin = cli1314.last[, c("GTmin")]

GTmin.comm = cor(comm.mv, GTmin)

GTmin.test = corr.test(comm.mv, GTmin, use="pairwise", method="spearman",
                    adjust="fdr", alpha=.05)

occor.r1 = GTmin.test$r 
occor.p1 = GTmin.test$p 
occor.r1[occor.p1>0.05|abs(occor.r1)<0.5] = 0 

write.csv(occor.r1, "../2.DX2013_145sites/data/61.GTmin.otu.comm.cor.test.p.r.csv")
