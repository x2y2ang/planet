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


comm.tmp = comm.mv[rownames(plant.cov), ]


envi.tmp = envi1314[rownames(plant.cov), ]


### < 8 samples
envi.4300 = subset(envi.tmp, Elevation == 4300)
envi.4300 = envi.4300[, 11:19]

plant4300 = plant.cov[rownames(envi.4300), ]

otu4300 = comm.tmp[rownames(envi.4300), ]

plant4300 = t(plant4300)
otu4300 = t(otu4300)


plant.RA4300 = plant.RA[rownames(envi.4300), ]
plant.RA4300 = t(plant.RA4300)

write.csv(plant4300, "../2.DX2013_145sites/data/plant.cov.4300.csv")
write.csv(otu4300, "../2.DX2013_145sites/data/otu4300.csv")
write.csv(envi.4300, "../2.DX2013_145sites/data/envi.4300.csv")
write.csv(plant.RA4300, "../2.DX2013_145sites/data/plant.RA4300.csv")


### 4400 m ####
envi.4400 = subset(envi.tmp, Elevation == 4400)
envi.4400 = envi.4400[, 11:19]

plant4400 = plant.cov[rownames(envi.4400), ]

otu4400 = comm.tmp[rownames(envi.4400), ]

plant4400 = t(plant4400)
otu4400 = t(otu4400)


plant.RA4400 = plant.RA[rownames(envi.4400), ]
plant.RA4400 = t(plant.RA4400)

write.csv(plant4400, "../2.DX2013_145sites/data/plant.cov.4400.csv")
write.csv(otu4400, "../2.DX2013_145sites/data/otu4400.csv")
write.csv(envi.4400, "../2.DX2013_145sites/data/envi.4400.csv")
write.csv(plant.RA4400, "../2.DX2013_145sites/data/plant.RA4400.csv")


### 4500 m ####
envi.4500 = subset(envi.tmp, Elevation == 4500)
envi.4500 = envi.4500[, 11:19]

plant4500 = plant.cov[rownames(envi.4500), ]

otu4500 = comm.tmp[rownames(envi.4500), ]

plant4500 = t(plant4500)
otu4500 = t(otu4500)

plant.RA4500 = plant.RA[rownames(envi.4500), ]
plant.RA4500 = t(plant.RA4500)

write.csv(plant4500, "../2.DX2013_145sites/data/plant.cov.4500.csv")
write.csv(otu4500, "../2.DX2013_145sites/data/otu4500.csv")
write.csv(envi.4500, "../2.DX2013_145sites/data/envi.4500.csv")
write.csv(plant.RA4500, "../2.DX2013_145sites/data/plant.RA4500.csv")




### 4600 m ####
envi.4600 = subset(envi.tmp, Elevation == 4600)
envi.4600 = envi.4600[, 11:19]

plant4600 = plant.cov[rownames(envi.4600), ]

otu4600 = comm.tmp[rownames(envi.4600), ]

plant4600 = t(plant4600)
otu4600 = t(otu4600)

plant.RA4600 = plant.RA[rownames(envi.4600), ]
plant.RA4600 = t(plant.RA4600)

write.csv(plant4600, "../2.DX2013_145sites/data/plant.cov.4600.csv")
write.csv(otu4600, "../2.DX2013_145sites/data/otu4600.csv")
write.csv(envi.4600, "../2.DX2013_145sites/data/envi.4600.csv")
write.csv(plant.RA4600, "../2.DX2013_145sites/data/plant.RA4600.csv")




### 4700 m ####
envi.4700 = subset(envi.tmp, Elevation == 4700)
envi.4700 = envi.4700[, 11:19]

plant4700 = plant.cov[rownames(envi.4700), ]

otu4700 = comm.tmp[rownames(envi.4700), ]

plant4700 = t(plant4700)
otu4700 = t(otu4700)

plant.RA4700 = plant.RA[rownames(envi.4700), ]
plant.RA4700 = t(plant.RA4700)

write.csv(plant4700, "../2.DX2013_145sites/data/plant.cov.4700.csv")
write.csv(otu4700, "../2.DX2013_145sites/data/otu4700.csv")
write.csv(envi.4700, "../2.DX2013_145sites/data/envi.4700.csv")
write.csv(plant.RA4700, "../2.DX2013_145sites/data/plant.RA4700.csv")




### 4800 m ####
envi.4800 = subset(envi.tmp, Elevation == 4800)
envi.4800 = envi.4800[, 11:19]

plant4800 = plant.cov[rownames(envi.4800), ]

otu4800 = comm.tmp[rownames(envi.4800), ]

plant4800 = t(plant4800)
otu4800 = t(otu4800)

plant.RA4800 = plant.RA[rownames(envi.4800), ]
plant.RA4800 = t(plant.RA4800)

write.csv(plant4800, "../2.DX2013_145sites/data/plant.cov.4800.csv")
write.csv(otu4800, "../2.DX2013_145sites/data/otu4800.csv")
write.csv(envi.4800, "../2.DX2013_145sites/data/envi.4800.csv")
write.csv(plant.RA4800, "../2.DX2013_145sites/data/plant.RA4800.csv")



### 4900 m ####
envi.4900 = subset(envi.tmp, Elevation == 4900)
envi.4900 = envi.4900[, 11:19]

plant4900 = plant.cov[rownames(envi.4900), ]

otu4900 = comm.tmp[rownames(envi.4900), ]

plant4900 = t(plant4900)
otu4900 = t(otu4900)

plant.RA4900 = plant.RA[rownames(envi.4900), ]
plant.RA4900 = t(plant.RA4900)

write.csv(plant4900, "../2.DX2013_145sites/data/plant.cov.4900.csv")
write.csv(otu4900, "../2.DX2013_145sites/data/otu4900.csv")
write.csv(envi.4900, "../2.DX2013_145sites/data/envi.4900.csv")
write.csv(plant.RA4900, "../2.DX2013_145sites/data/plant.RA4900.csv")



### 5000 m ####
envi.5000 = subset(envi.tmp, Elevation == 5000)
envi.5000 = envi.5000[, 11:19]

plant5000 = plant.cov[rownames(envi.5000), ]

otu5000 = comm.tmp[rownames(envi.5000), ]

plant5000 = t(plant5000)
otu5000 = t(otu5000)

plant.RA5000 = plant.RA[rownames(envi.5000), ]
plant.RA5000 = t(plant.RA5000)

write.csv(plant5000, "../2.DX2013_145sites/data/plant.cov.5000.csv")
write.csv(otu5000, "../2.DX2013_145sites/data/otu5000.csv")
write.csv(envi.5000, "../2.DX2013_145sites/data/envi.5000.csv")
write.csv(plant.RA5000, "../2.DX2013_145sites/data/plant.RA5000.csv")



### 5100 m ####
envi.5100 = subset(envi.tmp, Elevation == 5100)
envi.5100 = envi.5100[, 11:19]

plant5100 = plant.cov[rownames(envi.5100), ]

otu5100 = comm.tmp[rownames(envi.5100), ]

plant5100 = t(plant5100)
otu5100 = t(otu5100)

plant.RA5100 = plant.RA[rownames(envi.5100), ]
plant.RA5100 = t(plant.RA5100)

write.csv(plant5100, "../2.DX2013_145sites/data/plant.cov.5100.csv")
write.csv(otu5100, "../2.DX2013_145sites/data/otu5100.csv")
write.csv(envi.5100, "../2.DX2013_145sites/data/envi.5100.csv")
write.csv(plant.RA5100, "../2.DX2013_145sites/data/plant.RA5100.csv")


### 5200 m ####
envi.5200 = subset(envi.tmp, Elevation == 5200)
envi.5200 = envi.5200[, 11:19]

plant5200 = plant.cov[rownames(envi.5200), ]

otu5200 = comm.tmp[rownames(envi.5200), ]

plant5200 = t(plant5200)
otu5200 = t(otu5200)

plant.RA5200 = plant.RA[rownames(envi.5200), ]
plant.RA5200 = t(plant.RA5200)

write.csv(plant5200, "../2.DX2013_145sites/data/plant.cov.5200.csv")
write.csv(otu5200, "../2.DX2013_145sites/data/otu5200.csv")
write.csv(envi.5200, "../2.DX2013_145sites/data/envi.5200.csv")
write.csv(plant.RA5200, "../2.DX2013_145sites/data/plant.RA5200.csv")




### combine together ###
PMS = cbind(plant.cov,  comm.tmp) # envi.tmp[, c(11:19)],


?graph.data.frame
g <- graph.data.frame(PMS, directed=FALSE)

bipartite.mapping(g)

library(bipartite)
# NOT RUN {
data(Safariland)

a = networklevel(Safariland)

b = networklevel(Safariland, index="ALLBUTDD") #excludes degree distribution fits
# }



## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))

relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))

g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)
print(g, e=TRUE, v=TRUE)

## The opposite operation
as_data_frame(g, what="vertices")
as_data_frame(g, what="edges")

bipartite.mapping(g)

V(g)$type <- bipartite_mapping(g)$type 
plot(g)



# Clear environment
rm(list=ls())

# Load required packages for bipartite
library(permute)
library(lattice)
library(vegan)
library(statnet.common)
library(network)
library(sna)
library(bipartite)

# Read all the interaction network matrices into a list
webs <- list(Safariland, barrett1987, elberling1999, 
             memmott1999, motten1982, olesen2002aigrettes)

# Re-name the datasets according to the sites for each plant-pollinator network
webs.names <- c("Argentina", "Canada", "Sweden", "UK", "USA", "Azores") 
names(webs) <- webs.names


# View data (interaction matrix, for generating a network, or a web)
lapply(webs, head, n = 2L) # Only display the first two rows in the dataset

visweb(webs$Argentina)

plotweb(webs$Argentina, text.rot=90, col.low = "green", col.high = "blue")


# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest <- lapply(webs, networklevel, index = 'nestedness') 

# Calculate network metric links per species for all plant-pollinator sites
net.metrics.links <- lapply(webs, networklevel, index = 'links per species') 


# Time consuming step!

# Load environment (already saved objects)
#load("data/network_analysis_example.RData")

# Make null models for all sites using the r2dtable null
net.nulls.r2d <- lapply(webs, nullmodel, method = "r2dtable", N = 500) 

# Make null models for all sites using the vaznull null
net.nulls.vaz <- lapply(webs, nullmodel, method = "vaznull", N = 500) 

# Make null models for all sites using the swap.web null
net.nulls.swap <- lapply(webs, nullmodel, method = "swap.web", N = 500)

# Save null objects
#save(net.nulls.r2d, net.nulls.vaz, net.nulls.swap, file = "data/network_analysis_example.RData")