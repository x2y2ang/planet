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

write.csv(cli.pre.AI, "../2.DX2013_145sites/data/cli.pre.AI188sites.csv")

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

### combine together ###
PMS = cbind(plant.cov, envi.tmp[, c(11:19)], comm.tmp)

PMS[1:8, 1:8]

pers.cutoff <- 0.20
iter <- 200
tax.shuffle <- T
use.custom.cors <- F

zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}

neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}

pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}

b <- PMS


####
if(use.custom.cors == T) {
  custom.cor.mat <- occor.r
  custom.cor.mat <- as.matrix(custom.cor.mat)
  # Check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2] == dim(custom.cor.mat)[2])
}

####
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]

###
rowsums.orig <- rowSums(c)

###
zero.cutoff <- ceiling(pers.cutoff * dim(c)[1])

###
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]

###
d <- d[rowSums(d) > 0, ]


if(use.custom.cors == T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c, 2, zero) < (dim(c)[1]-zero.cutoff), 
                                       apply(c, 2, zero) < (dim(c)[1]-zero.cutoff)]
}

rel.d <- d / rowsums.orig
hist(rowSums(rel.d))

cor.mat.true <- cor(rel.d)

med.tax.cors <- vector()
if(use.custom.cors == F) {
  if(tax.shuffle) { 
    for(which.taxon in 1:dim(rel.d)[2]){
      
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        #Create empty matrix of same dimension as rel.d
        perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #For each otu
        for(j in 1:dim(rel.d)[2]){ 
         
          perm.rel.d[, j ] <- sample(rel.d[ ,j ]) 
        }
        
        perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]
        
        cor.mat.null <- cor(perm.rel.d)
        
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      
      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  } else { 
    for(which.taxon in 1:dim(rel.d)[2]){
      
      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        #Create duplicate matrix to shuffle abundances
        perm.rel.d <- rel.d 
        
        #For each taxon
        for(j in 1:dim(rel.d)[1]){ 
          which.replace <- which(rel.d[j, ] > 0 ) 
          # if the focal taxon is greater than zero, take it out of the replacement vector, so the focal abundance stays the same
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          
          #Replace the original taxon vector with a vector where the values greater than 0 have been randomly permuted 
          perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal]) 
        }
        
        cor.mat.null <- cor(perm.rel.d)
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  }
}


ifelse(use.custom.cors == T, {
  obs.exp.cors.mat <- custom.cor.mat.sub}, {
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }
)

diag(obs.exp.cors.mat) <- 0



connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("Negative Connectedness", "Positive Connectedness", 
                   "Negative Cohesion", "Positive Cohesion")

print(output)

write.csv(output$`Negative Cohesion`, "../2.DX2013_145sites/data/plant188.micro.soil.1314.Negative.cohesion.csv")
write.csv(output$`Positive Cohesion`, "../2.DX2013_145sites/data/plant188.micro.soil.1314.Positive.cohesion.csv")
