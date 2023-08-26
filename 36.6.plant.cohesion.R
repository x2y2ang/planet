
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


# # use 6 for plant abundance data and as majority selection criterion
# veg.tmp = t(veg1314)
# veg.mv <- veg.tmp[rowSums(veg.tmp)>= 6,]
# 
# # occurring in 30% of all samples 188*0.3 = 56
# dim(veg.tmp)
# veg.R1 = veg.tmp
# veg.R1[veg.R1 >0] <- 1
# veg.mv <- veg.tmp[which(rowSums(veg.R1) >= 56), ]  
# dim(veg.mv)


# #####  combine plant micro ######
# veg.mv = t(veg.mv)
# micro188 = comm.mv[rownames(veg.mv), ]
# 
# veg.micro = cbind(veg.mv, micro188)

####  for micro #####
otu = veg1314

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

b <- veg1314
occor <- rcorr(as.matrix(as.data.frame(lapply(veg1314, as.numeric))), 
      type = "spearman")
#occor<-rcorr(b, type = "spearman")
occor.r = occor$r 
occor.p = occor$P 
# Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
occor.p<-p.adjust(occor.p, method="BH")
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0 
diag(occor.r) <- 0   

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
          # ?滻????ֵ
          perm.rel.d[, j ] <- sample(rel.d[ ,j ]) 
        }
        # focal????ά??ԭ״?????? 
        perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]
        
        # ????????
        cor.mat.null <- cor(perm.rel.d)
        
        # ????ÿ?ζ???focal???ֵ????Ľ???
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      
      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  } else { ##?ڶ?????ģ?ͣ??????????????з??ȣ???????
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
        
        ###????????
        cor.mat.null <- cor(perm.rel.d)
        ###????ÿ?ζ???focal???ֵ????Ľ???
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  }
}

####ʵ??????-??ģ??????
ifelse(use.custom.cors == T, {
  obs.exp.cors.mat <- custom.cor.mat.sub}, {
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }
)

diag(obs.exp.cors.mat) <- 0


####???????????ص?��ͨ??
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)
###???ݶ??壬????cohesion
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

##???? 
output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("Negative Connectedness", "Positive Connectedness", 
                   "Negative Cohesion", "Positive Cohesion")

print(output)


write.csv(output$`Negative Cohesion`, "../2.DX2013_145sites/data/plant1314.Negative.cohesion.csv")
write.csv(output$`Positive Cohesion`, "../2.DX2013_145sites/data/plant1314.Positive.cohesion.csv")


#### plot ####
tmp = data.frame(output$`Negative Cohesion`)
cli185 = cli1314.last[rownames(tmp), ]

AI.coh = cbind(cli185, output$`Negative Cohesion`, output$`Positive Cohesion`)
colnames(AI.coh)[11:12] = c("Negative.Cohesion", "Positive.Cohesion")


plot(AI.coh$Elevation, abs(AI.coh$Negative.Cohesion))
plot(AI.coh$Aridity.index, abs(AI.coh$Negative.Cohesion))
plot(AI.coh$pre.AI, abs(AI.coh$Negative.Cohesion))

plot(AI.coh$Elevation, abs(AI.coh$Positive.Cohesion))
plot(AI.coh$Aridity.index, abs(AI.coh$Positive.Cohesion))

write.csv(AI.coh, "../2.DX2013_145sites/data/plant1314.cohesion.combine.with.AI.csv")


library(vegan)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


AI.coh$Elevation = factor(AI.coh$Elevation)

colnames(AI.coh)

p1 = ggplot(AI.coh, aes(pre.AI, abs(Negative.Cohesion), color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y ~ x + I(x^2), color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+ 
 # facet_wrap( ~ Group, scales = "free", ncol = 3)+ 
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
  xlab("Aridity")+ylab("Plant negative cohesion")+
  theme(legend.position = 'none')
  #ylim(0.29, 0.4)

p1


p2 = ggplot(AI.coh, aes(pre.AI, Positive.Cohesion, color = Elevation))+
  geom_point(size = 2)+
  stat_smooth(method="lm", formula = y ~ x + I(x^2), color = "black")+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+ 
  # facet_wrap( ~ Group, scales = "free", ncol = 3)+ 
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
  xlab("Aridity")+ylab("Plant positive cohesion")+
  theme(legend.position = 'none')
#ylim(0.29, 0.4)

p2

library(cowplot)
plot_grid(p1,p2,ncol = 2)
ggsave("../2.DX2013_145sites/figures/Elevation.aridity.envi.plant.cohesion.1314year.pdf",
       width = 6.5, height = 3)
