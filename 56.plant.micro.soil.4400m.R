
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")


library(Hmisc)
library(reshape2)
library(igraph)


####读取三组数据
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


veg.mv1 = t(veg.mv)
micro188 = comm.mv[rownames(veg.mv1), ]
envi188 = envi1314[rownames(veg.mv1),]
envi188.select = envi188[, c(11:19)]


# ele = 4400
ele44 = subset(envi188, Elevation == "4400")

veg44 = veg.mv1[rownames(ele44), ]
micro44 = micro188[rownames(ele44), ]
envi44 = envi188.select[rownames(ele44), ]


#### micro-plant spearman ####
corr44 <- rcorr(as.matrix(micro44), as.matrix(veg44), 
                         type = 'spearman')

# Matrix of correlation coefficient r-values and significance p-values
MP.r <- corr44$r
MP.p <- corr44$P

# Only the correlation coefficient of microbial abundance-plant abundance is retained
# Remove the correlation coefficient between microbes-microbes, plant-plant
MP.r1 <- MP.r[colnames(micro44),colnames(veg44)]
MP.p1 <- MP.p[colnames(micro44),colnames(veg44)]

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
MP.p1 <- p.adjust(MP.p1, method = 'BH')   

MP.r1[MP.p1>0.05|abs(MP.r1)<0.2] = 0 


#### plant and envi spearman correlation
pecorr44 <- rcorr(as.matrix(envi44), as.matrix(veg44), type = 'spearman')

r <- pecorr44$r
p <- pecorr44$P

r <- r[colnames(envi44),colnames(veg44)]
p <- p[colnames(envi44),colnames(veg44)]

p <- p.adjust(p, method = 'BH')  

r[p>0.05|abs(r)<0.15] = 0 


# ele = 4300
ele43 = subset(envi188, Elevation == "4300")

veg43 = veg.mv1[rownames(ele43), ]
micro43 = micro188[rownames(ele43), ]
envi43 = envi188.select[rownames(ele43), ]


mvcorr43 <- rcorr(as.matrix(micro43), as.matrix(veg43), 
                         type = 'spearman')

MP.r <- mvcorr43$r
MP.p <- mvcorr43$P

MP.r1 <- MP.r[colnames(micro43),colnames(veg43)]
MP.p1 <- MP.p[colnames(micro43),colnames(veg43)]
MP.p1 <- p.adjust(MP.p1, method = 'BH')   

# Choose a correlation coefficient with a significant p-value less than 0.05, i.e. p<0.05
MP.r1[MP.p1>0.05|abs(MP.r1)<0.2] = 0 


#### plant and envi spearman correlation
envi188.select_veg.shift_corr <- rcorr(as.matrix(envi43), as.matrix(veg43), type = 'spearman')

r <- envi188.select_veg.shift_corr$r
p <- envi188.select_veg.shift_corr$P

r <- r[colnames(envi43),colnames(veg43)]
p <- p[colnames(envi43),colnames(veg43)]

p <- p.adjust(p, method = 'BH')  

r[p>0.05|abs(r)<0.15] = 0 

# ele = 4600
ele46 = subset(envi188, Elevation == "4600")

veg46 = veg.mv1[rownames(ele46), ]
micro46 = comm.mv[rownames(ele46), ]
envi46 = envi188.select[rownames(ele46), ]


#### micro-plant spearman ####
genus_args_corr <- rcorr(as.matrix(micro46), as.matrix(veg46), 
                         type = 'spearman')

# Matrix of correlation coefficient r-values and significance p-values
MP.r <- genus_args_corr$r
MP.p <- genus_args_corr$P

# Only the correlation coefficient of microbial abundance-plant abundance is retained
# Remove the correlation coefficient between microbes-microbes, plant-plant
MP.r1 <- MP.r[colnames(micro188),colnames(veg46)]
MP.p1 <- MP.p[colnames(micro188),colnames(veg46)]

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
MP.p1 <- p.adjust(MP.p1, method = 'BH')   

MP.r1[MP.p1>0.05|abs(MP.r1)<0.2] = 0 


#### plant and envi spearman correlation
envi188.select_veg.shift_corr <- rcorr(as.matrix(envi46), as.matrix(veg46), type = 'spearman')

r <- envi188.select_veg.shift_corr$r
p <- envi188.select_veg.shift_corr$P

r <- r[colnames(envi46),colnames(veg46)]
p <- p[colnames(envi46),colnames(veg46)]

p <- p.adjust(p, method = 'BH')  

r[p>0.05|abs(r)<0.15] = 0 



# ele = 5000
ele50 = subset(envi188, Elevation == "5000")

veg50 = veg.mv1[rownames(ele50), ]
micro50 = comm.mv[rownames(ele50), ]
envi50 = envi188.select[rownames(ele50), ]


#### micro-plant spearman ####
genus_args_corr <- rcorr(as.matrix(micro50), as.matrix(veg50), 
                         type = 'spearman')

# Matrix of correlation coefficient r-values and significance p-values
MP.r <- genus_args_corr$r
MP.p <- genus_args_corr$P

# Only the correlation coefficient of microbial abundance-plant abundance is retained
# Remove the correlation coefficient between microbes-microbes, plant-plant
MP.r1 <- MP.r[colnames(micro188),colnames(veg50)]
MP.p1 <- MP.p[colnames(micro188),colnames(veg50)]

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
MP.p1 <- p.adjust(MP.p1, method = 'BH')   

MP.r1[MP.p1>0.05|abs(MP.r1)<0.2] = 0 


#### plant and envi spearman correlation
envi188.select_veg.shift_corr <- rcorr(as.matrix(envi50), as.matrix(veg50), type = 'spearman')

r <- envi188.select_veg.shift_corr$r
p <- envi188.select_veg.shift_corr$P

r <- r[colnames(envi50),colnames(veg50)]
p <- p[colnames(envi50),colnames(veg50)]

p <- p.adjust(p, method = 'BH')  

r[p>0.05|abs(r)<0.15] = 0 

