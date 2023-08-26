
pv.fun1 <- function(var1){
  if(var1[1]==var1[2]){return(0)}else{
    return((1-(min(var1[1],var1[2])/max(var1[1],var1[2]))))}
}
pv.fun <- function(var, k=0){
  var<-var[!is.na(var)]
  if(length(var)<3){return(NA)}
  if (min(var)<0){var <- var + abs(min(var)) + 0.01*(max(var)-min(var))}
  var <- var + k
  if(length(which(is.na(var)))>=(length(var)-1)){return(NA)}else{
    var<-as.numeric(na.omit(var))
    if (any(var<0)){var <- var + abs(min(var)) + 1}
    return(mean(as.numeric(combn(var, 2, FUN = pv.fun1, simplify = F))))}}


# 设置随机数种子，保证结果可复现
set.seed(123)

# 数据集，这里假设14个数据点存在一个向量中
data <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)

# 定义重复实验次数
num_repeats <- 100

# 初始化一个向量来存储每次实验的PV值
PV_values <- numeric(num_repeats)

# 开始蒙特卡洛模拟循环
for (i in 1:num_repeats) {
  # 随机抽取9个点
  sampled_data <- sample(data, size = 9, replace = FALSE)
  
  # 计算PV
  PV_values[i] <- pv.fun(sampled_data)
}

# 输出每次实验的PV值
print(PV_values)

# 输出PV的均值和标准差
print(paste("Mean PV:", mean(PV_values)))
print(paste("SD PV:", sd(PV_values)))



# 设置随机数种子，保证结果可复现
set.seed(123)

# 十个海拔的数据集
altitude_data <- list(
  altitude1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  altitude2 = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
  altitude3 = c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18),
  altitude4 = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23),
  altitude5 = c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28),
  altitude6 = c(20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33),
  altitude7 = c(25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38),
  altitude8 = c(30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43),
  altitude9 = c(35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48),
  altitude10 = c(40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53)
)

# 定义重复实验次数
num_repeats <- 1000

# 初始化一个列表，用于存储每个海拔的CV值
cv_values_list <- list()
asy_values_list <- list()
pv_values_list <- list()

# 开始双重循环
for (altitude_name in names(altitude_data)) {
  # 获取当前海拔的数据集
  data <- altitude_data[[altitude_name]]
  
  # 初始化一个向量来存储每次实验的CV值
  cv_values <- numeric(num_repeats)
  asy_values <- numeric(num_repeats)
  pv_values <- numeric(num_repeats)
  
  # 计算抽取数据的平均值和标准差
  Dmean <- mean(sampled_data)
  Dsd <- sd(sampled_data)
  
  # Asy functions
  Dmax = max(sampled_data)
  Dmin = min(sampled_data)

  # 内层循环进行蒙特卡洛模拟
  for (i in 1:num_repeats) {
    # 随机抽取9个点
    sampled_data <- sample(data, size = 9, replace = FALSE)
    # 计算CV
    cv_values[i] <- (Dsd / Dmean)
    asy_values[i] = ((Dmax - Dmean)/(Dmean - Dmin))
    pv_values[i] <- pv.fun(sampled_data)
    
  }
  
  # 将当前海拔的CV值存储到列表中
  cv_values_list[[altitude_name]] <- cv_values
  asy_values_list[[altitude_name]] <- asy_values
  pv_values_list[[altitude_name]] <- pv_values
}




### 自己的数据 试试看 ###



rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(vegan)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

veg188.cli.tmp = read.csv("../2.DX2013_145sites/data/DX1314.Plant.rich.AI.for.PV.csv",
                          row.names = 1)

Plant.micro.soil.EMF = read.csv("../2.DX2013_145sites/data/DX1314.Plant.micro.soil.EMF.AI.for.PV.csv",
                                row.names = 1)


### Plant.rich #### 4300 只有4个点，需要单独算
veg188.de4300 = subset(veg188.cli.tmp, Elevation != 4300)


# 当雄数据最多29个点抽取9个点的组合数
n <- 29
k <- 9
combinations <- choose(n, k)

# 输出结果
print(combinations)

# 当雄数据每个海拔15个点抽取9个点的组合数
n <- 15
k <- 9
combinations <- choose(n, k)

# 输出结果
print(combinations)
# 5005

# 定义重复实验次数
num_repeats <- 5000

plant.rich.tmp = data.frame(veg188.de4300[, "Plant.rich"])
ele = unique(veg188.de4300[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  data <- plant.rich.tmp[which(veg188.de4300$Elevation == ele.tmp),]
  

  # 内层循环进行蒙特卡洛模拟
  for (j in 1:num_repeats) {
    # 随机抽取9个点
    sampled_data <- sample(data, size = 9, replace = FALSE)
    # 计算CV
    # 计算抽取数据的平均值和标准差
    Dmean <- mean(sampled_data)
    Dsd <- sd(sampled_data)
    
    # Asy functions
    Dmax <- max(sampled_data)
    Dmin <- min(sampled_data)
    
    cv_values <- (Dsd / Dmean)
    asy_values <- ((Dmax - Dmean)/(Dmean - Dmin))
    pv_values<- pv.fun(sampled_data)
    ele.cv.asy.pv.tmp <- data.frame(cbind(ele.tmp, cv_values, asy_values, pv_values))
    
    if (j==1){
      ele.cv.asy.pv = ele.cv.asy.pv.tmp
    }
    else{
      ele.cv.asy.pv = rbind(ele.cv.asy.pv, ele.cv.asy.pv.tmp)
    }
    }


if (i==1)
{
  plant.rich.cv.asy.pv = ele.cv.asy.pv
}
else{
  plant.rich.cv.asy.pv = rbind(plant.rich.cv.asy.pv, ele.cv.asy.pv)
}
}

write.csv(plant.rich.cv.asy.pv, "../2.DX2013_145sites/data/69.1.plant.rich.cv.asy.pv.rep5000times.csv")




### Plant.rich43 #### 单独算 4300(只有4个点)
Plant.rich43 = subset(veg188.cli.tmp, Elevation == 4300)


# 定义重复实验次数
num_repeats <- 4

plant.rich.tmp = data.frame(Plant.rich43[, "Plant.rich"])
ele = unique(Plant.rich43[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  data <- plant.rich.tmp[which(Plant.rich43$Elevation == ele.tmp),]
  
  
  # 内层循环进行蒙特卡洛模拟
  for (j in 1:num_repeats) {
    # 随机抽取3个点
    sampled_data <- sample(data, size = 3, replace = FALSE)
    # 计算CV
    # 计算抽取数据的平均值和标准差
    Dmean <- mean(sampled_data)
    Dsd <- sd(sampled_data)
    
    # Asy functions
    Dmax <- max(sampled_data)
    Dmin <- min(sampled_data)
    
    cv_values <- (Dsd / Dmean)
    asy_values <- ((Dmax - Dmean)/(Dmean - Dmin))
    pv_values<- pv.fun(sampled_data)
    ele.cv.asy.pv.tmp <- data.frame(cbind(cv_values, asy_values, pv_values))
    
    if (j==1){
      ele.cv.asy.pv = ele.cv.asy.pv.tmp
    }
    else{
      ele.cv.asy.pv = rbind(ele.cv.asy.pv, ele.cv.asy.pv.tmp)
    }
  }
  
  
  if (i==1)
  {
    plant.rich43.cv.asy.pv = ele.cv.asy.pv
  }
  else{
    plant.rich43.cv.asy.pv = rbind(plant.rich43.cv.asy.pv, ele.cv.asy.pv)
  }
}

write.csv(plant.rich43.cv.asy.pv, "../2.DX2013_145sites/data/69.1.plant.rich4300m.cv.asy.pv.rep4times.csv")



### Micro.rich ####

# 定义重复实验次数
num_repeats <- 5000

Micro.rich.tmp = data.frame(Plant.micro.soil.EMF[, "Micro.rich"])
ele = unique(Plant.micro.soil.EMF[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  data <- Micro.rich.tmp[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  # 内层循环进行蒙特卡洛模拟
  for (j in 1:num_repeats) {
    # 随机抽取9个点
    sampled_data <- sample(data, size = 9, replace = FALSE)
    # 计算CV
    # 计算抽取数据的平均值和标准差
    Dmean <- mean(sampled_data)
    Dsd <- sd(sampled_data)
    
    # Asy functions
    Dmax <- max(sampled_data)
    Dmin <- min(sampled_data)
    
    cv_values <- (Dsd / Dmean)
    asy_values <- ((Dmax - Dmean)/(Dmean - Dmin))
    pv_values<- pv.fun(sampled_data)
    ele.cv.asy.pv.tmp <- data.frame(cbind(ele.tmp, cv_values, asy_values, pv_values))
    
    if (j==1){
      ele.cv.asy.pv = ele.cv.asy.pv.tmp
    }
    else{
      ele.cv.asy.pv = rbind(ele.cv.asy.pv, ele.cv.asy.pv.tmp)
    }
  }
  
  
  if (i==1)
  {
    Micro.rich.cv.asy.pv = ele.cv.asy.pv
  }
  else{
    Micro.rich.cv.asy.pv = rbind(Micro.rich.cv.asy.pv, ele.cv.asy.pv)
  }
}

write.csv(Micro.rich.cv.asy.pv, "../2.DX2013_145sites/data/69.2.Micro.rich.cv.asy.pv.rep5000times.csv")



### EMF ####

# 定义重复实验次数
num_repeats <- 5000

EMF.tmp = data.frame(Plant.micro.soil.EMF[, "EMF"])
ele = unique(Plant.micro.soil.EMF[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  data <- EMF.tmp[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  # 内层循环进行蒙特卡洛模拟
  for (j in 1:num_repeats) {
    # 随机抽取9个点
    sampled_data <- sample(data, size = 9, replace = FALSE)
    # 计算CV
    # 计算抽取数据的平均值和标准差
    Dmean <- mean(sampled_data)
    Dsd <- sd(sampled_data)
    
    # Asy functions
    Dmax <- max(sampled_data)
    Dmin <- min(sampled_data)
    
    cv_values <- (Dsd / Dmean)
    asy_values <- ((Dmax - Dmean)/(Dmean - Dmin))
    pv_values<- pv.fun(sampled_data)
    ele.cv.asy.pv.tmp <- data.frame(cbind(ele.tmp, cv_values, asy_values, pv_values))
    
    if (j==1){
      ele.cv.asy.pv = ele.cv.asy.pv.tmp
    }
    else{
      ele.cv.asy.pv = rbind(ele.cv.asy.pv, ele.cv.asy.pv.tmp)
    }
  }
  
  
  if (i==1)
  {
    EMF.cv.asy.pv = ele.cv.asy.pv
  }
  else{
    EMF.cv.asy.pv = rbind(EMF.cv.asy.pv, ele.cv.asy.pv)
  }
}

write.csv(EMF.cv.asy.pv, "../2.DX2013_145sites/data/69.3.EMF.cv.asy.pv.rep5000times.csv")



### AGB ####

# 定义重复实验次数
num_repeats <- 5000

AGB.tmp = data.frame(Plant.micro.soil.EMF[, "AGB"])
ele = unique(Plant.micro.soil.EMF[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  data <- AGB.tmp[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  # 内层循环进行蒙特卡洛模拟
  for (j in 1:num_repeats) {
    # 随机抽取9个点
    sampled_data <- sample(data, size = 9, replace = FALSE)
    # 计算CV
    # 计算抽取数据的平均值和标准差
    Dmean <- mean(sampled_data)
    Dsd <- sd(sampled_data)
    
    # Asy functions
    Dmax <- max(sampled_data)
    Dmin <- min(sampled_data)
    
    cv_values <- (Dsd / Dmean)
    asy_values <- ((Dmax - Dmean)/(Dmean - Dmin))
    pv_values<- pv.fun(sampled_data)
    ele.cv.asy.pv.tmp <- data.frame(cbind(ele.tmp, cv_values, asy_values, pv_values))
    
    if (j==1){
      ele.cv.asy.pv = ele.cv.asy.pv.tmp
    }
    else{
      ele.cv.asy.pv = rbind(ele.cv.asy.pv, ele.cv.asy.pv.tmp)
    }
  }
  
  
  if (i==1)
  {
    AGB.cv.asy.pv = ele.cv.asy.pv
  }
  else{
    AGB.cv.asy.pv = rbind(AGB.cv.asy.pv, ele.cv.asy.pv)
  }
}

write.csv(AGB.cv.asy.pv, "../2.DX2013_145sites/data/69.4.AGB.cv.asy.pv.rep5000times.csv")


### AGB43 #### 单独算 4300(只有4个点)
AGB43 = subset(Plant.micro.soil.EMF, Elevation == 4300)

AGB43.la = AGB43[rownames(Plant.rich43),]

# 定义重复实验次数
num_repeats <- 4

AGB.tmp = data.frame(AGB43.la[, "AGB"])
ele = unique(AGB43.la[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  data <- AGB.tmp[which(AGB43.la$Elevation == ele.tmp),]
  
  
  # 内层循环进行蒙特卡洛模拟
  for (j in 1:num_repeats) {
    # 随机抽取3个点
    sampled_data <- sample(data, size = 3, replace = FALSE)
    # 计算CV
    # 计算抽取数据的平均值和标准差
    Dmean <- mean(sampled_data)
    Dsd <- sd(sampled_data)
    
    # Asy functions
    Dmax <- max(sampled_data)
    Dmin <- min(sampled_data)
    
    cv_values <- (Dsd / Dmean)
    asy_values <- ((Dmax - Dmean)/(Dmean - Dmin))
    pv_values<- pv.fun(sampled_data)
    ele.cv.asy.pv.tmp <- data.frame(cbind(cv_values, asy_values, pv_values))
    
    if (j==1){
      ele.cv.asy.pv = ele.cv.asy.pv.tmp
    }
    else{
      ele.cv.asy.pv = rbind(ele.cv.asy.pv, ele.cv.asy.pv.tmp)
    }
  }
  
  
  if (i==1)
  {
    AGB43.cv.asy.pv = ele.cv.asy.pv
  }
  else{
    AGB43.cv.asy.pv = rbind(AGB43.cv.asy.pv, ele.cv.asy.pv)
  }
}

write.csv(AGB43.cv.asy.pv, "../2.DX2013_145sites/data/69.4.AGB4300m.cv.asy.pv.rep4times.csv")



### PLFA 只有3个重复 ####

# 当雄数据每个海拔最多6个点抽取3个点的组合数
n <- 6
k <- 3
combinations <- choose(n, k)

# 输出结果
print(combinations)
# 20

PLFA.Origin = read.csv("../2.DX2013_145sites/data/DX1314.PLFA.orgin.csv",
                       row.names = 1)

PLFA.LA = Plant.micro.soil.EMF[rownames(PLFA.Origin), ]

# 定义重复实验次数
num_repeats <- 20

PLFA.tmp = data.frame(PLFA.LA[, "PLFA"])
ele = unique(PLFA.LA[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  data <- PLFA.tmp[which(PLFA.LA$Elevation == ele.tmp),]
  
  
  # 内层循环进行蒙特卡洛模拟
  for (j in 1:num_repeats) {
    # 随机抽取3个点
    sampled_data <- sample(data, size = 3, replace = FALSE)
    # 计算CV
    # 计算抽取数据的平均值和标准差
    Dmean <- mean(sampled_data)
    Dsd <- sd(sampled_data)
    
    # Asy functions
    Dmax <- max(sampled_data)
    Dmin <- min(sampled_data)
    
    cv_values <- (Dsd / Dmean)
    asy_values <- ((Dmax - Dmean)/(Dmean - Dmin))
    pv_values<- pv.fun(sampled_data)
    ele.cv.asy.pv.tmp <- data.frame(cbind(ele.tmp, cv_values, asy_values, pv_values))
    
    if (j==1){
      ele.cv.asy.pv = ele.cv.asy.pv.tmp
    }
    else{
      ele.cv.asy.pv = rbind(ele.cv.asy.pv, ele.cv.asy.pv.tmp)
    }
  }
  
  
  if (i==1)
  {
    PLFA.cv.asy.pv = ele.cv.asy.pv
  }
  else{
    PLFA.cv.asy.pv = rbind(PLFA.cv.asy.pv, ele.cv.asy.pv)
  }
}

write.csv(PLFA.cv.asy.pv, "../2.DX2013_145sites/data/69.5.PLFA.cv.asy.pv.rep20times.csv")



### Soil.pca1 ####

# 定义重复实验次数
num_repeats <- 5000

Soil.pca1.tmp = data.frame(Plant.micro.soil.EMF[, "Soil.pca1"])
ele = unique(Plant.micro.soil.EMF[,"Elevation"])


i=1
for (i in 1:length(ele)){
  ele.tmp <- ele[i]
  data <- Soil.pca1.tmp[which(Plant.micro.soil.EMF$Elevation == ele.tmp),]
  
  
  # 内层循环进行蒙特卡洛模拟
  for (j in 1:num_repeats) {
    # 随机抽取9个点
    sampled_data <- sample(data, size = 9, replace = FALSE)
    # 计算CV
    # 计算抽取数据的平均值和标准差
    Dmean <- mean(sampled_data)
    Dsd <- sd(sampled_data)
    
    # Asy functions
    Dmax <- max(sampled_data)
    Dmin <- min(sampled_data)
    
    cv_values <- (Dsd / Dmean)
    asy_values <- ((Dmax - Dmean)/(Dmean - Dmin))
    pv_values<- pv.fun(sampled_data)
    ele.cv.asy.pv.tmp <- data.frame(cbind(ele.tmp, cv_values, asy_values, pv_values))
    
    if (j==1){
      ele.cv.asy.pv = ele.cv.asy.pv.tmp
    }
    else{
      ele.cv.asy.pv = rbind(ele.cv.asy.pv, ele.cv.asy.pv.tmp)
    }
  }
  
  
  if (i==1)
  {
    Soil.pca1.cv.asy.pv = ele.cv.asy.pv
  }
  else{
    Soil.pca1.cv.asy.pv = rbind(Soil.pca1.cv.asy.pv, ele.cv.asy.pv)
  }
}

write.csv(Soil.pca1.cv.asy.pv, "../2.DX2013_145sites/data/69.6.Soil.pca1.cv.asy.pv.rep5000times.csv")


### combine all ###

cv.asy.pv.all = rbind(plant.rich.cv.asy.pv, Micro.rich.cv.asy.pv, EMF.cv.asy.pv,
                      AGB.cv.asy.pv, PLFA.cv.asy.pv, Soil.pca1.cv.asy.pv)

Group = data.frame(c(rep("Plant.rich", 45000), rep("Micro.rich", 50000), 
                     rep("EMF", 50000), rep("AGB", 50000),
                     rep("PLFA", 50000), rep("Soil.pca1", 50000)))

cv.asy.pv.all1 = cbind(Group, cv.asy.pv.all)
colnames(cv.asy.pv.all1)[1] = c("Group")

write.csv(cv.asy.pv.all1, "../2.DX2013_145sites/data/69.7.DX1314.cv.asy.pv.rep.all.5000times.csv")


cv.asy.pv.all1.me = melt(cv.asy.pv.all1, id = c("Group", "ele.tmp"))


PV.mean = dcast(cv.asy.pv.all1.me, Group + ele.tmp ~ variable, mean)

PV.sd = dcast(cv.asy.pv.all1.me, Group + ele.tmp ~ variable, sd)

PV.length = dcast(cv.asy.pv.all1.me, Group + ele.tmp ~ variable, length)


PV.mean[is.na(PV.mean)] <- 0

PV.sd[is.na(PV.sd)] <- 0

PV.se.la = PV.sd[, 3:5]/sqrt(PV.length[, 3:5])

PV.se.la1 = cbind(PV.sd[, 1:2], PV.se.la)

PV.mean.me = melt(PV.mean, id = c("Group", "ele.tmp"))
PV.se.la.me = melt(PV.se.la1, id = c("Group", "ele.tmp"))

pv.last = cbind(PV.mean.me, PV.se.la.me$value)
colnames(pv.last)[5] = c("se")

write.csv(pv.last, "../2.DX2013_145sites/data/69.8.DX1314.cv.asy.pv.rep.all.5000times.mean.se.csv")

pv.last.la = subset(pv.last, variable == "pv_values")

write.csv(pv.last.la, "../2.DX2013_145sites/data/69.9.DX1314.pv.rep5000times.mean.se.csv")

pv.last.la$ele.tmp = factor(pv.last.la$ele.tmp)
ggplot(pv.last.la, aes(x=ele.tmp, y=value, color = ele.tmp))+ 
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se))+
  facet_wrap(variable ~ Group, scales="free_y", ncol=3)+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme_bw()+ 
  xlab("Elevation") + ylab("PV")+
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))+
  theme( panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5))+
  theme(legend.position = 'none')
