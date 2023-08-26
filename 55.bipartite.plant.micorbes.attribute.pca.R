
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(reshape2)
library(ggplot2)
library(FactoMineR)

BI.AI = read.csv('../2.DX2013_145sites/data/Bipartite.total.attribute.csv')
colnames(BI.AI)[1:2] = c("AI", "Elevation")
BI.AI.tmp = BI.AI[, 3:16]

BI.pca = PCA(BI.AI.tmp, scale.unit=T, graph=F)
BI.pca 
BI.pca$ind$coord
BI.pca$eig

# > BI.pca$eig
# eigenvalue percentage of variance cumulative percentage of variance
# comp 1 9.351236016            66.79454297                          66.79454
# comp 2 2.991018704            21.36441932                          88.15896

BI.pca.t <- BI.pca$ind$coord[,1:2]
colnames(BI.pca.t) = c("Bi.pca1", "Bi.pca2")

BI.LAST = cbind(BI.AI[, c("Elevation","AI")], 
                  BI.pca.t)

gp = c(rep("Above", 4),rep("Below", 5))


BI.LAST = cbind(BI.LAST, gp)
# 这里要用elevation

BI.LAST$Elevation <- as.factor(BI.LAST$Elevation)
BI.LAST$AI <- as.factor(BI.LAST$AI)

ggplot(BI.LAST, aes(x=Bi.pca1, y=Bi.pca2, color = AI)) + # shape = Elevation 
  scale_shape_manual(values= c(1:10))+
  geom_point(size=3)+
  theme_bw()+ 
  guides(col = guide_legend(nrow = 5))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"))+
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
  stat_ellipse(aes(fill=BI.LAST$gp))

ggsave("../2.DX2013_145sites/figures/Bipartite.plant.micorbes.attribute.pca.all.20230526.pdf",
       width = 6, height = 5)


ggplot(BI.LAST, aes(x=Bi.pca1, y=Bi.pca2, color = gp))+ 
  geom_point(alpha=.7, size=2)+
  geom_polygon(data = BI.LAST, alpha = 0.4, 
               aes(fill = factor(gp)),
               show.legend = T)+ 
  # facet_wrap( ~ Domain, scales="free_y", ncol=2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"))+
  # panel.grid.minor = element_blank(),
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black'))+
  theme(text = element_text(size = 10))+
  theme(plot.title = element_text(hjust = 0.5))


ggsave("../2.DX2013_145sites/figures/Bipartite.plant.micorbes.attribute.pca.all.20230609.pdf",
       width = 3.5, height = 3)


### 跨界网络，与plant-plant, microbe-microbe相互作用的关系？###


micro.coh.gp = read.csv("../2.DX2013_145sites/data/micro1314.coh.last.gp.csv",
                  row.names = 1)


micro.mean = dcast(micro.coh.gp, Elevation ~ variable, mean)
micro.mean1 = micro.mean[-1, ]

micro.coh.Bi = cbind(BI.AI[, 1:2], BI.pca.t, micro.mean1[, 2:5])

colnames(micro.coh.Bi)[5:8] = c("Micro.neg", "Micro.neg.pos", "Micro.pos",
                                "Micro.total") 

cor.test(micro.coh.Bi$Bi.pca1, micro.coh.Bi$Micro.total, method = "pearson")


lm.micro = lm(micro.coh.Bi$Bi.pca1 ~ micro.coh.Bi$Micro.total, micro.coh.Bi)
summary(lm.micro)

qu.micro = lm(micro.coh.Bi$Bi.pca1 ~ micro.coh.Bi$Micro.total + 
                I((micro.coh.Bi$Micro.total)^2), micro.coh.Bi)
AIC(lm.micro, qu.micro)

summary(qu.micro)

ggplot(micro.coh.Bi, aes(x=Micro.total, y=Bi.pca1, color = AI))+ 
  geom_point(size=2)+
  geom_smooth(method = "lm", formula = y ~ x+I(x^2))+
  theme_bw()


plant.coh.all = read.csv("../2.DX2013_145sites/data/plant1314.188sites.coverage.cohesion.all.melt.gp.csv",
                   row.names = 1)


plant.mean = dcast(plant.coh.all, Elevation ~ variable, mean)

plant.mean1 = plant.mean[-1, ]

plant.coh.Bi = cbind(BI.AI[, 1:2], BI.pca.t, plant.mean1[, 2:5])

cor.test(plant.coh.Bi$Bi.pca1, plant.coh.Bi$Total.coh, method = "kendall")


plot(plant.coh.Bi$Total.coh, plant.coh.Bi$Bi.pca1)

lm.plant = lm(plant.coh.Bi$Bi.pca1 ~ plant.coh.Bi$Total.coh, plant.coh.Bi)
summary(lm.plant)

qu.plant = lm(plant.coh.Bi$Bi.pca1 ~ plant.coh.Bi$Total.coh + 
                I((plant.coh.Bi$Total.coh)^2), plant.coh.Bi)
AIC(lm.plant, qu.plant)

summary(qu.plant)

plant.coh.Bi$AI = as.factor(plant.coh.Bi$AI)
colnames(plant.coh.Bi)

ggplot(plant.coh.Bi, aes(x=Total.coh, y=Bi.pca1, color = AI))+ 
  geom_point(size=2)+
  geom_smooth(method = "lm", formula = y ~ x+I(x^2))+
  theme_bw()+
  scale_shape_manual(values= c(1:10))+
  # guides(col = guide_legend(nrow = 5))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"))+
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
        axis.title.y=element_text(size=12))


## combine plant and micro ##
plant.micro = cbind(plant.coh.Bi[, c(1:3, 8)], micro.coh.Bi[, c(8)])
colnames(plant.micro)[4:5] = c("Plant.total.coh", "Micro.total.coh")

plant.micro.me = melt(plant.micro, id = c("AI", "Elevation", "Bi.pca1"))

AI = data.frame(rep(c(BI.AI$AI), 2))

plant.micro.me1 = cbind(AI, plant.micro.me)
colnames(plant.micro.me1)[1] = c("Aridity")

library(RColorBrewer)
ggplot(plant.micro.me1, aes(x=Bi.pca1, y=value, color = AI))+ 
  geom_point(size=2)+
  geom_smooth(method = "lm", formula = y ~ x+I(x^2), color = "black")+
  facet_wrap( ~ variable, ncol = 2, scales  = "free_y") +
  theme_bw()+
  # scale_shape_manual(values= c(1:10))+
  scale_shape_manual(values = c(1:19))+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"))+
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
        axis.title.y=element_text(size=12))

ggsave("../2.DX2013_145sites/final.figures2023/Fig.S.Bipartition.plant.micro.coh.relationship.pdf",
       width = 6.5, height = 3)




ggplot(plant.micro.me1, aes(x=Bi.pca1, y=value, color = AI))+ 
  geom_point(size=2)+
  geom_smooth(method = "lm", formula = y ~ x, lty = 2, color = "black")+
  facet_wrap( ~ variable, ncol = 2, scales  = "free_y") +
  theme_bw()+
  # scale_shape_manual(values= c(1:10))+
  scale_shape_manual(values = c(1:19))+
  scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"))+
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
        axis.title.y=element_text(size=12))

ggsave("../2.DX2013_145sites/final.figures2023/Fig.S.lm.Bipartition.plant.micro.coh.relationship.pdf",
       width = 6.5, height = 3)
