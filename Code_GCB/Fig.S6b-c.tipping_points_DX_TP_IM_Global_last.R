
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

library(vegan)
library(ggplot2)
library(splitstackshape)
library(segmented)
library(reshape2)
library(ggpubr)
library(RColorBrewer)

CV.Asy = read.csv("../data/GCB/Tipping_point_DX_InnerMongolia_World_drylands_sites.csv")
colnames(CV.Asy)[1] = c("name")


var.list = c("DX", "TP1", "TP2", "IM", "Global ")
CV.Asy$gp=factor(CV.Asy$gp,levels = var.list)
CV.Asy=CV.Asy[order(CV.Asy$gp),]


unique(CV.Asy$name)
var.list = c("Plant.div", "Micro.rich", "AGB", "PLFA", "Soil.pca1", "EF")
CV.Asy$name=factor(CV.Asy$name,levels = var.list)
CV.Asy=CV.Asy[order(CV.Asy$name),]

var.list = c("CV", "Asy")
CV.Asy$Variable=factor(CV.Asy$Variable,levels = var.list)
CV.Asy=CV.Asy[order(CV.Asy$Variable),]

ggplot(CV.Asy, aes(name, value)) + 
  # geom_boxplot(color = "#6b7eb9", size = 0.5) + 
  # geom_point(aes(color = gp), size = 2.5) + # , alpha = 0.7
  facet_wrap( ~ Variable, scales="free_y", ncol=2)+
  scale_color_manual(values = c(
    "#E31A1C", "#FB9A99","#FF7F00", "#B2DF8A", "#00BA38"))+
  
  scale_y_continuous(breaks=c(0.4, 0.5, 0.6, 0.7,
                              0.8, 0.9, 1))+
  theme_bw()+   
  theme(legend.background=element_rect(colour="Black",size=0.5))+
  theme(strip.background = element_blank())+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))+
  theme(text = element_text(size = 12))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))+
  xlab("") + ylab("Aridity")+
  coord_cartesian(ylim = c(0.5,0.9))+
  geom_boxplot(size = 0.5, 
               alpha = 0.3, width = 0.7)+ 
  geom_jitter(aes(color = gp), size = 2.2, alpha = 0.8)


ggsave("../GCB_Figure/Fig.S6b-c.tipping_points_DX_TP_IM_Global_legend.pdf",
       height=3.5, width=6.5)

