
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

CV.Asy = read.csv("../2.DX2013_145sites/data/Tipping_point_DX_InnerMongolia_World_drylands.csv")
colnames(CV.Asy)[1] = c("name")
CV.Asy = CV.Asy[,-2]


var.list = c("CV", "Asy")
CV.Asy$Variable=factor(CV.Asy$Variable,levels = var.list)
CV.Asy=CV.Asy[order(CV.Asy$Variable),]



###### plot ########
ggboxplot(CV.Asy, x="name", y="value", 
          palette = "jco", color = "name",
          add = "jitter", 
          facet.by = "Variable",
          short.panel.labs = FALSE,
          ncol=2)+ #  scales="free_y",
  scale_color_manual(values = c("#C83B3B","#7FB77E", 
                                "#5B7DB1", "#905E96",
                                "#FFB200", "#7F3B08"))+
  scale_y_continuous(breaks=c(0.4,  0.5,  0.6,  0.7,
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
  xlab("") + ylab("Aridity tipping point")+
  theme(legend.position = 'none')+ 
  coord_cartesian(ylim = c(0.5,0.9))


ggsave("../2.DX2013_145sites/last_figures_and_tables/15.tipping.point.DX.TP.World.Drylands.pdf",
       height=3.5, width=6)

