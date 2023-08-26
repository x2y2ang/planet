
rm(list=ls())

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# library
library(vegan)
library(FactoMineR)
library(lavaan)
library(semPlot)
library(AICcmodavg)
library(data.table)
library(MASS)


all.micro = read.csv("../2.DX2013_145sites/data/Micro.rich.inter.biomass.combine.AI.group0.69.csv", 
                     row.names = 1)


### A.above0.69 ###
A.above0.69 = subset(all.micro, Group == "Above")

lm.micro.div = lm(Micro.mass ~ Micro.rich, A.above0.69)
qu.micro.div = lm(Micro.mass ~ Micro.rich + I(Micro.rich^2), A.above0.69)
AIC(lm.micro.div, qu.micro.div)
# df      AIC
# lm.micro.div  3 653.7363
# qu.micro.div  4 655.3997

summary(lm.micro.div)
# Adjusted R-squared:  0.02116 , p-value: 0.07741


### B.below0.69 ###
B.below0.69 = subset(all.micro, Group == "Below")

lm.micro.div = lm(Micro.mass ~ Micro.rich, B.below0.69)
qu.micro.div = lm(Micro.mass ~ Micro.rich + I(Micro.rich^2), B.below0.69)
AIC(lm.micro.div, qu.micro.div)
# df      AIC
# lm.micro.div  3 1119.666
# qu.micro.div  4 1121.427

summary(lm.micro.div)
# Adjusted R-squared:  -0.008347 , p-value: 0.9353


envi.qu.F1 = all.micro$Group %in% c("Above", "Below")



micro.div.PLFA = 
  ggplot(all.micro, aes(x=Micro.rich, y=Micro.mass, color = Group))+ 
  
  geom_point(size = 2)+
  # (method="lm", color = "black")+
  # scale_colour_manual(values = c(brewer.pal(10,"Paired")))+
  
  stat_smooth(method="lm",formula=y~x, size=1, lty =2,
              se = T, data = all.micro[envi.qu.F1,])+
  
  theme_bw()+  
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
  xlab("Micro.rich")+ylab("PLFA")

micro.div.PLFA


##########  Micro.net and PLFA ###########

### A.above0.69 ###
lm.Micro.net = lm(Micro.mass ~ Micro.net, A.above0.69)
qu.Micro.net = lm(Micro.mass ~ Micro.net + I(Micro.net^2), A.above0.69)
AIC(lm.Micro.net, qu.Micro.net)
# df      AIC
# lm.Micro.net  3 653.2041
# qu.Micro.net  4 654.2791

summary(lm.Micro.net)
# Adjusted R-squared:  0.02626, p-value: 0.05649


### B.below0.69 ###
lm.Micro.net = lm(Micro.mass ~ Micro.net, B.below0.69)
qu.Micro.net = lm(Micro.mass ~ Micro.net + I(Micro.net^2), B.below0.69)
AIC(lm.Micro.net, qu.Micro.net)
# df      AIC
# lm.Micro.net  3 1113.631
# qu.Micro.net  4 1113.500

summary(qu.Micro.net)
# Adjusted R-squared:  0.04947, p-value: 0.01859

net.lm.F1 = all.micro$Group %in% c("Above")
net.qu.T1 = all.micro$Group %in% c("Below")

Micro.net.PLFA = 
  ggplot(all.micro, aes(x=Micro.net, y=Micro.mass, color = Group))+ 
  
  geom_point(size = 2)+
  
  stat_smooth(method="lm",formula=y~x,size=1, lty=2,
              se = T, data = all.micro[net.lm.F1,])+
  
  stat_smooth(method="lm",formula=y~x+I(x^2),size=1, 
              se = T, data = all.micro[net.qu.T1,])+
  theme_bw()+  
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
  xlab("Micro.net")+ylab("PLFA")

Micro.net.PLFA

library(cowplot)

plot_grid(micro.div.PLFA, Micro.net.PLFA)


ggsave("../2.DX2013_145sites/final.figures2023/Figure.Micro.rich.net.PLFA.above.below0.69.pdf",
       width = 7.5, height = 3)




lm.Micro.net.all = lm(Micro.mass ~ Micro.net, all.micro)
qu.Micro.net.all = lm(Micro.mass ~ Micro.net + I(Micro.net^2), all.micro)
AIC(lm.Micro.net.all, qu.Micro.net.all)


ggplot(all.micro, aes(x=Micro.net, y=Micro.mass))+ 
  geom_point(size = 2)+
  stat_smooth(method="lm",formula=y~x+I(x^2),size=1)+
  theme_bw()+  
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
  xlab("Micro.net")+ylab("PLFA")

