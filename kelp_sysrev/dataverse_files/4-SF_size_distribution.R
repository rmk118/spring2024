##### WORKSPACE SETUP #####
## load packages
library(plyr)
library(tidyverse)
library(ggpubr)
library(car)
library(ggplot2)+theme_set(theme_bw()+
                             theme(strip.background = element_rect(fill="white"),
                                   axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
                                   axis.text.x = element_text(colour = "black", face = "bold", size = 12),
                                   legend.text = element_text(size = 8, face ="bold", colour ="black"),
                                   legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
                                   axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
                                   legend.title = element_text(size = 14, colour = "black", face = "bold"),
                                   legend.key=element_blank(),
                                   #axis.ticks = element_blank()
                                   ))


## set the position dodge for ggplot
pd = position_dodge(0.2) # move them .05 to the left and right

## make color vector for site
site.colors = c("seagreen1", "#39B54A", "darkgreen", "plum", "#A30262")

## read in data
metadata = read.csv("field_condition_metadata.csv")

## add julian date 
metadata$jul.date = paste0(as.numeric(metadata$year), "-", 
                                  as.numeric(metadata$month), "-", 
                                  as.numeric(metadata$day))
metadata$jul.date = as.Date(metadata$jul.date)
metadata$jul.date = as.POSIXlt(metadata$jul.date)$yday

## order the sites
metadata$site_id = factor(metadata$site_id, levels=c("lhp", "tb","scp","gw","cp"))

##### PLOT WIDTH BY SITE ######
ggplot(metadata, aes(y=width_cm, x = site_id, fill=site_id, color=site_id))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_colour_manual(values = site.colors)+
  scale_fill_manual(values=site.colors)+
  labs(y="Thallus width (cm)", x="Sampling Site")+
  annotate("text", x=1, y=30, label="a", size=4, color="black")+
  annotate("text", x=2, y=30, label="b", size=4, color="black")+
  annotate("text", x=3, y=30, label="b", size=4, color="black")+
  annotate("text", x=4, y=30, label="a", size=4, color="black")+
  annotate("text", x=5, y=30, label="a", size=4, color="black")

##### ANOVA OF WIDTH BY SITE #####
kruskal.test(width_cm~site_id, data=metadata)

pairwise.wilcox.test(metadata$width_cm, metadata$site_id,
                     p.adjust.method = "BH")



