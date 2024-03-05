##### WORKSPACE SETUP #####
library(phyloseq)
library(tidyverse)
library(data.table)
library(vegan)
library(car)
library(plyr)
library(ggpubr)
library(FSA)
library(ggh4x)
library(ggplot2)+theme_set(theme_bw()+
                             theme(strip.background = element_rect(fill="white"),
                                   axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
                                   axis.text.x = element_text(colour = "black", face = "bold", size = 12),
                                   legend.text = element_text(size = 8, face ="bold", colour ="black"),
                                   legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
                                   axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
                                   legend.title = element_text(size = 14, colour = "black", face = "bold"),
                                   legend.key=element_blank(),
                                   axis.ticks = element_blank()
                             ))

dephyloseq = function(phylo_obj){
  
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  
  ## how many metadata columns you have 
  metacols = ncol(meta)+1
  
  ## get out the otu table 
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  
  ## get out the taxonomy file 
  tax = as.data.frame(phylo_obj@tax_table)
  
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="asv_id")
  
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "asv_id", values_to="asv_abundance")
  
  ## Join the metadata and otu table with the taoxnomy table 
  mot = full_join(mo, tax)
  
  ## Specify the output for the dephyloseq funciton 
  output = mot
}

## make color vector for site
salinity.colors = c("cyan", "dodgerblue2", "navy")

## read in data
sf = readRDS("SF_2021_and_2022_coveraged_rarified001percent.RDS")

## add julian day as a variable
sf@sam_data$julday = as.Date(paste0(sf@sam_data$year, "-", sf@sam_data$month, "-", sf@sam_data$day))
sf@sam_data$julday = as.POSIXlt(sf@sam_data$julday)$yday


##### GET ALPHA DIVERSITY STATS #####
lab = subset_samples(sf, aquarium_day %in% c("D-5"))

## extract data
meta = as.data.frame(as.matrix(lab@sam_data))
otu = as.data.frame(t(as.matrix(lab@otu_table)))

## format for vegan
otu = as.matrix(otu)

## get species shannon for each sample
meta$rich = specnumber(otu)
meta$shannon = diversity(otu)


meta$updated_sample_type = factor(meta$updated_sample_type, levels=c("meristem", "water", "airline"))

##### RUN ANOVA FOR shannon ##### 
meta.meri = subset(meta, meta$updated_sample_type=="meristem")
meta.water = subset(meta, meta$updated_sample_type=="water")
meta.air = subset(meta, meta$updated_sample_type=="airline")

## meristm
aov.meri= aov(shannon~aquarium_salinity, data=meta.meri)
summary(aov.meri)
#plot(aov.rich)
leveneTest(aov.meri)

## water
aov.water= aov(shannon~aquarium_salinity, data=meta.water)
summary(aov.water)
#plot(aov.rich)
leveneTest(aov.water)

## airline
aov.air= aov(shannon~aquarium_salinity, data=meta.air)
summary(aov.air)
#plot(aov.rich)
leveneTest(aov.air)

##### PLOT ASV shannon #####

## plot shannon
bp.meri=ggplot(meta.meri, aes(x=aquarium_salinity,
                              y = shannon, color=aquarium_salinity, fill=aquarium_salinity))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  labs(y="Shannon–Wiener Index", x=" ")+
  theme(legend.position = "none")+
  ylim(0,4)+
  annotate("text", x = 3, y = 0.5, label = "Lab Saccharina
F(2,29)=1.294, p=0.289", size=4, color="black")
bp.meri

bp.water=ggplot(meta.water, aes(x=aquarium_salinity,
                                y = shannon, color=aquarium_salinity, fill=aquarium_salinity))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  labs(y="Shannon–Wiener Index", x=" ")+
  theme(legend.position = "none")+
  ylim(0,4)+
  annotate("text", x = 3, y = 0.5, label = "Lab Water
F(2,27)=0.343, p=0.713", size=4, color="black")
bp.water

ggarrange(bp.meri, bp.water, ncol=2, nrow=1, labels=c("A", "B"))
ggsave("lab_shannon_main.pdf", width=10, height=4.6, units="in", path = path)

bp.air=ggplot(meta.air, aes(x=aquarium_salinity,
                            y = shannon, color=aquarium_salinity, fill=aquarium_salinity))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  labs(y="Shannon–Wiener Index", x=" ")+
  theme(legend.position = "none")+
  annotate("text", x = 1, y = 0.5, label = "Lab Airline
F(2,12)=0.291, p=0.753", size=4, color="black")
bp.air
ggsave("lab_shannon_airline.pdf", width=6, height=6, units="in", path = path)


