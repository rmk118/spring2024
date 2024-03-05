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
library(MASS)
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


setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Freshet Saccharina/publication_spring_freshet/metadata_and_imput_data")
path = "C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Freshet Saccharina/publication_spring_freshet/output/alpha-diversity/"

## read in data
sf = readRDS("SF_2021_and_2022_coveraged_rarified001percent.RDS")
sf
metadata = read.csv("metadata_2023_06_21.csv")
metadata = metadata %>% column_to_rownames(var="illumina_id")

sf = phyloseq(sample_data(metadata), otu_table(sf@otu_table), tax_table(sf@tax_table))

sf@sam_data$read_depth_filtered = sample_sums(sf)

## add julian day as a variable
sf@sam_data$julday = as.Date(paste0(sf@sam_data$year, "-", sf@sam_data$month, "-", sf@sam_data$day))
sf@sam_data$julday = as.POSIXlt(sf@sam_data$julday)$yday


##### GET ALPHA DIVERSITY STATS #####
field = subset_samples(sf, lab_or_field_sampled=="field")

## extract data
meta = as.data.frame(as.matrix(field@sam_data))
otu = as.data.frame(t(as.matrix(field@otu_table)))

## format for vegan
otu = as.matrix(otu)

## get species shan for each sample
meta$rich = specnumber(otu)
meta$shan = diversity(otu, index="shannon")
meta$simp = diversity(otu, index="simpson")
meta$shan = diversity(otu, index="shannon")
meta$pielou = meta$shan/(log(meta$shan))


##### RUN LINEAR REGRESSION MODELS FOR ALPHA DIVERSITY ##### 

meta.meri = subset(meta, meta$updated_sample_type=="meristem")
meta.water = subset(meta, meta$updated_sample_type=="water")
meta.rock = subset(meta, meta$updated_sample_type=="rock")

lm.meri = lm(shan~as.numeric(mean.conduct)*as.numeric(mean.temp)*as.numeric(julday), data=meta.meri)
summary(lm.meri)
#plot(lm.meri)

lm.water = lm(shan~as.numeric(mean.conduct)*as.numeric(mean.temp)*as.numeric(julday), data=meta.water)
summary(lm.water)
#plot(lm.water)

lm.rock = lm(shan~as.numeric(mean.conduct)*as.numeric(mean.temp)*as.numeric(julday), data=meta.rock)
summary(lm.rock)
#plot(lm.rock)

##### ALPHA DIVERSITY PLOTS #####
## plot shan
meri=ggplot(meta.meri, aes(x=as.numeric(mean.conduct), y = shan))+
  geom_point(cex=2, alpha=0.3)+
  geom_smooth(se=F, linewidth=1, method="lm", color="black")+
  labs(y="Shannon–Wiener Index", x="Salinity")+
  theme(legend.position = "none")+
  annotate("text", x = 27, y = 5, label = "Field Saccharina", size=4, color="black")+
  ylim(0,6)
meri

water=ggplot(meta.water, aes(x=as.numeric(mean.conduct), y = shan))+
  geom_point(cex=2, alpha=0.3)+
  geom_smooth(se=F, linewidth=1, method="lm", color="black")+
  labs(y="Shannon–Wiener Index", x="Salinity")+
 theme(legend.position = "none")+
  annotate("text", x = 27, y = 1, label = "Field Water", size=4, color="black")+
  ylim(0,6)
water
  

rock=ggplot(meta.rock, aes(x=as.numeric(mean.conduct), y = shan))+
  geom_point(cex=2, alpha=0.3)+
  geom_smooth(se=F, linewidth=1, method="lm", color="black")+
  labs(y="Shannon–Wiener Index", x="Salinity")
rock

ggsave("field_shan_rock.pdf", width=6, height=6, units="in", path = path)

##### ARRANGE PLOTS AND SAVE ######

ggarrange(meri, water, ncol=2, nrow=1, labels=c("E", "F"))
ggsave("field_shan_main.pdf", width=10, height=4.6, units="in", path = path)



