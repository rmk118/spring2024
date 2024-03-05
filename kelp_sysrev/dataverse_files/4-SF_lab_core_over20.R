##### SET UP ####
## load packages
library(indicspecies)
library(tidyverse)
library(data.table)
library(phyloseq)
library(plyr)
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
                                   # axis.ticks = element_blank()
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

salinity.colors = c("cyan","dodgerblue2","navy")

## read in data
sf = readRDS("SF_2021_and_2022_filtered001percent_notrarefied.rds")

metadata = read.csv("metadata_2023_06_21.csv")
metadata = metadata %>% column_to_rownames(var="illumina_id")

nonrare = phyloseq(sample_data(metadata), otu_table(sf@otu_table), tax_table(sf@tax_table))

nonrare@sam_data$read_depth_filtered = sample_sums(nonrare)

coredf = read.csv(paste0(indvalpath, "over20_core_asv.csv"))


##### SUBSET LAB DATA TO ONLY KEEP DAY FOUR, MERISTEMS, AND CORE ASVS #####
lab = subset_samples(nonrare, updated_sample_type=="meristem" & paper_aquarium_day %in% c("day4"))
lab = subset_taxa(lab, asv_sequence %in%c(coredf$asv_sequence))


##### GET MEAN CORE RA ON DAY 0 #####
fieldd0 = subset_samples(nonrare, updated_sample_type=="meristem" & paper_aquarium_day %in% c("day0"))
fieldd0 = subset_taxa(fieldd0, asv_sequence %in%c(coredf$asv_sequence))
d0 = dephyloseq(fieldd0)

## summarize data
d0sum = ddply(d0, c("Row.names","read_depth_filtered"),
                  summarise,
                  coresum = sum(asv_abundance)) %>%
  mutate(corera = as.numeric(coresum)/as.numeric(read_depth_filtered))

mean(d0sum$corera)
sd(d0sum$corera)


##### PLOT THE LAB DATA #####
## get out of phyloseq
labdf = dephyloseq(lab)

## make presence/abscence column to sum when summarizing 
labdf$presabs = ifelse(labdf$asv_abundance>0, "1", "0")

## replace field day 0 slainity with a salintiy interval 
#labdf$aquarium_salinity = ifelse(labdf$paper_aquarium_day=="day0", "field", labdf$aquarium_salinity)

#### dotplot of EACH CORE ASV  by salinity #####
meriybasv = ddply(labdf, c("Row.names","aquarium_salinity","read_depth_filtered", "Genus", "asv_id"),
                  summarise,
                  corecount = sum(as.numeric(presabs)),
                  coresum = sum(asv_abundance)) %>%
  mutate(corera = as.numeric(coresum)/as.numeric(read_depth_filtered))

coregenusplot=ggplot(meriybasv, aes(x=aquarium_salinity, 
                                    y = corera, 
                                    color=aquarium_salinity,
                                    fill=aquarium_salinity))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  facet_wrap(.~Genus)+
  theme(legend.position = "none")
coregenusplot


##### KRUSKALL WALLIS FOR EACH ASV #####
per = subset(meriybasv, meriybasv$asv_id=="ASV1")
aper=aov(corera~as.factor(aquarium_salinity), data=per)
leveneTest(corera~as.factor(aquarium_salinity), data=per)
#plot(aper)
summary(aper)

gam = subset(meriybasv, meriybasv$asv_id=="ASV15")
agam=aov(corera~as.factor(aquarium_salinity), data=gam)
leveneTest(corera~as.factor(aquarium_salinity), data=gam)
#plot(agam)
summary(agam)

lit = subset(meriybasv, meriybasv$asv_id=="ASV3")
alit=aov(corera~as.factor(aquarium_salinity), data=lit)
leveneTest(corera~as.factor(aquarium_salinity), data=lit)
#plot(alit)
summary(alit)
TukeyHSD(alit)

mar = subset(meriybasv, meriybasv$asv_id=="ASV4")
amar=aov(corera~as.factor(aquarium_salinity), data=mar)
leveneTest(corera~as.factor(aquarium_salinity), data=mar)
#plot(amar)
summary(amar)

coc = subset(meriybasv, meriybasv$asv_id=="ASV5")
acoc=aov(corera~as.factor(aquarium_salinity), data=coc)
leveneTest(corera~as.factor(aquarium_salinity), data=coc)
#plot(acoc)
summary(acoc)

rob = subset(meriybasv, meriybasv$asv_id=="ASV7")
kruskal.test(corera~as.factor(aquarium_salinity), data=rob)
pairwise.wilcox.test(rob$corera, rob$aquarium_salinity, p.adjust.method = "BH")

##### DO RELATIVE ABUNDANCE OF WHOLE CORE BY SLAINITY - PLOT #####

mericore = ddply(labdf, c("Row.names","aquarium_salinity", "read_depth_filtered"),
               summarise,
               corecount = sum(as.numeric(presabs)),
               coresum = sum(asv_abundance)) %>%
  mutate(corera = as.numeric(coresum)/as.numeric(read_depth_filtered))

coreall=ggplot(meriybasv, aes(x=aquarium_salinity, 
                                    y = corera, 
                                    color=aquarium_salinity,
                                    fill=aquarium_salinity))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  theme(legend.position = "none")

coreall

kruskal.test(corera~aquarium_salinity, data=meriybasv)
pairwise.wilcox.test(meriybasv$corera, meriybasv$aquarium_salinity, p.adjust.method = "BH")

##### arrange plots ######
ggarrange(coreall,coregenusplot, 
          labels=c('C', "D"),
          widths=c(0.7, 1))

ggsave("lab_core_ra_count_over20.pdf", path=path, width = 17.9, height=5.9, units="in")



