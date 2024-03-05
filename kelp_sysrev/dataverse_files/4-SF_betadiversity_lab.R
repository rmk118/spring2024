##### WORKSPACE SETUP #####
## load packages
library(phyloseq)
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(ggpubr)
library(car)
library(plyr)
library(microbiome)
library(rstatix)
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

sf = readRDS("SF_2021_and_2022_coveraged_rarified001percent.RDS")

sf@sam_data$read_depth_filtered = sample_sums(sf)


## add julian day as a variable
sf@sam_data$julday = as.Date(paste0(sf@sam_data$year, "-", sf@sam_data$month, "-", sf@sam_data$day))
sf@sam_data$julday = as.POSIXlt(sf@sam_data$julday)$yday

salinity.colors = c("cyan","dodgerblue2","navy")

##### NMDS PLOT FOR ALL LAB SAMPLES ######
alllab = subset_samples(sf, labnotes_aquarium_day =="D-5")

lab.ord = ordinate(alllab, "NMDS","bray")
lab.ord

plot_ordination(alllab, lab.ord, color="updated_sample_type")+
  geom_point(size=3)+
  scale_color_manual(values=c("#B0C4DE","#DAA520","#0000FF"))+
  labs(color = "")

ggsave("lab_all_samples_NMDS.pdf", 
       width = 6, height = 6, units="in", path=labpath)

##### RUN A PERMANOVA - MERISTEMS ####
## only keep meristems
meri = subset_samples(sf, updated_sample_type=="meristem" & labnotes_aquarium_day %in% c("D-5"))
## remove empty taxa
meri = subset_taxa(meri, taxa_sums(meri) > 0)
# get data frames
otu.meri = as.data.frame(t(as.matrix(meri@otu_table)))
meta.meri = as.data.frame(as.matrix(meri@sam_data))
## merge the data
meri.mo = merge(meta.meri, otu.meri, by=0)
## count colum numbers
metacols = ncol(meta.meri)+1

## run permanova with salinity
perm.meristems= adonis2(meri.mo[,-c(1:metacols)] ~ as.factor(aquarium_salinity), 
                            data=meri.mo, method = "bray")
perm.meristems

## post hoc test
sal.pairwise = pairwise.adonis(meri.mo[,-c(1:metacols)], as.factor(meri.mo$aquarium_salinity))
sal.pairwise

## run permanova with Julian Day
date.meristems= adonis2(meri.mo[,-c(1:metacols)] ~ as.numeric(julday), data=meri.mo, method = "bray")
date.meristems

## betadispersion test
bdt <- phyloseq::distance(meri, method = "bray")
sample_df <- data.frame(sample_data(meri))
# salinity
sal <- betadisper(bdt, sample_df$aquarium_salinity) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bsal=permutest(sal) 
bsal



##### RUN A PERMANOVA - WATER ####
## only keep meristems
water = subset_samples(sf, updated_sample_type=="water" & 
                         aquarium_salinity!="NA" & labnotes_aquarium_day %in% c("D-5"))
## remove empty taxa
water = subset_taxa(water, taxa_sums(water) > 0)
# get data frames
otu.water = as.data.frame(t(as.matrix(water@otu_table)))
meta.water = as.data.frame(as.matrix(water@sam_data))
## merge the data
water.mo = merge(meta.water, otu.water, by=0)
## count colum numbers
metacols = ncol(meta.water)+1
## run permanova with
perm.water= adonis2(water.mo[,-c(1:metacols)] ~ as.factor(aquarium_salinity),
                            data=water.mo, method = "bray")
perm.water

## post hoc test
sal.pairwise = pairwise.adonis(water.mo[,-c(1:metacols)], as.factor(water.mo$aquarium_salinity))
sal.pairwise


## run permanova with Julian Day
date.water= adonis2(water.mo[,-c(1:metacols)] ~ as.numeric(julday), data=water.mo, method = "bray")
date.water


## betadispersion test
bdt <- phyloseq::distance(water, method = "bray")
sample_df <- data.frame(sample_data(water))
# salinity
sal <- betadisper(bdt, sample_df$aquarium_salinity) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bsal=permutest(sal) 
bsal

##### RUN A PERMANOVA - airline ####
## only keep meristems
air = subset_samples(sf, updated_sample_type=="airline" & 
                        aquarium_salinity!="NA" & labnotes_aquarium_day %in% c("D-5"))
## remove empty taxa
air = subset_taxa(air, taxa_sums(air) > 0)

# get data frames
otu.air = as.data.frame(t(as.matrix(air@otu_table)))
meta.air = as.data.frame(as.matrix(air@sam_data))

## merge the data
air.mo = merge(meta.air, otu.air, by=0)

## run permanova with salinity
date.air= adonis2(air.mo[,-c(1:metacols)] ~ as.factor(aquarium_salinity), data=air.mo, method = "bray")
date.air

## betadispersion test
bdt <- phyloseq::distance(air, method = "bray")
sample_df <- data.frame(sample_data(air))
# salinity
sal <- betadisper(bdt, sample_df$aquarium_salinity) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bsal=permutest(sal) 
bsal

##### NMDS PLOT - MERSITEM #####
meri.ord = ordinate(meri, "NMDS","bray")
meri.ord

ms=plot_ordination(meri, meri.ord, color="aquarium_salinity")+
  geom_point(size=3)+
  scale_color_manual(values=c(salinity.colors))+
  labs(color = "Salinity")


# Extract the legend. Returns a gtable
sal <- get_legend(ms)
# Convert to a ggplot and print
salleg=as_ggplot(sal)
# remove legend for real plot
ms = ms + theme(legend.position = "none")


md=plot_ordination(meri, meri.ord, color="julday")+
  geom_point(size=3)+
  scale_color_gradient(low="lightgray", high="black")+
  labs(color = "Julian day")+
  annotate("text", x = -2.3, y = 2, label = "stress = 0.1959144 (Bray)", size=4, color="black")+
  annotate("text", x=-2.3, y=1.7, label="lab Saccharina", size=4, color="black")


# Extract the legend. Returns a gtable
day <- get_legend(md)
# Convert to a ggplot and print
dayleg=as_ggplot(day)
# remove legend for real plot
md = md + theme(legend.position = "none")



##### NMDS PLOT - water #####
water.ord = ordinate(water, "NMDS","bray")
water.ord

ws=plot_ordination(water, water.ord, color="aquarium_salinity")+
  geom_point(size=3)+
  scale_color_manual(values=c(salinity.colors))+
  annotate("text", x = -0.5, y = 0.9, label = "stress = 0.2156382 (Bray)", size=4, color="black")+
  annotate("text", x=-0.5, y=0.8, label="lab water (day 4)", size=4, color="black")+
  theme(legend.position = "none")


wd=plot_ordination(water, water.ord, color="julday")+
  geom_point(size=3)+
  scale_color_gradient(low="lightgray", high="black")+
  annotate("text", x = -0.5, y = 0.9, label = "stress = 0.2156382 (Bray)", size=4, color="black")+
  annotate("text", x=-0.5, y=0.8, label="lab water", size=4, color="black")+
  theme(legend.position = "none")

##### NMDS PLOT - air #####
air.ord = ordinate(air, "NMDS","bray")
air.ord

as=plot_ordination(air, air.ord, color="aquarium_salinity")+
  geom_point(size=3)+
  scale_color_manual(values=c(salinity.colors))+
  theme(legend.position = "none")

as

ad=plot_ordination(air, air.ord, color="julday")+
  geom_point(size=3)+
  scale_color_gradient(low="lightgray", high="black")+
  annotate("text", x = -0.7, y = 0.9, label = "stress = 0.1519201 (Bray)", size=4, color="black")+
  annotate("text", x=-0.7, y=0.7, label="lab airline", size=4, color="black")+
  theme(legend.position = "none")

ggsave(filename = ("lab_airline_d4_NMDS.pdf"), 
       width=5.3, height=5.3, units="in", 
       path=labpath)


##### ARRANGE AND SAVE NMDS PLOTS #####
ggarrange(ms, ws, salleg, 
         # md, wd, ad, dayleg,
         nrow=1, ncol=3, 
         labels=c("A", "B"," "),
         widths = c(1, 1, 0.5))

ggsave(filename = ("lab_main_text_d4_NMDS.pdf"), 
       width=17, height=5.3, units="in", 
       path=labpath)




##### BETA DIVERSITY FOR EACH SAMPLE GROUP #####
## extract data
lab = subset_samples(sf, lab_or_field_sampled=="lab" & labnotes_aquarium_day=='D-5')

## make list to cycle though
betalist.loop = c(unique(lab@sam_data$updated_sample_type))


## dummy dataframe
beta.df = NULL

## for loop to calculate mean betadiversity for each group
for(i in betalist.loop){
  for(j in i){
    
    #### https://microbiome.github.io/tutorials/Betadiversity.html
    b30 = as.data.frame(divergence(subset_samples(lab, aquarium_salinity=="30"),
                     apply(abundances(subset_samples(lab, aquarium_salinity=="30")),1, median))) 
    
    b20 = as.data.frame(divergence(subset_samples(lab, aquarium_salinity=="20"),
                     apply(abundances(subset_samples(lab, aquarium_salinity=="20")),1, median)))
    
    b10 = as.data.frame(divergence(subset_samples(lab, aquarium_salinity=="10"),
                     apply(abundances(subset_samples(lab, aquarium_salinity=="10")),1, median)))
    
    
    ## fix names
    names(b30) = c("divergence")
    names(b20) = c("divergence")
    names(b10) = c("divergence")
    
    ### join together
    all.sals = rbind(b30, b20, b10)
    
    ## extract rownames
    all.sals$sample_id = rownames(all.sals)
    
    ## join with dummy dataframe
    beta.df = rbind.fill(beta.df, all.sals)
    

  }
}

beta.df


## summarize by caluclating mean (multiple calculations for same sample)
beta.df = ddply(beta.df, c("sample_id"),
                summarise,
                divergence = mean(divergence))


### get metadata out of phyloseq
meta = as.data.frame(sf@sam_data)

## merge
beta.df = left_join(beta.df, meta)


##### separate by substrate #####
meribeta = subset(beta.df, beta.df$updated_sample_type=="meristem")
waterbeta = subset(beta.df, beta.df$updated_sample_type=="water")
airbeta = subset(beta.df, beta.df$updated_sample_type=="airline")


##### calculate mean betadiversity per sample, ANOVA and plot - meristem #####
## anova
am = oneway.test(divergence~as.factor(aquarium_salinity), data=meribeta, var.equal = F)
am
posthoc_anova(divergence~as.factor(aquarium_salinity), data=meribeta, method = "Games-Howell")

## plot
meribp=ggplot(meribeta, aes(x=as.factor(aquarium_salinity), y = divergence, 
                     color=as.factor(aquarium_salinity), fill=as.factor(aquarium_salinity)))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  labs(y="Bray-Curtis dissimilarity", x=" ")+
  theme(legend.position = "none")+
  ylim(0,1)+
  annotate("text", x = 2.5, y = 0.1, label = "Lab S. latissima
Welsh ANOVA F(2,15.082)=6.182, p=0.011", size=4, color="black")+
  annotate("segment", x=1.1, xend=1.9, y=0.5, yend=0.5)+
  annotate("segment", x=2.1, xend=2.9, y=0.5, yend=0.5)+
  annotate("segment", x=1.1, xend=2.9, y=0.3, yend=0.3)+
  annotate("text", x=1.5, y=0.55, label="p=0.552")+
  annotate("text", x=2.5, y=0.55, label="p=0.011")+
  annotate("text", x=2, y=0.35, label="p=0.239")
meribp

##### calculate mean betadiversity per sample, ANOVA and plot - water #####
aw = aov(divergence~as.factor(aquarium_salinity), data=waterbeta)
summary(aw)
leveneTest(divergence~as.factor(aquarium_salinity), data=waterbeta)
#plot(aw)
TukeyHSD(aw)

## plot
waterbp=ggplot(waterbeta, aes(x=as.factor(aquarium_salinity), y = divergence, 
                            color=as.factor(aquarium_salinity), fill=as.factor(aquarium_salinity)))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  labs(y="Bray-Curtis dissimilarity", x=" ")+
  theme(legend.position = "none")+
  ylim(0,1)+
  annotate("text", x = 2.5, y = 0.1, label = "Lab Water
ANOVA F(2,27)=9.481, p<0.001", size=4, color="black")+
  annotate("segment", x=1.1, xend=1.9, y=0.5, yend=0.5)+
  annotate("segment", x=2.1, xend=2.9, y=0.5, yend=0.5)+
  annotate("segment", x=1.1, xend=2.9, y=0.3, yend=0.3)+
  annotate("text", x=1.5, y=0.55, label="p=0.403")+
  annotate("text", x=2.5, y=0.55, label="p<0.001")+
  annotate("text", x=2, y=0.35, label="p<0.001")
waterbp

##### calculate mean betadiversity per sample, ANOVA and plot - airline #####
aa = aov(divergence~as.factor(aquarium_salinity), data=airbeta)
summary(aa)
leveneTest(divergence~as.factor(aquarium_salinity), data=airbeta)
#plot(aa)
TukeyHSD(aa)

## plot
airbp=ggplot(airbeta, aes(x=as.factor(aquarium_salinity), y = divergence, 
                              color=as.factor(aquarium_salinity), fill=as.factor(aquarium_salinity)))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  labs(y="Bray-Curtis dissimilarity", x=" ")+
  theme(legend.position = "none")+
  ylim(0,1)+
  annotate("text", x = 2.5, y = 0.1, label = "Airline
ANOVA F(2,12)=2.5601, p=0.119", size=4, color="black")
airbp


##### SAVE PLOTS #####
ggarrange(meribp, waterbp, nrow=1, ncol=2, labels = c("A", "B"))

ggsave("lab_betadispersion_salinity.pdf", 
       width = 10, height = 4.6, units="in", path=labpath)

ggarrange(airbp, nrow=1, ncol=1, labels = c("A", "B"))
ggsave("lab_betadispersion_airline.pdf", 
       width = 6, height = 6, units="in", path=labpath)
