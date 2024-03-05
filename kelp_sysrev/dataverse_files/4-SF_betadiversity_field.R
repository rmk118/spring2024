##### WORKSPACE SETUP #####
## load pacakges
library(phyloseq)
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(ggpubr)
library(MASS)
library(plyr)
library(microbiome)
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
site.colors = c("seagreen1", "#39B54A", "darkgreen", "plum", "#A30262")


## tell R where to get data

## read in data
sf = readRDS("SF_2021_and_2022_coveraged_rarified001percent.RDS")

sf@sam_data$read_depth_filtered = sample_sums(sf)

## add julian day as a variable
sf@sam_data$julday = as.Date(paste0(sf@sam_data$year, "-", sf@sam_data$month, "-", sf@sam_data$day))
sf@sam_data$julday = as.POSIXlt(sf@sam_data$julday)$yday

##### RUN AND ANOVA AND MAKE NMDS FOR ALL SUBSTRATES #####
sf = subset_samples(sf, lab_or_field_sampled=="field")

## run permanova
## remove empty taxa
all.field = subset_taxa(sf, taxa_sums(sf) > 0)
# get data frames
otu.all.field = as.data.frame(t(as.matrix(all.field@otu_table)))
meta.all.field = as.data.frame(as.matrix(all.field@sam_data))
## merge the data
all.field.mo = merge(meta.all.field, otu.all.field, by=0)
## count colum numbers
metacols = ncol(meta.all.field)+1
## run permanova with marginal effects 
perm.all.field= adonis2(all.field.mo[,-c(1:metacols)] ~ as.character(updated_sample_type), data=all.field.mo,  method = "bray")
perm.all.field

## post hoc test
all.pairwise = pairwise.adonis(all.field.mo[,-c(1:metacols)], as.character(all.field.mo$updated_sample_type))
all.pairwise

## betadispersion tests
bdt <- phyloseq::distance(all.field, method = "bray")
sample_df <- data.frame(sample_data(all.field))

# salinity
samps <- betadisper(bdt, sample_df$updated_sample_type) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bsam=permutest(samps) 
bsam

#### NMDS 
field.ord = ordinate(all.field, "NMDS","bray")
field.ord

plot_ordination(all.field, field.ord, color="updated_sample_type")+
  geom_point(size=3)+
  scale_color_manual(values=c("#800000","#DAA520","#C0C0C0","#0000FF"))+
  labs(color = "")

ggsave("field_all_samples_NMDS.pdf", 
       width = 6, height = 6, units="in", path=fieldpath)


##### RUN A MARGINAL EFFECTS PERMANOVA - MERISTEMS ####
## only keep meristems
meri = subset_samples(sf, updated_sample_type=="meristem" & mean.temp!="0" & mean.conduct !="0")
## remove empty taxa
meri = subset_taxa(meri, taxa_sums(meri) > 0)
# get data frames
otu.meri = as.data.frame(t(as.matrix(meri@otu_table)))
meta.meri = as.data.frame(as.matrix(meri@sam_data))
## merge the data
meri.mo = merge(meta.meri, otu.meri, by=0)
## count colum numbers
metacols = ncol(meta.meri)+1
## run permanova with marginal effects 
marginal.meristems= adonis2(meri.mo[,-c(1:metacols)] ~ as.numeric(mean.conduct)+as.numeric(mean.temp)+as.numeric(julday)+as.character(site_id), 
                            data=meri.mo, by = "margin", method = "bray")
marginal.meristems

## run permanova with salinity
sal.meristems= adonis2(meri.mo[,-c(1:metacols)] ~ as.numeric(mean.conduct), data=meri.mo, method = "bray")
sal.meristems

## run permanova with temperature
temp.meristems= adonis2(meri.mo[,-c(1:metacols)] ~ as.numeric(mean.temp), data=meri.mo, method = "bray")
temp.meristems

## run permanova with Julian Day
date.meristems= adonis2(meri.mo[,-c(1:metacols)] ~ as.numeric(julday), data=meri.mo, method = "bray")
date.meristems

## run permanova with site 
site.meristems= adonis2(meri.mo[,-c(1:metacols)] ~ as.character(site_id), data=meri.mo, method = "bray")
site.meristems
site.pw=pairwise.adonis(meri.mo[,-c(1:metacols)], meri.mo$site_id)
site.pw

## betadispersion tests
bdt <- phyloseq::distance(meri, method = "bray")
sample_df <- data.frame(sample_data(meri))

# salinity
sal <- betadisper(bdt, sample_df$mean.conduct) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bsal=permutest(sal) 
bsal

# temp
temp <- betadisper(bdt, sample_df$mean.temp) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
btmp=permutest(temp) 
btmp

# julian date
jd <- betadisper(bdt, sample_df$julday) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bjd=permutest(jd) 
bjd

##### RUN A MARGINAL EFFECTS PERMANOVA - WATER ####
## only keep water
water = subset_samples(sf, updated_sample_type=="water" & mean.temp!="0" & mean.conduct !="0")
## remove empty taxa
water = subset_taxa(water, taxa_sums(water) > 0)
# get data frames
otu.water = as.data.frame(t(as.matrix(water@otu_table)))
meta.water = as.data.frame(as.matrix(water@sam_data))
## merge the data
water.mo = merge(meta.water, otu.water, by=0)
## count colum numbers
metacols = ncol(meta.water)+1
## run permanova with marginal effects 
marginal.water= adonis2(water.mo[,-c(1:metacols)] ~ as.numeric(mean.conduct)+as.numeric(mean.temp)+as.numeric(julday)+as.character(site_id), 
                            data=water.mo, by = "margin", method = "bray")
marginal.water

## run permanova with salinity
sal.water= adonis2(water.mo[,-c(1:metacols)] ~ as.numeric(mean.conduct), data=water.mo, method = "bray")
sal.water

## run permanova with temperature
temp.water= adonis2(water.mo[,-c(1:metacols)] ~ as.numeric(mean.temp), data=water.mo, method = "bray")
temp.water

## run permanova with Julian Day
date.water= adonis2(water.mo[,-c(1:metacols)] ~ as.numeric(julday), data=water.mo, method = "bray")
date.water

## run permanova with site 
site.water= adonis2(water.mo[,-c(1:metacols)] ~ as.character(site_id), data=water.mo, method = "bray")
site.water
site.pw=pairwise.adonis(water.mo[,-c(1:metacols)], water.mo$site_id)
site.pw



## betadispersion tests
bdt <- phyloseq::distance(water, method = "bray")
sample_df <- data.frame(sample_data(water))

# salinity
sal <- betadisper(bdt, sample_df$mean.conduct) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bsal=permutest(sal) 
bsal

# temp
temp <- betadisper(bdt, sample_df$mean.temp) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
btmp=permutest(temp) 
btmp

# julian date
jd <- betadisper(bdt, sample_df$julday) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bjd=permutest(jd) 
bjd


##### RUN A MARGINAL EFFECTS PERMANOVA - ROCK ####
## only keep rock
rock = subset_samples(sf, updated_sample_type=="rock" & mean.temp!="0" & mean.conduct !="0")
## remove empty taxa
rock = subset_taxa(rock, taxa_sums(rock) > 0)
# get data frames
otu.rock = as.data.frame(t(as.matrix(rock@otu_table)))
meta.rock = as.data.frame(as.matrix(rock@sam_data))
## merge the data
rock.mo = merge(meta.rock, otu.rock, by=0)
## count colum numbers
metacols = ncol(meta.rock)+1
## run permanova with marginal effects 
marginal.rock= adonis2(rock.mo[,-c(1:metacols)] ~ as.numeric(mean.conduct)+as.numeric(mean.temp)+as.numeric(julday)+site_id, 
                       data=rock.mo, by = "margin", method = "bray")
marginal.rock

## run permanova with salinity
sal.rock= adonis2(rock.mo[,-c(1:metacols)] ~ as.numeric(mean.conduct), data=rock.mo, method = "bray")
sal.rock

## run permanova with temperature
temp.rock= adonis2(rock.mo[,-c(1:metacols)] ~ as.numeric(mean.temp), data=rock.mo, method = "bray")
temp.rock

## run permanova with julian day
date.rock= adonis2(rock.mo[,-c(1:metacols)] ~ as.numeric(julday), data=rock.mo, method = "bray")
date.rock

## run permanova with site 
site.rock= adonis2(rock.mo[,-c(1:metacols)] ~ as.character(site_id), data=rock.mo, method = "bray")
site.rock
site.pw=pairwise.adonis(rock.mo[,-c(1:metacols)], rock.mo$site_id)
site.pw


## betadispersion tests
bdt <- phyloseq::distance(rock, method = "bray")
sample_df <- data.frame(sample_data(rock))

# salinity
sal <- betadisper(bdt, sample_df$mean.conduct) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bsal=permutest(sal) 
bsal

# temp
temp <- betadisper(bdt, sample_df$mean.temp) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
btmp=permutest(temp) 
btmp

# julian date
jd <- betadisper(bdt, sample_df$julday) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bjd=permutest(jd) 
bjd

##### NMDS PLOT - MERSITEM #####
meri.ord = ordinate(meri, "NMDS","bray")
meri.ord

meri@sam_data$mean.conduct = as.numeric(meri@sam_data$mean.conduct)
meri@sam_data$mean.temp = as.numeric(meri@sam_data$mean.temp)

ms=plot_ordination(meri, meri.ord, color="mean.conduct")+
  geom_point(size=3)+
  scale_color_gradient(low="cyan", high="navy")+
  labs(color = "Salinity")+
  annotate("text", x = -2.2, y = 2.1, label = "Field S. latissima
Marginal PERMANOVA, salinity 
F = 1.9115, R^2 = 0.00656, p = 0.023", size=4, color="black")

ms
# Extract the legend. Returns a gtable
sal <- get_legend(ms)
# Convert to a ggplot and print
salleg=as_ggplot(sal)
# remove legend for real plot
ms = ms + theme(legend.position = "none")


##### NMDS PLOT - water #####
water.ord = ordinate(water, "NMDS","bray")
water.ord

water@sam_data$mean.conduct = as.numeric(water@sam_data$mean.conduct)
water@sam_data$mean.temp = as.numeric(water@sam_data$mean.temp)

ws=plot_ordination(water, water.ord, color="mean.conduct")+
  geom_point(size=3)+
  scale_color_gradient(low="cyan", high="navy")+
  theme(legend.position = "none")+
  annotate("text", x = 1.2, y = -0.95, label = "Field Water
Marginal PERMANOVA, salinity 
F = 3.1194, R^2 = 0.02957, p = 0.001", size=4, color="black")
ws

##### NMDS PLOT - rock #####
rock.ord = ordinate(rock, "NMDS","bray")
rock.ord

rock@sam_data$mean.conduct = as.numeric(rock@sam_data$mean.conduct)
rock@sam_data$mean.temp = as.numeric(rock@sam_data$mean.temp)

rs=plot_ordination(rock, rock.ord, color="mean.conduct")+
  geom_point(size=3)+
  scale_color_gradient(low="cyan", high="navy")+
  theme(legend.position = "none")+
  annotate("text", x = -0.7, y = 0.9, label = "Field Rock
Marginal PERMANOVA, salinity 
F = 0.7986, R^2 = 0.008, p = 0.857", size=4, color="black")
rs


ggsave(filename = ("field_rock_NMDS.pdf"), 
       width=5.3, height=5.3, units="in", 
       path=fieldpath)

salleg
##### ARRANGE AND SAVE NMDS PLOTS #####
ggarrange(ms, ws, salleg,
         nrow=1, ncol=3, 
         labels=c("C", "D"," "),
         widths = c(1, 1, 0.5))

ggsave(filename = ("field_substrates_abitic_NMDS.pdf"), 
       width=17, height=5.3, units="in",  
       path=fieldpath)

##### BETA DIVERSITY FOR EACH SAMPLE GROUP #####
## extract data
field = subset_samples(sf, lab_or_field_sampled=="field")

## make list to cycle though
field@sam_data$betalist = paste0(field@sam_data$updated_sample_type, "-", 
                                 field@sam_data$sampling_round, "-", 
                                 field@sam_data$site_id)
betalist.loop = c(unique(field@sam_data$betalist))


## dummy dataframe
beta.df = NULL

## for loop to calculate mean betadiversity for each group
for(i in betalist.loop){
  for(j in i){
    
    #### https://microbiome.github.io/tutorials/Betadiversity.html
    betasub = as.data.frame(divergence(subset_samples(field, betalist %in% c(j)),
                                   apply(abundances(subset_samples(field, betalist %in% c(j))),1, median))) 

    
    ## fix name
    names(betasub) = c("divergence")
    
    ## extract rownames
    betasub$illumina_id = rownames(betasub)
    
    ## join with dummy dataframe
    beta.df = rbind.fill(beta.df, betasub)
    
    
  }
}

beta.df

##### MERGE DISPERSION WITH METADATA #####
metadata = read.csv("metadata_2023_06_21.csv")

beta.df.join = left_join(beta.df, metadata)


## calculate julian day
## add julian day as a variable
beta.df.join$julday = as.Date(paste0(beta.df.join$year, "-", beta.df.join$month, "-", beta.df.join$day))
beta.df.join$julday = as.POSIXlt(beta.df.join$julday)$yday

## summarize by caluclating mean (multiple calculations for same sample)
beta.df = ddply(beta.df.join, c("illumina_id", "julday","mean.conduct", "updated_sample_type", "mean.temp"),
                summarise,
                divergence = mean(divergence))


## separate out by substrate
meri = subset(beta.df, beta.df$updated_sample_type=="meristem")
water = subset(beta.df, beta.df$updated_sample_type=="water")
rock = subset(beta.df, beta.df$updated_sample_type=="rock")

##### RUN CORRELLATION TEST TO SEE IF BETADISPERION IS CORRELATED WITH SALINITY #####
cor.test(as.numeric(meri$mean.conduct), meri$divergence)
cor.test(as.numeric(water$mean.conduct), water$divergence)
cor.test(as.numeric(rock$mean.conduct), rock$divergence)

##### RUN LINEAR MODELS FOR SAL, TEMP, AND JUL DAY ######
lm.meri = lm(divergence~as.numeric(meri$mean.conduct)*as.numeric(meri$mean.temp)*as.numeric(meri$julday), data=meri)
summary(lm.meri)
#plot(lm.meri)

lm.water = lm(divergence~as.numeric(water$mean.conduct)*as.numeric(water$mean.temp)*as.numeric(water$julday), data=water)
summary(lm.water)

lm.rock = lm(divergence~as.numeric(rock$mean.conduct)*as.numeric(rock$mean.temp)*as.numeric(rock$julday), data=rock)
summary(lm.rock)
plot(lm.rock)

#### BETADIVERSITY BY SALINITY PLOTS #####
betameri=ggplot(meri, aes(x=as.numeric(mean.conduct), 
                          y= divergence))+
  geom_point(cex=2, alpha=0.3)+
  geom_smooth(se=F, linewidth=1, method="lm", color="black")+
  labs(x="Salintiy", y = "Bray-Curtis dissimilariy")+
  ylim(0,1)+
  scale_color_manual(values=c(site.colors))+
  annotate("text", x = 27, y = 0.95, label = "Field S. latissima
r=-0.230474, p=<0.001", size=4, color="black")
betameri

betawater=ggplot(water, aes(x=as.numeric(mean.conduct), y= divergence))+
  geom_point(cex=2, alpha=0.3)+
  geom_smooth(se=F, linewidth=1, method="lm", color="black")+
  labs(x="Salintiy", y = " ")+
  ylim(0,1)+
  scale_color_manual(values=c(site.colors))+
  theme(legend.position = "none")+
  annotate("text", x = 27, y = 0.95, label = "Field Water
r=0.05676368, p=0.565", size=4, color="black")
betawater


ggarrange(betameri, betawater, ncol=2, nrow=1, labels=c("E", "F"))
ggsave("field_betadiv_main.pdf", width=10, height=4.6, units="in", path = fieldpath)

ggplot(rock, aes(x=as.numeric(mean.conduct), y= divergence))+
  geom_point(cex=2, alpha=0.3)+
  geom_smooth(se=F, linewidth=1, method="lm", color="black")+
  labs(x="Salintiy", y = "Bray-Curtis dissimilariy")+
  ylim(0,1)

ggsave("field_betadiv_rock.pdf", width = 6, height = 6, units="in", path = fieldpath)



