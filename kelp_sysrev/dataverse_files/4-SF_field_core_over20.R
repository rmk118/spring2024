##### SET UP ####
## load packages
library(indicspecies)
library(tidyverse)
library(data.table)
library(phyloseq)
library(plyr)
library(ggpubr)
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

## set working directory and filepath 
setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Freshet Saccharina/publication_spring_freshet/metadata_and_imput_data")
indvalpath = "C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Freshet Saccharina/publication_spring_freshet/output/indval_output/"
path = "C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Freshet Saccharina/publication_spring_freshet/output/indval_output/"

## read in data
sf = readRDS("SF_2021_and_2022_filtered001percent_notrarefied.rds")

metadata = read.csv("metadata_2023_06_21.csv")
metadata = metadata %>% column_to_rownames(var="illumina_id")

nonrare = phyloseq(sample_data(metadata), otu_table(sf@otu_table), tax_table(sf@tax_table))

nonrare@sam_data$read_depth_filtered = sample_sums(nonrare)

field = subset_samples(nonrare, lab_or_field_sampled=="field")


##### INDVAL SUBSTRATE CALCULATIONS at ASV LEVEL ######
## get the metadata and otu data together for indval
metadata = as.data.frame(as.matrix(field@sam_data)) %>% rownames_to_column(var="rownames")
otu = as.data.frame(t(as.matrix(field@otu_table))) %>% rownames_to_column(var="rownames")
tax = as.data.frame(field@tax_table) %>% rownames_to_column(var="ASV")

## merge metadata and otu 
metaasv = left_join(metadata, otu)

## only keep salinity over or equal to 25 to find the core 
metaasv$mean.conduct = as.numeric(metaasv$mean.conduct)
saldf = subset(metaasv, metaasv$mean.conduct>=20)

## get metadata size
medatada.cols = ncol(metadata)

## Get otu table where samples are rows and ASVs are columns 
otu.sal = saldf[,-c(1:medatada.cols)]

## indval comparison 
substrate <- as.character(saldf$updated_sample_type) 

## run indval
indval <- multipatt(otu.sal, substrate, duleg = TRUE, control = how(nperm=999))

# Get indval statistic
indval.16.str <- as.data.frame(indval$str)
indval.16.str$rn <- rownames(indval.16.str)

# get p-value
indval.stat <- as.data.frame(indval$sign) #get dataframe of indval statistic
indval.stat$rn <- rownames(indval.stat) # make column of ASVs

# Prevalence as dataframe
indval.prev <- as.data.frame(indval$A)
# extract rownames into column
setDT(indval.prev, keep.rownames = TRUE)[] 
## rename columns 
colnames(indval.prev) <- paste0("prev.", colnames(indval.prev))
names(indval.prev)[names(indval.prev) == 'prev.rn'] <- 'rn'

# Fidelity as dataframe
indval.fid <- as.data.frame(indval$B) 
# extract rownames into column
setDT(indval.fid, keep.rownames = TRUE)[] 
## rename columns 
colnames(indval.fid) <- paste0("fid.", colnames(indval.fid))
names(indval.fid)[names(indval.fid) == 'fid.rn'] <- 'rn'

# Join statistics together (you could do multi_join but this might crash your computer)
str.and.stat = full_join(indval.16.str, indval.stat,
                         by="rn")
prev.and.fid = full_join(indval.prev, indval.fid,
                         by="rn")
indval_table = full_join(str.and.stat, prev.and.fid,
                         by="rn")

## rename columns to join with taxonomy
names(indval_table)[names(indval_table) == 'rn'] <- 'ASV'

## merge with taxonomy
indval_table= inner_join(indval_table, tax)

coresig = subset(indval_table, indval_table$index != "NaN" & indval_table$p.value<=0.05)
    
    


##### FILTER INDVAL OUTPUT TO ONLY KEEP MERISTEM ASSOCIATED BY 0.8 #####

## filter indval output
sig = subset(coresig, coresig$s.meristem=="1" & 
                   coresig$stat>=0.7)



##### FORMAT THE SINGIFICANT ASSOCITIONS FOR ANALYSIS/PLOTTING #####

## get asv list by threshold
#m09 = c(unique(m09$asv_sequence))
m08 = c(unique(sig$asv_sequence))

## make dataframe for entire salinity dataset
meri = subset_samples(field, updated_sample_type=="meristem")
meri = dephyloseq(meri)

## use significant indicator lists to subset the meristem dataframe
#meri09 = subset(meri, meri$asv_sequence %in% c(m09))
meri08 = subset(meri, meri$asv_sequence %in% c(m08))

## get the dataframe of the core for the paper 
coredf = ddply(meri08, c("asv_id", "Order", "Genus", "asv_sequence"),
               summarise,
               ncore = length(year))
write.csv(coredf, paste0(indvalpath, "over20_core_asv.csv"))

## make presence/abscence column to sum when summarizing 
meri08$presabs = ifelse(meri08$asv_abundance>0, "1", "0")

#### dotplot of EACH CORE ASV  by salinity #####
meriybasv = ddply(meri08, c("Row.names","mean.conduct", "mean.temp", "read_depth_filtered", "Genus", "asv_id"),
               summarise,
               corecount = sum(as.numeric(presabs)),
               coresum = sum(asv_abundance)) %>%
  mutate(corera = as.numeric(coresum)/as.numeric(read_depth_filtered))


coregenus=ggplot(meriybasv, aes(x=as.numeric(mean.conduct), y = corera))+
  geom_point(cex=2, alpha=0.5)+
  geom_smooth(se=F, linewidth=1, method="lm")+
  labs(x="Salintiy", y = "")+
  #annotate("text", x =14, y = 0.2, label = "r = 0.2290, p < 0.001", size=4, color="black")+
  ylim(0,1)+
  facet_wrap(.~Genus)
coregenus

##### CORR TEST FOR EACH ASV #####
ver = subset(meriybasv, meriybasv$asv_id=="ASV1")
cor.test(as.numeric(ver$mean.conduct), ver$corera)

rgam = subset(meriybasv, meriybasv$asv_id=="ASV15")
cor.test(as.numeric(gam$mean.conduct), gam$corera)

cau = subset(meriybasv, meriybasv$asv_id=="ASV3")
cor.test(as.numeric(cau$mean.conduct), cau$corera)

mar = subset(meriybasv, meriybasv$asv_id=="ASV4")
cor.test(as.numeric(mar$mean.conduct), mar$corera)

coc = subset(meriybasv, meriybasv$asv_id=="ASV5")
cor.test(as.numeric(coc$mean.conduct), coc$corera)

rob = subset(meriybasv, meriybasv$asv_id=="ASV7")
cor.test(as.numeric(rob$mean.conduct), rob$corera)

#### summarise the data to count the relative abundance of core per #####

meri08 = ddply(meri08, c("Row.names","mean.conduct", "mean.temp", "read_depth_filtered"),
               summarise,
               corecount = sum(as.numeric(presabs)),
               coresum = sum(asv_abundance)) %>%
  mutate(corera = as.numeric(coresum)/as.numeric(read_depth_filtered))


##### CORRELATION TEST OF CORE BY SALINITY #####
cor.test(as.numeric(meri08$mean.conduct), meri08$corera)


#### dotplot of ALL CORE relative abundance by salinity #####
corera=ggplot(meri08, aes(x=as.numeric(mean.conduct), y = corera))+
  geom_point(cex=2, alpha=0.5)+
  geom_smooth(se=F, linewidth=1, method="lm")+
  labs(x="Salintiy", y = "Relative abundance of core ASV
as percent of total reads")+
  annotate("text", x =14, y = 0.8, label = "r = 0.2290, p < 0.001", size=4, color="black")+
  ylim(0,1)





##### ARRANGE PLOTS ####
ggarrange(corera,coregenus, 
          labels=c('A', "B"),
          widths=c(0.7, 1))

ggsave("field_core_ra_count_over20.pdf", path=path, width = 17.9, height=5.9, units="in")




