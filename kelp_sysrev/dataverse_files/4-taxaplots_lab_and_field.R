###### SET UP #######
## load packages
library(phyloseq)
library(tidyverse)
library(plyr)
library(qualpalr)
library(ggpubr)
library(ggh4x)
library(ggplot2); theme_set(theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size=10, face="bold"),
        strip.text = element_text(color="black", size=10),
        legend.text=element_text(size=10),
        axis.line = element_line(colour = "black")))

## load functions
dephyloseq = function(phylo_obj){
  
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  
  ## how many metadata columns you have 
  metacols = ncol(meta)+1
  
  ## get out the otu table 
  ## if your metadta is empty after running this, you need to use 
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  #otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  
  ## get out the taxonomy file 
  tax = as.data.frame(phylo_obj@tax_table)
  
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="ASVid")
  
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "ASVid", values_to="asv_abundance")
  
  ## Join the metadata and otu table with the taoxnomy table 
  mot = full_join(mo, tax)
  
  ## Specify the output for the dephyloseq funciton 
  output = mot
}

## tell R whereto get data
## tell R where to get data
setwd("")
path=""

## read in field data
sf = readRDS("SF_2021_and_2022_filtered001percent_notrarefied.rds")

field = subset_samples(sf, lab_or_field_sampled=="field")

## make field groups
field@sam_data$Group = paste0(field@sam_data$updated_sample_type, field@sam_data$site_id)

lab = subset_samples(sf, paper_aquarium_day %in% c("day1", 'day4'))

## make lab groups
lab@sam_data$Group = paste0(lab@sam_data$aquarium_salinity, lab@sam_data$updated_sample_type, lab@sam_data$paper_aquarium_day)

##### FORMAT DATA FOR PLOTS #####
full = merge_phyloseq(field, lab)

## summarize at rank 6
full = tax_glom(full, taxrank = "Genus")

full@sam_data$rd_filtered = sample_sums(full)

## get data out of phyloseq and inta a dataframe
fulldf = dephyloseq(full)

## caluclate relative abundance of each rank 6 within each sample
fulldf$relativeabundance = as.numeric(fulldf$asv_abundance)/as.numeric(fulldf$rd_filtered)

## make plotnames
fulldf$plotnames = paste0(fulldf$Order, ";", fulldf$Genus)


##### use groups to prep data for plot calculations #####

# make list of groups
grouplist = c(unique(fulldf$Group))

## summarize data by taxaplot geoup type. 
full.sum = ddply(fulldf, c("Group", "plotnames"),
                     summarise,
                     sum = sum(relativeabundance))

## sort data by relative abundance. This is how the loop will pick the mos tabundant taxa
sorted = full.sum[order(-full.sum$sum),]

##### CALULCUATE MORE ABUNDANT TAXA ######
## make empty dataframe to store output from the loop
top.df = NULL

## start loop
for(i in grouplist) {
  for(j in i) {
    
    ## subset dataframe by samples
    #!# Remeber to change te substrate to your group! 
    sample = subset(sorted, sorted$Group %in% c(j))
    
    ## get top 15 genera
    top = sample[c(1:15),]
    
    ## save list of top  abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)
    
    ## close loop 
  }
}

## add identifier for top 15 taxa
top.df$place = "top_15"

###### FORMAT TOP TAXA OUTPUT #####
## join the top taxa and existing dataframe
alldata = full_join(fulldf, top.df)

## make the empty "place" cells say bottom. This workes because we used full_join
alldata$place = replace(alldata$place, is.na(alldata$place), "bottom")

## replace plot_names that have bottom taxa as their "place" with Other
alldata[alldata$place == "bottom",]$plotnames <- "Others"

###### GET COLORS FOR TAXAPLOT #####

# 1. find out how many colors you need
numcol <- length(unique(alldata$plotnames))

# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(15)

# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")

# 4. Extract hex colors
hex = as.data.frame(newpal$hex)
colnames(hex) <- c("taxa_color")

# 5. Get list of taxa
tops = as.data.frame(c(unique(alldata$plotnames)))
colnames(tops) <- c("plotnames")

# 6. Join color list and taxa names
topcolors = cbind(tops, hex)

# 7. for the "others" plot name, replace that with grey 90 (this is just an astetic thing)
topcolors[topcolors$plotnames == "Others",]$taxa_color <- "grey90"

# 8. Make an object R can pull form for the colors
plotcolors <- topcolors$taxa_color
names(plotcolors) <- topcolors$plotnames

##### ORDER THE TAXA SO OTHERS ARE AT THE BOTTOM #####
## order by decreasing relative abundance
alldata = alldata[order(-alldata$relativeabundance),]

## get list of factors in order
natural.genus.order = as.list(c(unique(alldata$plotnames)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")

## set plot_names levels
plot.order = unlist(plot.order)

## order dataframe by relative abundance
alldata$plotnames = factor(alldata$plotnames, levels=c(plot.order))


##### MAKE THE TAXAPLOTS ######
field.groups = c(unique(field@sam_data$updated_sample_type))

## field taxaplot loop
for (i in field.groups){
  for (j in i){
    sub.df = subset(alldata, alldata$updated_sample_type== c(j) & lab_or_field_sampled=="field")
    
    myplot=ggplot(sub.df, aes(x=as.character(Row.names), y=relativeabundance, fill=plotnames))+
      geom_bar(stat = "identity")+
      scale_fill_manual(values=plotcolors)+
      facet_nested(.~updated_sample_type+site_id+as.numeric(month), scales="free")+
      labs(y="Relative Abundance", x="Sample", fill="Taxa")+
      labs(x=" ", y="Genus relative abundance in samples", fill="Order; Genus")+
      guides(fill=guide_legend(ncol=2))
    
    myplot
    ## save plot
    ggsave(myplot, filename=paste(j,"_field_taxaplot",".pdf",sep=""), width=20, heigh=7, 
           path=path)
  }
}


## lab taxaplot loop 
lab.groups = c(unique(lab@sam_data$updated_sample_type))

for (i in lab.groups){
  for (j in i){
    
    sub.df = subset(alldata, alldata$updated_sample_type== c(j) & paper_aquarium_day %in% c("day1", "day4"))
    
    myplot=ggplot(sub.df, aes(x=as.character(Row.names), y=relativeabundance, fill=plotnames))+
      geom_bar(stat = "identity")+
      scale_fill_manual(values=plotcolors)+
      facet_nested(.~updated_sample_type+aquarium_salinity+paper_aquarium_day, scales="free")+
      labs(y="Relative Abundance", x="Sample", fill="Taxa")+
      labs(x=" ", y="Genus relative abundance in samples", fill="Order; Genus")+
      guides(fill=guide_legend(ncol=2))
    
    myplot
    ## save plot
    ggsave(myplot, filename=paste(j,"_lab_taxaplot",".pdf",sep=""), width=20, heigh=7, 
           path=path)
  }
}

###### LAB MERISTEM TAXA BUBBLE PLOTS ######
## keep taxa in any top 15 group with lab meristems
labmeris = subset(top.df, top.df$Group %in% c("30meristemday4","30meristemday1",
                                              "20meristemday1","20meristemday4",
                                              "10meristemday1","10meristemday4"))
labtaxa = c(unique(labmeris$plotnames))

full = merge_phyloseq(field, lab)
full = tax_glom(full, taxrank = "Genus")
full@sam_data$rd_fult = sample_sums(full)
full = dephyloseq(full)

## only keep top meristem taxa 
full$plotnames = paste0(full$Order, ";", full$Genus)
full = subset(full, full$plotnames %in% c(labtaxa))

## caluclate relative abundace
full$ra = as.numeric(full$asv_abundance)/as.numeric(full$rd_fult)

labdf = subset(full, full$lab_or_field_sampled=="lab")
fielddf = subset(full, full$lab_or_field_sampled=="field")

##### SUMMARIZE TO MAKE BUBBLE PLOTS ######
labdf = ddply(labdf, c("paper_aquarium_day", "aquarium_salinity", "updated_sample_type", "plotnames", "Order"),
              summarise,
              meanra = mean(ra)) 

fielddf = ddply(fielddf, c("updated_sample_type", "month", "site_id", "plotnames", "Order"),
              summarise,
              meanra = mean(ra))

####### BUBBLE PLOTS #####

## meristem
labdf$presabs = ifelse(labdf$meanra >0, "pres", "abs")
labmeri = subset(labdf, labdf$updated_sample_type=="meristem")
labmeri$samp.id = paste0("lab",rownames(labmeri))

fielddf$presabs = ifelse(fielddf$meanra >0, "pres", "abs")
fieldmeri = subset(fielddf, fielddf$updated_sample_type=="meristem")
fieldmeri$samp.id = paste0("field",rownames(fieldmeri))

## merge for plot
meri = full_join(labmeri, fieldmeri)
meri$plot.x.axis = paste0(meri$aquarium_salinity, meri$site_id)

ggplot(meri, aes(x=plot.x.axis, y=plotnames, size=meanra, 
                 fill=plotnames, alpha=presabs))+
  geom_point(pch=21)+
  scale_alpha_manual(values=c(0.1, 1))+
  scale_fill_manual(values=plotcolors, guide = "none")+
  facet_nested(Order~updated_sample_type+paper_aquarium_day+month, scales="free", space="free")+
  theme(strip.text.y = element_text(angle = 0))

ggsave(filename="meristem_bubble_taxaplot.pdf", width=20, heigh=9.5, 
       path=path)
