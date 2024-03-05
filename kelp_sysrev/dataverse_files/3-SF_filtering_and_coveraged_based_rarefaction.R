##### WORKSPACE SETUP #####
library(phyloseq)
library(tidyverse)
library(dplyr)
library(data.table)
library(readr)
library(devtools)
library(metagMisc)
library(iNEXT)
library(vegan)

## tell R where to get the files

## read in files
sf = readRDS("unfiltered_phyloseq_SF.rds")
sf

##### REMOVE UNWANTED TAXA FROM TAXONOMY FILE #####
sf = subset_taxa(sf, Kingdom!="Eukaryota"&
                   Kingdom!="Unassigned" &
                   Order!="Chloroplast" &
                   Family!="Mitochondria" &
                   Family!="Chloroplast" &
                   Species !="Chloroplast" &
                   Genus!="Pseudomonas")


##### FILTERING - REMOVE SAMPLES WITH LESS THAN 1000 READS #####
## Remove samples with less than N reads 
sample_sums(sf)
plot(sort(sample_sums(sf)))

## add number of total sequences in a sample (Read depth)
# e.g. if one sample had 3 ASV with 1, 5, and 7 reads respectively, the Read_depth of that sample would be 13
# this is the same as colSums for this dataset because of the orientation of the OTU table (rows=asv, columns = samples)
# rowSums would sum the total number of sequences each ASV has. So if across 10 samples that ASV was detected 50 time, the rowSum would be 50
sf@sam_data$read_depth_noofftargets = sample_sums(sf) 

## check which samples have less than 1000 reads
which(sf@sam_data$read_depth_noofftargets < 1000) 
# this one gets rid of all the controls and nothing else
# decided to remove samples with less than 1000 reads
sf.pruned <- prune_samples(sample_sums(sf) >= 1000, sf)
## write file to know which reds were lost here
sf.below1000 <- prune_samples(sample_sums(sf) < 1000, sf)
sf.below1000 = as.matrix(sf.below1000@sam_data)
write_rds(sf.below1000, "SF_samples_less_than_1000.csv")

##### FILTERING - REMOVE INDIVIDUAL ASVS WITH LESS THAN 100 READS #####
## extract OTU dataframe from phyloseq object
otu.pruned <- as.data.frame(as.matrix(otu_table(sf.pruned)))

## remove OTU (rows) with less than 100 reads accross whole dataset but keep all samples
## make sure asv sequence is rownames and sample id is column name
otu.pruned$rowsum = rowSums(otu.pruned)

## remove low frequency asvs (less than 0.001%)
total_asvs = sum(otu.pruned$rowsum)
otu.pruned$total = total_asvs
otu.pruned$percent_abundance = otu.pruned$rowsum/otu.pruned$total
otu.pruned = subset(otu.pruned, otu.pruned$percent_abundance>0.00001)

## remove rowsum column
otu = subset(otu.pruned, select=-c(rowsum, total, percent_abundance))


##### MAKE ALL THE CELLS IN THE OTU TABLE WITH VALUES 5 OR LESS 0 #####
otu <- mutate_all(otu, funs(ifelse(. < 5, 0, .)))


##### REMOVE ASVs FOUND IN 2 SAMPLES OR LESS #####

# has sample ID as column name and asv is as row name. Needs to be this way to use richness function 
# function to calculate richness, sums along a row (OTU)
richness = function(x){return(sum(x>0))}

## calculate richness on entire dataframe
otu$richness = apply(otu,1,richness) # use all columns of otu dataframe
summary(otu$richness)

## remove OTU (rows) with richness 2 or less (found in two samples or less) but keep all samples (columns)
otu = subset(otu, otu$richness>2)
## check that it worked
summary(otu$richness)
## remove richness column
otu = subset(otu, select=-c(richness))


##### CREATE AND READ BACK IN FILTERED BUT NOT RAREFIED PHYLOSEQ OBJECT ######

## format for phyloseq
otu_mat = as.matrix(otu)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat = sf.pruned@tax_table
TAX = tax_table(tax_mat)  

## get out dataframe
samples = as.data.frame(as.matrix(sf.pruned@sam_data))

## get metadata ready for phyloseq
samples = sample_data(samples)

# make new phyloseq object
sf.prerarefaction = phyloseq(OTU, TAX, samples)
# check that the phyloseq object was made correctly
sf.prerarefaction
rank_names(sf.prerarefaction)
sample_variables(sf.prerarefaction)
#phyltax = tax_table(sf.prerarefaction)
#print(phyltax)
#head(otu_table(sf.prerarefaction))

##### REMOVE SAMPLES WITH DESMERESTIA - LAB CULTURE ROUND 5, AQUARIA 10 C #####
sf.prerarefaction = subset_samples(sf.prerarefaction, desmerestia!="y")

sf.prerarefaction

##### FINAL DENOISING AND SAVE FILTERED DATA #####

## remove taxa with empty otus (from denoising)
sf.prerarefaction = prune_taxa(taxa_sums(sf.prerarefaction)>0, sf.prerarefaction)

## get final read depth
sf.prerarefaction@sam_data$read_depth_filtered = sample_sums(sf.prerarefaction)

## remove samples with highest sample counts
histogram(sf.prerarefaction@sam_data$read_depth_filtered, breaks=100)
sf.prerarefaction = subset_samples(sf.prerarefaction, read_depth_filtered<200000)
histogram(sf.prerarefaction@sam_data$read_depth_filtered, breaks=100)

write_rds(sf.prerarefaction, "SF_2021_and_2022_filtered001percent_notrarefied.rds")


#sf_otu = as.data.frame(sf@otu_table)
#write.csv(sf_otu, "filtered_otu_table_notrarefied.csv")


##### 6. look at minimum, mean, and maximum sample counts, if desired #####
## read in data
sf_prerare = readRDS("SF_2021_and_2022_filtered001percent_notrarefied.rds")

smin <- min(sample_sums(sf_prerare)) #742
meanreads <-mean(sample_sums(sf_prerare)) #25,997.17
smax <- max(sample_sums(sf_prerare)) #482,900
totalreads <- sum(sample_sums(sf_prerare)) #24,073,379

get_sample(sf_prerare)
sample_sums(sf_prerare)


##### normalizing data using coverage based rarefaction ASV LEVEL #####

### phyloseq_coverage_raref(physeq, coverage = NULL, iter = 1, replace = F, correct_singletons = FALSE, seeds = NULL, multithread = F, drop_lowcoverage = F, ...)

### Samples standardized by size will have different degrees of completeness. When we compare samples with the same coverage, we are making sure that samples are equally complete and that the unsampled species constitute the same proportion of the total individuals in each community (Chao, Jost, 2012).
### i.e. a seagrass sample with 10,000 total reads will have a different coverage than a seawater sample with 10,000 reads if seawater samples have many more species


### check if taxa are rows in the phyloseq object
taxa_are_rows(sf_prerare)
## mine is TRUE

### prepare otu table for function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(sf_prerare)),
  correct_singletons = T)

#View(as.data.frame(otu_table(laby_prerare)))
### check if read counts are correct (samples should show "numeric" in the second column)

SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)


### Due to the stochasticity introduced in random subsampling results could be slightly different.
### So you have to average diversity estimates or sample dissimilarities across multiple rarefactions.

### run coverage-based rarefaction (Chao & Jost, 2012) correcting for singletons (Chiu & Chao 2016)
### 1,000 iterations
all_16S_COVERAGE_RAREF_1000 <- phyloseq_coverage_raref(physeq=sf_prerare, coverage = 0.8, iter = 1000, replace = F, correct_singletons = TRUE, drop_lowcoverage = F) 
saveRDS(all_16S_COVERAGE_RAREF_1000, "SF_COVERAGE_RAREF_1000.rds")

all_16S_COVERAGE_RAREF_1000=readRDS("SF_COVERAGE_RAREF_1000.rds")
### Average otu tables from all iterations to get a final robust table ###
subset_phylo_objects_1000 <- all_16S_COVERAGE_RAREF_1000[c(1:1000)]

### first, extract otu tables from phyloseq objects
### this is how you do it for a single phyloseq object:
### y <- as.data.frame(t(phyloseq::otu_table(all_16S_COVERAGE_RAREF)))
### now do it for the list of phyloseq objects
otu_tables_1000 <- lapply(subset_phylo_objects_1000, function(z) as.data.frame(t(phyloseq::otu_table(z))))

### average all matrices to get the mean abundance across all iterations
average_otu_tables_1000 <- Reduce("+",otu_tables_1000)/length(otu_tables_1000)
### IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round <- average_otu_tables_1000 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))


### add SampleID column back
average_otu_tables_1000_round$sample_id<- rownames(average_otu_tables_1000) 
average_otu_tables_1000_round <- average_otu_tables_1000_round %>% 
  select(sample_id, everything())
write.csv(average_otu_tables_1000_round, "SF_average_otu_tables_rarefied001percent.csv", quote=F, row.names=F )

##### save as final phyloseq object #####
## metadata
metadata = as.matrix(sf_prerare@sam_data)
metadata = as.data.frame((metadata))
samples = sample_data(metadata)

## TRANSPOSED asv table 
otu_1000 = average_otu_tables_1000_round[,-1]
# rename rows
rownames(otu_1000) <- average_otu_tables_1000_round[,1]

## as matrix an transpose
otu_matrix = t(as.matrix(otu_1000))

# format for phyloseq
OTU = otu_table(otu_matrix, taxa_are_rows = TRUE)

## taxonomy
TAX = sf_prerare@tax_table

## make the phyloseq object
sf_rare = phyloseq(OTU, TAX, samples)
sf_rare
rank_names(sf_rare)
sample_variables(sf_rare)

sf_rare@sam_data$read_depth_rarefied = sample_sums(sf_rare)

## save phyloseq object as RDS
write_rds(sf_rare, "SF_2021_and_2022_coveraged_rarified001percent.RDS")






