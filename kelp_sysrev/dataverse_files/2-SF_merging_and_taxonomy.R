## https://github.com/benjjneb/dada2/issues/95
## with lots of help from Bia Segovia and Geneveive Lajoie

##### WORKSPACE SET UP #####
## load packages
library(tidyverse)
library(phyloseq)
library(dada2)
library(data.table)
library(zoo)


## set working directory 
setwd("/parfreylab/schenk/spring_freshet_output")


##### READ IN DATA #####
## load processed sequence files 
st2021 = readRDS("seqtab_nochim_SF_2021.rds")
st2022.lib1 = readRDS("seqtab_nochim_seqtab_SF_2022_lib1.RDS")
st2022.lib2 = readRDS("seqtab_nochim_seqtab_SF_2022_lib2.RDS")


##### MERGE DATASETS ####
## merge both datatsets 
st.all <- mergeSequenceTables(st2021, st2022.lib1, st2022.lib2, tryRC=TRUE)
## remove bimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE,verbose = FALSE)
dim(seqtab.nochim)


# Collapse ASVs that only differ by end base pair
seqtab.nochim<-collapseNoMismatch(seqtab.nochim, verbose=T) # Check Issue 716 from dada2 github
# https://github.com/benjjneb/dada2/issues/626
dim(seqtab.nochim)

## save seqtab file as RDS
write_rds(seqtab.nochim, "seqtab_nochim_SF_2021_and_2022.rds")

##### ASSIGN TAXONOMY #####
## SILVA needs to be downloaded in the directory to run
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v138_train_set.fa", 
                       multithread=TRUE, tryRC=TRUE)

# Save taxonomy table
write.table(data.frame("row_names"=rownames(taxa),taxa),"SF_taxonomy_SILVAv138_16s.txt", 
            row.names=FALSE, quote=F, sep="\t")

## add species information to taxonomy data
taxa_sp <- addSpecies(taxa, "silva_species_assignment_v138.fa")

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa_sp[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa_sp[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa_sp[NAs,1] <- "Unassigned" #apply new label to identified indices
colnames(taxa_sp) <- c("domain", "phylum", "class", "order", "family", "genus", "species") 


# Save taxonomy table
write.table(data.frame("row_names"=rownames(taxa_sp),taxa_sp),"SF_taxonomy_noNAs_SILVAv138_16s.txt", 
            row.names=FALSE, quote=F, sep="\t")


##### READ IN FILES #####
## read in sequence table
seqtab.nochim = readRDS("seqtab_nochim_SF_2021_and_2022.rds")
## read in taxonomy table
tax = read.table(file="SF_taxonomy_noNAs_SILVAv138_16s.txt", sep='\t', header = T)
## read in metadata
metadata = read.csv("sf_metadata_sequencing.csv")


##### PROPAGATE THE TAXONOMY INTO EMPTY LOWER RANKS #####
# propagates taxonomy from left
tax_propagated = tax %>%
  t() %>% #transpose (moves taxonomy from column names to row names)
  na.locf() %>% #fill the NAs with the values from the cell to the left (higher taxonomic rank)
  t() # transpose back to have column names be taxonomy and row names be ASV


# make into dataframe
tax_propagated = as.data.frame(tax_propagated)

# extract rownames into column
setDT(tax_propagated, keep.rownames = TRUE)[] 

## add ASV ID marker
tax_propagated$asv_id = paste0("ASV", tax_propagated$rn)
tax_propagated = tax_propagated[,-c(1)]

# change column names
colnames(tax_propagated) <- c("asv_sequence","Kingdom","Phylum","Class","Order","Family","Genus", "Species","asv_id") 

## re-order taxonomy file
tax_propagated = tax_propagated[,c(2:9,1)]

##### MATCH OTU SEQUENCES WITH ASV ID #####
## turn seqtab dataset and get sequences out
seqtab.nochim = as.data.frame(t(seqtab.nochim)) 

## get out sequence from rownmaes
seqtab.nochim = seqtab.nochim %>% rownames_to_column(var="asv_sequence")

## merge taxonomy and otu table
taxotu = full_join(tax_propagated, seqtab.nochim)
ncols = ncol(taxotu)

## separate out taxonomy and otu table but keep asv id in both
tax_propagated = taxotu[,c(1:9)]
seqtab.nochim = taxotu[,c(8,10:1104)]

## make asv id rownames
tax_propagated = tax_propagated %>% column_to_rownames(var="asv_id")
seqtab.nochim = seqtab.nochim %>% column_to_rownames(var="asv_id")

##### Get LIST OF COLUMN NAMES OF SEQTAB TO MATCH TO METADATA #####
seqnames = as.data.frame(colnames(seqtab.nochim))

## rename seqnames column
names(seqnames)[1] = "illumina_sample_id"

## merge seqmanes and metadata
metadata = full_join(seqnames, metadata)



##### CREATE PHYLOSEQ OBJECT #####
##
otu_mat = as.matrix(seqtab.nochim)

tax_mat = as.matrix(tax_propagated)

metadata = metadata %>% column_to_rownames(var="illumina_sample_id")

## transform into phyloseq-ready objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)  
samples = sample_data(metadata)

## make the phyloseq object
SF = phyloseq(OTU, TAX, samples)

## check the phyloseq object is formatted correctly
SF
rank_names(SF)
sample_variables(SF)

#phyltax = tax_table(SF)
#print(phyltax)

#head(SF@otu_table)


#####  ADD INITIAL SAMPLE SUMS BEFORE FILTERING ####
SF@sam_data$unfiltered_read_depth = sample_sums(SF)



## SAVE PHYLOSEQ OBJECT
write_rds(SF, "unfiltered_phyloseq_SF.rds")

## GO TO 3-SF_filtering_and_coveraged_based_rarefaction.R


