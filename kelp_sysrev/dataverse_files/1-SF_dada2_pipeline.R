##### WORKSPACE SET UP #####
## load packages
library(dada2)
library(readr)
library(stringr)
library(tidyverse)

## set path on cluster where your reads are
path = ""


## list the files in the path
list.files(path)

## set sample names
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) #change the pattern to match all your R1 files
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1) 
sample.names <- sapply(strsplit(basename(fnRs), "_S"), `[`, 1) 
#change the delimiter in quotes and the number at the end of this command to decide how to split up the file name, 
#and which element to extract for a unique sample name

head(sample.names)


##### IDENFITY PRIMER SEQUENCE AND LENGTH #####
FWD = nchar(c("GTGYCAGCMGCCGCGGTAA"))
REV = nchar(c("GGACTACHVGGGTWTCTAAT"))

##### INSPECT READ QUALITY PROFILES #####
#plotQualityProfile(fnFs[1:9])
#plotQualityProfile(fnFs[60:63])
#plotQualityProfile(fnRs[1:9])
#plotQualityProfile(fnRs[60:63])


##### FILTER AND TRIM #####
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R2_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

##### SANITY CHECK AND FILTERING #####

## last sanity check that sample match 
any(duplicated(c(fnFs, fnRs)))
any(duplicated(c(filtFs, filtRs)))

length(fnFs)
length(fnRs)

head(fnFs)
head(fnRs)
head(filtFs)
head(filtRs)


## filter out basic things, like N in sequencing and set minimal length for sequence
filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
              truncLen=c(FWnumber, RVnumber), ## cut the sequences at this length. FW will usually be longer than RV
              maxN=0, ## must be 0. DADA2 can't have N in sequences 
              maxEE=c(5,5), ## max error. Stay below 5
              truncQ=2, 
              rm.phix=TRUE,
              compress=TRUE, 
              multithread=TRUE,
              matchIDs=TRUE,
              trimLeft = c(19, 20)) ## trim off primer length



##### LEARN THE ERROR RATES #####

## make sure all files you are asking for actually exist
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]


## learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

## plot error rates
#plotErrors(errF, nominalQ=TRUE)


##### SAMPLE INFERENCE #####
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]
# dada-class: object describing DADA2 denoising results
# 6 sequence variants were inferred from 28 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16


##### MERGE PAIRED READS #####
## aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
## inspect the merger data.frame from the first sample
head(mergers[[1]])


##### CONSTRUCT SEQUENCING TABLE #####
## construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

## inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##### REMOVE CHIMERAS #####
## The core dada method corrects substitution and indel errors, but chimeras remain. Remove them here 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 650 bimeras out of 7781 input sequences.
dim(seqtab.nochim)
## outputs fraction of sequence that are NOT chimeras
sum(seqtab.nochim)/sum(seqtab)


##### REMOVE SEQUENCES THAT ARE TOO SHORT/LONG AND SAVE PROCESSED SEQUENCE TABLE #####
## inspect distribution of sequence lengths
## target V4 region from 515f and 806r primers is ~253 in length
table(nchar(getSequences(seqtab.nochim))) 
## see plot of sequence lengths 
#plot(table(nchar(getSequences(seqtab.nochim))))



## save files as RDS
setwd("")
write_rds(seqtab.nochim, "")


## go to 2-SF_merging_and_taxonomy.R for next steps 