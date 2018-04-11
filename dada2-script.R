# Script to prepare sequences for analysis with DADA2
# Based on DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial_1_4.html
# Data: Miseq-test
# Mona Parizadeh - March 2018

#### Package preparation ####
# Install dada2 and get ready
#source("https://bioconductor.org/biocLite.R")
#biocLite("dada2")
biocLite("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2")
library(dada2)
packageVersion("dada2")

#require(parallel)
#detectCores() #12


path <- "~/../../data/users/mona/miseq-test-2016/dada2-analysis/processed_seqs/idemp-demultiplexed/"
list.files(path)


#####should add r1 r2 to reads1 and reads2 (in idemp)

#### Filter and trim ####
# List of orward and reverse fastq filenames have format
fnFs <- sort(list.files(path, pattern="reads1.fastq_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="reads2.fastq_", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)

# Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:6]) # forward reads 
plotQualityProfile(fnRs[1:6]) # reverse reads 
# The reverse reads are of significantly worse quality, especially at the end, which is common in Illumina sequencing

# Perform filtering and trimming 
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter the forward and reverse reads
#detectCores() 12: in order to choose multithread, so that it will not use all the cores.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(180,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, verbose = TRUE, matchIDs = TRUE , multithread=10) 
# matchIDs: match the paired reads by their Illumina ID
# On Windows set multithread=FALSE
# maxN = 0; After truncation (various at the end), sequences with more than maxN Ns will be discarded. Note that dada does not allow Ns.
# truncLen: Cut where the forward and reverse reads quality drop off, based on the plots (For paired-end reads consider the length of your amplicon when choosing  truncLen as your reads must overlap after truncation in order to merge them later).
head(out)
write.table(out, "out.txt", sep = "\t", row.names = F, quote = F)

#### Calculate error rates ####
errF <- learnErrors(filtFs, multithread=10)  
#Convergence after  8  rounds.
#116611380  total bases in  647841  reads used for learning the error model.

errR <- learnErrors(filtRs, multithread=10) 
#Convergence after  5  rounds.
#142525020  total bases in  647841  reads used for learning the error model.

plotErrors(errF, nominalQ=TRUE)
#The red line shows the error rates expected under the nominal definition of the Q-value.

#### Dereplicate the filtered fastq files ####
#Dereplication combines all identical sequencing reads into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence, to reduce computation time by eliminating redundant comparisons
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#### Infer the sequence variants in each sample ####
dadaFs <- dada(derepFs, err=errF, multithread=10)
dadaRs <- dada(derepRs, err=errR, multithread=10)

dadaFs[[1]] # 15 sequence variants were inferred from 490 input unique sequences
dadaRs[[1]] # 16 sequence variants were inferred from 517 input unique sequences
dadaFs[[3]] # 588 sequence variants were inferred from 47660 input unique sequences
dadaRs[[3]] # 726 sequence variants were inferred from 43351 input unique sequences

#### Merge the denoised forward and reverse reads ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE) # removes paired reads that did not exactly overlapped
# verbose = TRUE: to print a summary of the function results
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
# If the samples 0 paired-reads after merging, it means reads don't overlap. That is, the primers were separated by too much to meet in the middle.

#### Construct sequence table ####
seqtab <- makeSequenceTable(mergers) #The sequences being tabled vary in length.
dim(seqtab) #53 73260
#seqtab is a matrix w/ rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. 
# The lengths of the merged sequences should all fall within the expected range for this V4 amplicon.
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#### Remove chimeric sequences ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=10, verbose=TRUE)
#A bimera is a two-parent chimera, in which the left side is made up of one parent sequence, and the right-side made up of a second parent sequence.
#Identified 46129 bimeras out of 73260 input sequences.
dim(seqtab.nochim) #253 27131
sum(seqtab.nochim)/sum(seqtab) #0.8274735
#Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though). If most of your reads were removed as chimeric, upstream processing may need to be revisited. In almost all cases this is caused by primer sequences with ambiguous nucleotides that were not removed prior to beginning the DADA2 pipeline.

#### Track reads through the pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

#### Assign taxonomy ####
# Download RDP taxonomic training data formatted for DADA2 (rdp_train_set_16.fa.gz) from https://zenodo.org/record/801828#.Wi31VbQ-cWo
# Check for other datasets formatted for DADA2: https://benjjneb.github.io/dada2/training.html

# Assign - RDP
taxaRDP <- assignTaxonomy(seqtab.nochim, "~/../../data/users/mona/miseq-test/rdp_train_set_16.fa.gz", multithread=10)
unname(head(taxaRDP))
#     [,6]              
#[1,] "Methylobacterium"
#[2,] "Methylobacterium"
#[3,] "Methylobacterium"
#[4,] "Quadrisphaera"   
#[5,] "Buchnera"        
#[6,] "Methylobacterium"
write.table(taxaRDP, "taxa-rdp.miseq-test.txt", sep = "\t", row.names = F, quote = F)

# Assign - Silva
taxaSilva <- assignTaxonomy(seqtab.nochim, "~/../../data/users/mona/miseq-test/silva_nr_v132_train_set.fa.gz", multithread=10)
unname(head(taxaSilva))
#     [,6]              
#[1,] "Methylobacterium"
#[2,] "Methylobacterium"
#[3,] "Methylobacterium"
#[4,] "Quadrisphaera"   
#[5,] "Buchnera"        
#[6,] "Methylobacterium"
write.table(taxaSilva, "taxa-silva.miseq-test.txt", sep = "\t", row.names = F, quote = F)

# Assign taxonomy at species level - silva
taxaSpeciesSilva <- addSpecies(taxa,"~/../../data/users/mona/miseq-test/silva_species_assignment_v128.fa.gz")
head(taxaSpeciesSilva)
#[1,] "Methylobacterium"
#[2,] "Methylobacterium"
#[3,] "Methylobacterium"
#[4,] "Quadrisphaera"   
#[5,] "Buchnera"        
#[6,] "Methylobacterium"
write.table(taxaSpeciesSilva, "taxa-species-silva.miseq-test.txt", sep = "\t", row.names = F, quote = F)

#### Evaluating DADA2’s accuracy on the Sy.S.D.MVM.06.16.6A1.PS2.E1 (for example) community ####
unqs.Sy.S.D.MVM.06.16.6A1.PS2.E1 <- seqtab.nochim["Sy.S.D.MVM.06.16.6A1.PS2.E1.fastq.gz", ]
unqs.Sy.S.D.MVM.06.16.6A1.PS2.E1 <- sort(unqs.Sy.S.D.MVM.06.16.6A1.PS2.E1[unqs.Sy.S.D.MVM.06.16.6A1.PS2.E1>0], decreasing=TRUE)
cat("DADA2 inferred", length(unqs.Sy.S.D.MVM.06.16.6A1.PS2.R1), "sample sequences present in the Sy.S.D.MVM.06.16.6A1.PS2.E1 community.\n")
#DADA2 inferred 793 sample sequences present in the Sy.S.D.MVM.06.16.6A1.PS2.E1 community.

#### Phyloseq ####
biocLite('phyloseq')
library(phyloseq)
packageVersion("phyloseq") #‘1.20.0’

#install.packages("ggplot2")
library(ggplot2)
packageVersion("ggplot2") #‘2.2.1’

# Make a data.frame holding the sample data
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "y."), `[`, 2) # split names on the letter y. and keep the second part
subject <- sapply(strsplit(subject, ".fastq.gz"), `[`, 1)

# Host
host <- substr(subject,1,1) # (x, start, stop)
host$"S" <- "Soil"
host$"L" <- "Phyllosphere"
length(host$"S")

# Treatment
treatment <- substr(subject,17,99)
#####neonic <- substr((treatment %% 2 == 0), 1, 1)

day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out