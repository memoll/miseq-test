## Data preparation for analysis with DADA2
## Mona Parizadeh - March 2018

# As mentioned in DADA2 tutorial, for the analysis we'll need demultiplexedpaired-end sequences in fastq format that have been split (or “demultiplexed”) by sample and also the barcodes/adapters have already been removed.

## 1. Trim adapters - bbduk ####
#/data/apps/bbmap/bbduk.sh –Xmx1g in=file_R#_001.fastq out=adapter-trimmed/file-trimmed_R#.fastq ref=/data/apps/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe


## 2. Make mapping_file.tsv - excel #### 
# #SampleID - BarcodeSequence	- LinkerPrimerSequence -	Barcode	- ReversePrimer	- Description

# Problem: in our Miseq SOP we're using barcodes w/ various lengths
# Strategy: making fake barcodes, by adding a part of the primer to the end of each barcode

# 3. Write the mapping_file.txt - R ####
mapping_file <- read.table("miseq-test-tables/miseq-run-test-2016 - mapping-file.tsv")
colnames(mapping_file) <- (c("#SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Barcode", "ReversePrimer", "Description"))
# Add primers to barcodes and cut them to the length of the longest barcode
mapping_file$BarcodeSequence <- strtrim(paste0(mapping_file$BarcodeSequence,mapping_file$LinkerPrimerSequence),max(nchar(as.character(mapping_file$BarcodeSequence)))); nchar(mapping_file$BarcodeSequence) # for BarcodeSequence (F)
mapping_file$Barcode <- strtrim(paste0(mapping_file$Barcode,mapping_file$Reverse), max(nchar(as.character(mapping_file$Barcode)))); nchar(mapping_file$Barcode) # for Barcode (R)
# Put new barcodes together and replace all the characters of the barcodes other than ATCG with "."
mapping_file$BarcodeSequence <- gsub("[^ATCG]",".", paste0(mapping_file$BarcodeSequence,mapping_file$Barcode)) 
mapping_file$Barcode <- NULL # clear the barcode (R) column
#mapping_file
write.table(mapping_file, "miseq-test-tables/mapping-miseq-run-test-2016.txt", sep = "\t", row.names = F, quote = F)

# 4. Validate the mapping file - qiime ####
# validate_mapping_file.py -m miseq-test-tables/mapping-miseq-run-test-2016.txt -o validate-mapping

# 5. Make mapping file - idemp (barcode tab identifier) ####
idemp_mapping_file <- mapping_file[,c(2,1)]; colnames(idemp_mapping_file) <- NULL
write.table(idemp_mapping_file, "miseq-test-tables/idemp-mapping-file.txt", sep = "\t", row.names = F, quote = F)

# 6. Extract barcodes - qiime ####
# Based of out sequencing protocol, R1 is reverse and R2 is forward; so for this step bc1 = R -> --bc1_len = 17 and bc2 = F -> --bc2_len = 14
# extract_barcodes.py -f rawdata/adapter-trimmed/Mona-run-test_S1_L001-trimmed_R2.fastq -r rawdata/adapter-trimmed/Mona-run-test_S1_L001-trimmed_R1.fastq -c barcode_paired_end -m miseq-test-tables/mapping-miseq-run-test-2016.txt --attempt_read_reorientation --bc1_len 17 --bc2_len 14 -o dada2-analysis/processed_seqs

# 7. Trim the rest of the primers from reads1 and reads2 - :####

# @param id : identifier
# @param seq.read : character vector containing reads
# @param seq.qual : character vector containing quality
# @return
fastq_file <- read_file("reads1.1000.fastq")
cut_left <- function(fastq_file){
  seq.read <- subseq(fastq_file ((1:length(fastq_file) %% 4 + 1) %in% c(2,4)), start=20, end=nchar(fastq_file))
}


library(Biostrings)
fastq_file <- readDNAStringSet("reads1.1000.fastq", format("fastq"))
fastq_file_19 <- subseq(fastq_file, start=20, end=nchar(fastq_file))
fastq_file_19[246]

write(fastq_file_19, "", sep = "\n")


#install.packages("readr")
library(readr)
#library(tools)
fastq_file <- "reads1.fastq"
fastq_data <- read_lines(fastq_file, skip = 0); length(fastq_file)
sel_line <- (1:length(fastq_data) %% 4 + 1) %in% c(2,4) # selects the lines that we need to manipulate, 2: reads, 4: quality.
trim_size <- max(nchar(as.character(mapping_file$LinkerPrimerSequence)))
prm_trimmed <- ifelse(sel_line, substring(fastq_data[10], trim_size + 1, nchar(fastq_data)), fastq_data) # removes the beginning of each read and quality line and keep fromn the max length of the primer to the end of the read and quality lines
write_lines(prm_trimmed, paste0(file_path_sans_ext(fastq_file)),".prm-trimmed.fastq")


# sp = split(x, line_nb)
# read <- substring(sp$'2', 31, nchar(sp$'2'))
# qual <- substring(sp$'4', 31, nchar(sp$'4'))
#substring(sp$`2`,30,length(sp$`2`))
#read <- sapply(strsplit(sp$'2',substring(sp$'2', 1, 30) ,2), `[`, 2)
#qual <- sapply(strsplit(sp$'4',substring(sp$'4', 1, 30) ,2), `[`, 2)
#table(read==read2)



# 8. Demultiplex - idemp:####
#idemp -b idemp-mapping-file.txt -I1 processed_seqs/barcodes.fastq -R1 processed_seqs/reads2.fastq -R2 processed_seqs/reads1.fastq -m 1 -o idemp-demultiplexed

# 12. Verify the number of demultiplexed samples files (Number of files in idemp-demultiplexed - 2)/2 :
# Number of files: cd idemp-demultiplexed > ls -1 | wc -l     










