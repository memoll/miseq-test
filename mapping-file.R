## Mona Parizadeh - March 2018
# How to prepare mapping file?

# 1. In excel: #SampleID - BarcodeSequence	- LinkerPrimerSequence -	Barcode	- ReversePrimer	- Description; save as .tsv ####

# 2. Import mapping_file.tsv in R:####
mapping_file <- read.table("miseq-test-tables/miseq-run-test-2016 - mapping-file.tsv")
colnames(mapping_file) <- (c("#SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Barcode", "ReversePrimer", "Description"))
#mapping_file

# Verify the length of the longest barcode - we'll need it to make fake barcodes later
max_BC_seq <- max(nchar(as.character(mapping_file$BarcodeSequence)))
max_BC <- max(nchar(as.character(mapping_file$Barcode)))

# 3. Add primers to the end of each barcode sequence - in order to make fake barcodes with the same length later ####
#BC_seq_prm <-  paste0(mapping_file$BarcodeSequence,mapping_file$LinkerPrimerSequence)
#mapping_file$BarcodeSequence <- BC_seq_prm
mapping_file$BarcodeSequence <-  paste0(mapping_file$BarcodeSequence,mapping_file$LinkerPrimerSequence)
#BC_prm <- paste0(mapping_file$Barcode,mapping_file$Reverse)
#mapping_file$Barcode <- BC_prm
mapping_file$Barcode <- paste0(mapping_file$Barcode,mapping_file$Reverse)

# 4. Cut all the new barcodes to the same length -length of the longest barcode-primer ####

# BarcodeSequence column:
#nchar(mapping_file$BarcodeSequence)
BC_seq_trimmed <- strtrim(mapping_file$BarcodeSequence,max_BarcodeSequence)
nchar(BarcodeSequence_trimmed)
mapping_file$BarcodeSequence <- BarcodeSequence_trimmed

# Barcode column:
#nchar(mapping_file$Barcode)
Barcode_trimmed <- strtrim(mapping_file$Barcode,max_Barcode)
nchar(Barcode_trimmed)
mapping_file$Barcode <- Barcode_trimmed

# 5. Add Barcode (reverse/r2) to BarcodeSequence (forward/r2) ####
BarcodeSequence_r1_r2 <- paste0(mapping_file$BarcodeSequence,mapping_file$Barcode)

# If all the characters of the barcodes other than ATCG with "."
BarcodeSequence_r1_r2 <- gsub("[^ATCG]",".",BarcodeSequence_r1_r2)
mapping_file$BarcodeSequence <- BarcodeSequence_r1_r2

# Remove the Barcode (reverse/r2) column
mapping_file$Barcode <- NULL
mapping_file

# 6. Export the mapping_file as a .txt file ####
write.table(mapping_file, "R_projects/miseq-test/mapping-miseq-run-test-2016.txt", sep = "\t", row.names = F, quote = F)

# 7. Validate the mapping file - qiime:####
#bash: validate_mapping_file.py -m mapping-miseq-run-test-2016.txt -o validate-mapping

# check if there are reverses of each barcodes of F and R in R1 and R2, and for the paired barcodes as well (reversed and not reversed)
# In R: reverse("ACAGCCACCCATCGAAAGAGCAACATCCTAA") and grep it in R1 and R2 - bash

# 8. Make mapping file for idemp - barcode tab identifier:####
idemp_mapping_file <- mapping_file[,c(2,1)]
colnames(idemp_mapping_file) <- NULL
write.table(idemp_mapping_file, "R_projects/miseq-test/idemp-mapping-file.txt", sep = "\t", row.names = F, quote = F)

# 9. Extract barcodes - qiime:####
# bash: extract_barcodes.py -f rawdata/adapter-trimmed/Mona-run-test_S1_L001-trimmed_R2.fastq -r rawdata/adapter-trimmed/Mona-run-test_S1_L001-trimmed_R1.fastq -c barcode_paired_end -m mapping-miseq-run-test-2016.txt --attempt_read_reorientation --bc1_len 17 --bc2_len 14 -o processed_seqs
# Based of out sequencing protocol, R1 is reverse and R2 is forward; so for this step bc1 = R -> --bc1_len = 17 and bc2 = F -> --bc2_len = 14

# 10. Trim the rest of primers from reads1 and reads2 - :####
# bash: 
# /data/apps/bbmap/bbduk.sh in=processed_seqs/reads1.fastq out=processed_seqs/primer-trimmed/reads1.primerTrimmed.fastq ftl=36
# reads1 : 18740998 - 18717004 = 23994
# ftl = 30 -> 18717239
# /data/apps/bbmap/bbduk.sh in=processed_seqs/reads2.fastq out=processed_seqs/primer-trimmed/reads2.primerTrimmed.fastq ftl=30
# reads2 : 18740998 - 18717327 = 23671
#ftl=36: 18717176


#install.packages("readr")
library(readr)

fastq_file = "/data/users/mona/miseq-test-2016/dada2-analysis/processed_seqs/reads1.fastq"
x = read_lines(fastq_file, skip = 0); length(x)
sel_line = (1:length(x) %% 4 + 1) %in% c(2,4) #selects the lines that we need to manipulate, 2: reada, 4: quality.
x2 = ifelse(sel_line , substring(x, 37, nchar(x)), x)
write_lines(x2,paste0(fastq_file,".trim30"))

# sp = split(x, line_nb)
# read <- substring(sp$'2', 31, nchar(sp$'2'))
# qual <- substring(sp$'4', 31, nchar(sp$'4'))
#substring(sp$`2`,30,length(sp$`2`))
#read <- sapply(strsplit(sp$'2',substring(sp$'2', 1, 30) ,2), `[`, 2)
#qual <- sapply(strsplit(sp$'4',substring(sp$'4', 1, 30) ,2), `[`, 2)
#table(read==read2)



# 11. Demultiplex - idemp:####
#idemp -b idemp-mapping-file.txt -I1 processed_seqs/barcodes.fastq -R1 processed_seqs/reads2.fastq -R2 processed_seqs/reads1.fastq -m 1 -o idemp-demultiplexed

# 12. Verify the number of demultiplexed samples files (Number of files in idemp-demultiplexed - 2)/2 :
# Number of files: cd idemp-demultiplexed > ls -1 | wc -l     










