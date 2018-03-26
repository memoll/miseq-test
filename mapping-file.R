## Mona Parizadeh - March 2018

# How to prepare mapping file?

# 1. In excel: #SampleID - BarcodeSequence	- LinkerPrimerSequence -	Barcode	- ReversePrimer	- Description; save as .tsv ####

# 2. Import mapping_file.tsv in R:####
mapping_file <- read.table("Google Drive/mapping-miseq-run-test-2016.tsv")
colnames(mapping_file) <- (c("#SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Barcode", "ReversePrimer", "Description"))
mapping_file

# Verify the length of the longest barcode - we'll need it to make fake barcodes later
max_BarcodeSequence <- max(nchar(as.character(mapping_file$BarcodeSequence)))
max_Barcode <- max(nchar(as.character(mapping_file$Barcode)))

# 3. Add primers to the end of each barcode sequence - in order to make fake barcodes with the same length later ####
BarcodeSequence_primer <-  paste0(mapping_file$BarcodeSequence,mapping_file$LinkerPrimerSequence)
mapping_file$BarcodeSequence <- BarcodeSequence_primer
Barcode_primer <- paste0(mapping_file$Barcode,mapping_file$Reverse)
mapping_file$Barcode <- Barcode_primer

# 4. Cut all the new barcodes to the same length -length of the longest barcode ####

# length of the shortest barcode
# BarcodeSequence column:
#nchar(mapping_file$BarcodeSequence)
BarcodeSequence_trimmed <- strtrim(mapping_file$BarcodeSequence,max_BarcodeSequence)
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
#In Terminal: validate_mapping_file.py -m mapping-miseq-run-test-2016.txt -o validate-mapping

# check if there are reverses of each barcodes of F and R in R1 and R2, and for the paired barcodes as well (reversed and not reversed)
# In R: reverse("ACAGCCACCCATCGAAAGAGCAACATCCTAA")

# 8. Make mapping file for idemp - barcode tab identifier:####
idemp_mapping_file <- mapping_file[,c(2,1)]
colnames(idemp_mapping_file) <- NULL
write.table(idemp_mapping_file, "R_projects/miseq-test/idemp-mapping-file.txt", sep = "\t", row.names = F, quote = F)

# 9. Extract barcodes - qiime:####
# In terminal: extract_barcodes.py -f rawdata/adapter-trimmed-new/Mona-run-test_S1_L001-trimmed_R1.fastq -r rawdata/adapter-trimmed-new/Mona-run-test_S1_L001-trimmed_R2.fastq -c barcode_paired_end -m mapping-miseq-run-test-2016.txt --attempt_read_reorientation --bc1_len 17 --bc2_len 14 -o processed_seqs

# 10. Demultiplex - idemp:####
#idemp -b idemp-mapping-file.txt -I1 processed_seqs/barcodes.fastq -R1 processed_seqs/reads1.fastq -R2 processed_seqs/reads2.fastq -m 1 -o idemp-demultiplexed





